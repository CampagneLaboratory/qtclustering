/*
 * Copyright (C) 2007 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.clustering;

import it.unimi.dsi.fastutil.ints.Int2BooleanMap;
import it.unimi.dsi.fastutil.ints.Int2BooleanAVLTreeMap;
import it.unimi.dsi.mg4j.util.ProgressLogger;

import java.util.List;
import java.util.ArrayList;

import org.apache.log4j.Logger;

/**
 * Base class for clustering implementationsx.
 */
public abstract class AbstractQTClusterer implements Clusterer {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER =
            Logger.getLogger(AbstractQTClusterer.class);

    /**
     * Indicate that clusters cannot overlap.
     */
    protected boolean clustersCannotOverlap = true;

    protected final int clusterCount;
    protected final int[][] clusters;
    protected final int[] clusterSizes;
    protected final int instanceCount;

    protected final Int2BooleanMap jVisited;
    /**
     * The time interval for a new log in milliseconds.
     *
     * @see ProgressLogger#DEFAULT_LOG_INTERVAL
     */
    protected long logInterval = ProgressLogger.DEFAULT_LOG_INTERVAL;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     */
    public AbstractQTClusterer(final int numberOfInstances) {
        super();
        clusterCount = numberOfInstances;
        // Data elements are available only if k<clusterSizes[l]
        clusters = new int[clusterCount][numberOfInstances];
        // number of instances in each cluster.
        clusterSizes = new int[clusterCount];
        instanceCount = numberOfInstances;

        resetTmpClusters();
        jVisited = new Int2BooleanAVLTreeMap();
    }

    /**
     * Initialize instance index in clusters to value that is not supported.
     */
    protected final void resetTmpClusters() {
        for (final int[] cluster : clusters) {
            for (int i = 0; i < cluster.length; ++i) {
                cluster[i] = Integer.MAX_VALUE;
            }
        }
        for (int i = 0; i < clusterSizes.length; ++i) {
            clusterSizes[i] = 0;
        }
    }

    /**
     * Get the the progress logging interval.
     *
     * @return the logging interval in milliseconds.
     */
    public final long getLogInterval() {
        return logInterval;
    }

    /**
     * Set the the progress logging interval.
     *
     * @param interval the logging interval in milliseconds.
     */
    public final void setLogInterval(final long interval) {
        this.logInterval = interval;
    }

    /**
     * If clustersCannotOverlap is true, then clustering will produce clusters
     * that do not overlap. If clustersCannotOverlap is false, overlapping is
     * allowed, and some instances will be part of several clusters.
     *
     * @return whether or not clusters can overlap
     */
    public final boolean isClustersCannotOverlap() {
        return clustersCannotOverlap;
    }

    /**
     * Indicate that clusters cannot overlap. If clustersCannotOverlap is true,
     * then clustering will produce clusters that do not overlap. If
     * clustersCannotOverlap is false, overlapping is allowed, and some
     * instances will be part of several clusters.
     *
     * @param cannotOverlap Indicates whether or not clusters can overlap
     */
    public final void setClustersCannotOverlap(final boolean cannotOverlap) {
        this.clustersCannotOverlap = cannotOverlap;
    }

    /**
     * Add an instance to a cluster.
     *
     * @param instanceIndex Index of the instance to add to the cluster.
     * @param clusterIndex  Index of the cluster where to add the instance.
     * @return true if instance appended to cluster, false otherwise
     */
    public final boolean addToCluster(final int instanceIndex,
                                      final int clusterIndex) {
        if (instanceIndex == Integer.MAX_VALUE) {
            throw new IllegalArgumentException("instance Index =="
                + Integer.MAX_VALUE + " is not supported");
        }

        if (clusterIndex >= clusters.length) {
            throw new IllegalArgumentException("Cluster index ("
            + clusterIndex + ") must be < cluser length ("
            + clusters.length + ")");
        }

        // return immediately if instance already in cluster;
        for (int i = 0; i < clusterSizes[clusterIndex]; ++i) {
            if (clusters[clusterIndex][i] == instanceIndex) {
                return false;
            }
        }

        final int lastInstanceIndex = clusterSizes[clusterIndex];
        assert lastInstanceIndex < clusters[clusterIndex].length;

        clusters[clusterIndex][lastInstanceIndex] = instanceIndex;
        clusterSizes[clusterIndex] += 1;
        return true;
    }

    /**
     * Performs the actual clustering.
     *
     * @param result           A list that should be used to store the results.
     * @param calculator       The {@link SimilarityDistanceCalculator} that
     *                         should be used when clustering
     * @param qualityThreshold
     * @param ignoreList
     * @param instances
     * @param progressLogger   A {@link ProgressLogger} that should used to
     *                         update clustering progress.
     */
    protected abstract void cluster(List<int[]> result,
                         SimilarityDistanceCalculator calculator,
                         float qualityThreshold,
                         Int2BooleanMap ignoreList,
                         int instances,
                         ProgressLogger progressLogger) throws Exception;

    /**
     * Returns the list of clusters produced by clustering.
     *
     * @return A list of integer arrays, where each array represents a cluster
     *         and contains the index of the instance that belongs to a given cluster.
     */
    public final List<int[]> getClusters() {
        final List<int[]> result = new ArrayList<int[]>();
        final int resultClusterCount = clusterCount;

        for (int l = 0; l < resultClusterCount; ++l) {
            final int[] trimmedCluster = new int[clusterSizes[l]]; // NOPMD
            System.arraycopy(clusters[l], 0, trimmedCluster, 0, clusterSizes[l]);
            result.add(trimmedCluster);
        }
        return result;
    }

    /**
     * Groups instances into clusters. Returns the indices of the instances
     * that belong to a cluster as an int array in the list result.
     *
     * @param calculator       The distance calculator to
     * @param qualityThreshold The QT clustering algorithm quality threshold.
     * @return The list of clusters.
     */
    public final List<int[]> cluster(
            final SimilarityDistanceCalculator calculator,
            final float qualityThreshold) {
        final ProgressLogger clusterProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "instances clustered");
        clusterProgressLogger.displayFreeMemory = true;
        clusterProgressLogger.expectedUpdates = instanceCount;
        clusterProgressLogger.start("Starting to cluster");

        final List<int[]> result = new ArrayList<int[]>();
        // set of instances to ignore.
        // Map returns true if instance must be ignored.
        final Int2BooleanAVLTreeMap ignoreList = new Int2BooleanAVLTreeMap();

        try {
            cluster(result, calculator, qualityThreshold, ignoreList,
                    instanceCount, clusterProgressLogger);
        } catch (Exception e) {
            if (e instanceof RuntimeException) {
                throw (RuntimeException) e;
            } else {
                throw new RuntimeException(e);
            }
        }
        clusterProgressLogger.done();
        return result;
    }
}
