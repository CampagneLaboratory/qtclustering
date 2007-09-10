/*
 * Copyright (C) 2005-2007 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import org.apache.log4j.Logger;
import org.apache.commons.lang.ArrayUtils;

import java.util.List;
import java.util.ArrayList;

/**
 * Implements the QT Clustering algorithm. (QT stands for Quality
 * Threshold/diameter of the cluster)
 * See  http://en.wikipedia.org/wiki/Data_clustering#Types_of_clustering
 * and http://www.genome.org/cgi/content/full/9/11/1106 Heyer LJ et al 1999.
 * This implementation should be fairly efficient. Memory needed by the
 * algorithm is allocated only at the beginning of the clustering process
 * (not as clustering proceeds). The implementation makes it easy to plugin
 * a new distance similarity measure or a new linkage method
 * (just implement the interface
 * {@link edu.cornell.med.icb.clustering.SimilarityDistanceCalculator}).
 *
 * @author Fabien Campagne
 *         Date: Oct 2, 2005
 *         Time: 6:23:51 PM
 */
public final class RecursiveQTClusterer extends AbstractQTClusterer {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER =
            Logger.getLogger(RecursiveQTClusterer.class);
    private final int[][] clusters;
    private final int[] clusterSizes;
    private final Int2BooleanMap jVisited;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     */
    public RecursiveQTClusterer(final int numberOfInstances) {
        super(numberOfInstances);

        // Data elements are available only if k<clusterSizes[l]
        clusters = new int[clusterCount][numberOfInstances];
        // number of instances in each cluster.
        clusterSizes = new int[clusterCount];
        resetTmpClusters();
        jVisited = new Int2BooleanAVLTreeMap();

    }

    /**
     * Groups instances into clusters. Returns the indices of the instances
     * that belong to a cluster as an int array in the list result.
     *
     * @param calculator       The distance calculator to
     * @param qualityThreshold The QT clustering algorithm quality threshold.
     * @return The list of clusters.
     */
    public List<int[]> cluster(
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

        cluster(result, calculator, qualityThreshold, ignoreList,
                instanceCount, clusterProgressLogger);

        clusterProgressLogger.done();
        return result;
    }

    /**
     * Performs the actual clustering.
     *
     * @param result A list that should be used to store the results.
     * @param calculator The
     * {@link edu.cornell.med.icb.clustering.SimilarityDistanceCalculator}
     * that should be used when clustering
     * @param qualityThreshold The QT clustering algorithm quality threshold (d)
     * @param ignoreList
     * @param instances
     * @param progressLogger A {@link it.unimi.dsi.mg4j.util.ProgressLogger}
     * that should used to update clustering progress.
     */
    private void cluster(final List<int[]> result,
                         final SimilarityDistanceCalculator calculator,
                         final float qualityThreshold,
                         final Int2BooleanMap ignoreList,
                         final int instances,
                         final ProgressLogger progressLogger) {
        resetTmpClusters();
        if (instances <= 1) { // one instance -> one cluster
            if (instances > 0) {
                final int[] singletonCluster = new int[instances];

                // find the remaining instance:
                for (int i = 0; i < instanceCount; ++i) {
                    if (!ignoreList.get(i)) {
                        singletonCluster[0] = i;
                    }
                }

                result.add(singletonCluster);
                progressLogger.update();
            }
            return;
        }

        final ProgressLogger innerLoopProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "inner loop iterations");
        innerLoopProgressLogger.displayFreeMemory = false;

        final ProgressLogger outerLoopProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "outer loop iterations");
        outerLoopProgressLogger.displayFreeMemory = false;
        outerLoopProgressLogger.expectedUpdates = instanceCount;
        outerLoopProgressLogger.start();

        for (int i = 0; i < instanceCount; ++i) { // i : cluster index
            // ignore instance i if part of previously selected clusters.
            if (!clustersCannotOverlap || !ignoreList.get(i)) {
                boolean done = false;
                jVisited.clear();
                addToCluster(i, i);
                jVisited.put(i, true);
                while (!done) {
                    double distance_i_j = Double.MAX_VALUE;
                    int minDistanceInstanceIndex = -1;
                    innerLoopProgressLogger.expectedUpdates = instanceCount;
                    innerLoopProgressLogger.start();

                    // find instance j such that distance i,j minimum
                    for (int j = 0; j < instanceCount; ++j) {
                        if (!ignoreList.get(j) && !jVisited.get(j)) {
                            final double newDistance =
                                    calculator.distance(clusters[i],
                                            clusterSizes[i], j);

                            if (newDistance < distance_i_j) {
                                distance_i_j = newDistance;
                                minDistanceInstanceIndex = j;
                            }
                        }
                        innerLoopProgressLogger.update();
                    }

                    // grow clusters until min distance between new instance
                    // and cluster reaches quality threshold:
                    if (distance_i_j > qualityThreshold) {
                        done = true;
                    } else {
                        if (clustersCannotOverlap
                                && ignoreList.get(minDistanceInstanceIndex)) {
                            done = true;
                        } else {
                            final boolean added =
                                    addToCluster(minDistanceInstanceIndex, i);
                            jVisited.put(minDistanceInstanceIndex, true);
                            if (!added) {
                                done = true;
                                LOGGER.info(String.format("Could not add instance minDistanceInstanceIndex=%d to cluster %d, distance was %f\n", minDistanceInstanceIndex, i, distance_i_j));

                            } else {
                                innerLoopProgressLogger.expectedUpdates = instanceCount;
                            }
                        }
                    }
                }
                innerLoopProgressLogger.done();
            }
            outerLoopProgressLogger.update();
        }
        outerLoopProgressLogger.done();

        // identify cluster with maximum cardinality
        int maxCardinality = 0;
        int selectedClusterIndex = -1;
        for (int l = 0; l < clusterSizes.length; ++l) {
            final int size = clusterSizes[l];
            if (LOGGER.isTraceEnabled() && size > 0) {
                LOGGER.trace("potential cluster " + l + ": "
                        + ArrayUtils.toString(getClusters().get(l)));
            }

            if (size > maxCardinality) {
                maxCardinality = clusterSizes[l];
                selectedClusterIndex = l;
            }
        }
        final int[] selectedCluster =
                getClusters().get(selectedClusterIndex).clone();
        result.add(selectedCluster);

        int instancesLeft = instances;
        for (final int ignoreInstance : selectedCluster) {
            // mark instances of the selected cluster so that they are
            // ignored in subsequent passes.
            ignoreList.put(ignoreInstance, true);
            progressLogger.update();
            --instancesLeft;
        }

        // recurse.
        cluster(result, calculator, qualityThreshold, ignoreList,
                instancesLeft, progressLogger);
    }

    /**
     * Initialize instance index in clusters to value that is not supported.
     */
    private void resetTmpClusters() {
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
     * Add an instance to a cluster.
     *
     * @param instanceIndex Index of the instance to add to the cluster.
     * @param clusterIndex  Index of the cluster where to add the instance.
     * @return true if instance appended to cluster, false otherwise
     */
    private boolean addToCluster(final int instanceIndex,
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
}