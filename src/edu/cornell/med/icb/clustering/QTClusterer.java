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

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.mg4j.util.ProgressLogger;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Implements the QT Clustering algorithm. (QT stands for Quality
 * Threshold/diameter of the cluster)
 * See  http://en.wikipedia.org/wiki/Data_clustering#Types_of_clustering
 * and http://www.genome.org/cgi/content/full/9/11/1106 Heyer LJ et al 1999.
 * (just implement the interface
 * {@link edu.cornell.med.icb.clustering.SimilarityDistanceCalculator}).
 *
 */
public final class QTClusterer extends AbstractQTClusterer {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER = Logger.getLogger(QTClusterer.class);

    /**
     * The list of clusters.
     */
    private final IntArrayList[] clusters;

    /**
     *  A "temporary" cluster list used to find the maximum cardinality.
     */
    private final IntArrayList[] tmpClusters;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     * i.e., |G| where G is the set of instances
     */
    public QTClusterer(final int numberOfInstances) {
        super(numberOfInstances);
        clusters = new IntArrayList[instanceCount];
        tmpClusters = new IntArrayList[instanceCount];
        for (int i = 0; i < instanceCount; i++) {
            clusters[i] = new IntArrayList();     // NOPMD
            tmpClusters[i] = new IntArrayList();  // NOPMD
        }
    }

    /**
     * Groups instances into clusters. Returns the indices of the instances
     * that belong to a cluster as an int array in the list result.
     *
     * @param calculator The
     * {@link edu.cornell.med.icb.clustering.SimilarityDistanceCalculator}
     * that should be used when clustering
     * @param qualityThreshold The QT clustering algorithm quality threshold (d)
     * @return The list of clusters.
     */
    public List<int[]> cluster(
            final SimilarityDistanceCalculator calculator,
            final float qualityThreshold) {
        final ProgressLogger clusterProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "instances clustered");
        clusterProgressLogger.displayFreeMemory = true;
        clusterProgressLogger.expectedUpdates = instanceCount;
        clusterProgressLogger.start("Starting to cluster "
                + instanceCount + " instances");

        // reset cluster results
        clusterCount = 0;
        // instanceList is the set "G" to cluster
        final LinkedList<Integer> instanceList = new LinkedList<Integer>();
        for (int i = 0; i < instanceCount; i++) {
            clusters[i].clear();

            // set each node in the instance list to it's
            // original position in the source data array
            instanceList.add(i);
        }

        final ProgressLogger innerLoopProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "inner loop iterations");
        innerLoopProgressLogger.displayFreeMemory = false;

        final ProgressLogger outerLoopProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "outer loop iterations");
        outerLoopProgressLogger.displayFreeMemory = false;

        // loop over the instances until they have all been added to a cluster
        while (!instanceList.isEmpty()) {
            // cluster the remaining instances to find the maximum cardinality
            for (int i = 0; i < instanceCount; i++) {
                tmpClusters[i].clear();
            }

            outerLoopProgressLogger.expectedUpdates = instanceList.size();
            outerLoopProgressLogger.start("Entering outer loop for "
                    + instanceList.size() + " iterations");
            // foreach i in G (instance list)
            for (int i = 0; i < instanceList.size(); i++) {
                @SuppressWarnings("unchecked")
                final LinkedList<Integer> notClustered =
                        (LinkedList<Integer>) instanceList.clone();

                // add the first instance to the next cluster
                tmpClusters[i].add(notClustered.remove(i));

                innerLoopProgressLogger.expectedUpdates = notClustered.size();
                innerLoopProgressLogger.start("Entering inner loop for "
                        +  notClustered.size() + " iterations");

                // cluster the remaining instances to find the maximum
                // cardinality find instance j such that distance i,j minimum
                boolean done = false;
                while (!done && !notClustered.isEmpty()) {
                    // find the node that has minimum distance between the
                    // current cluster and the instances that have not yet
                    // been clustered
                    double minDistance = Double.MAX_VALUE;
                    int minDistanceInstanceIndex = 0;
                    int instanceIndex = 0;
                    for (final int instance : notClustered) {
                        final double newDistance =
                                calculator.distance(tmpClusters[i].elements(),
                                        tmpClusters[i].size(), instance);

                        if (newDistance <= minDistance) {
                            minDistance = newDistance;
                            minDistanceInstanceIndex = instanceIndex;
                        }
                        instanceIndex++;
                    }

                    // grow clusters until min distance between new instance
                    // and cluster reaches quality threshold
                    // if (diameter(Ai U {j}) > d)
                    if (minDistance > qualityThreshold) {
                        done = true;
                    } else {
                        tmpClusters[i].add(notClustered.remove(minDistanceInstanceIndex));
                    }
                    innerLoopProgressLogger.lightUpdate();
                }
                innerLoopProgressLogger.stop("Inner loop completed.");
                outerLoopProgressLogger.lightUpdate();
            }
            outerLoopProgressLogger.stop("Outer loop completed.");

            // identify cluster (set C) with maximum cardinality
            int maxCardinality = 0;
            int selectedClusterIndex = -1;
            for (int i = 0; i < instanceCount; i++) {
                final int size = tmpClusters[i].size();
                if (LOGGER.isTraceEnabled() && size > 0) {
                    LOGGER.trace("potential cluster " + i + ": "
                            + ArrayUtils.toString(tmpClusters[i]));
                }
                if (size > maxCardinality) {
                    maxCardinality = size;
                    selectedClusterIndex = i;
                }
            }

            final IntArrayList selectedCluster =
                    tmpClusters[selectedClusterIndex];

            if (LOGGER.isTraceEnabled()) {
                LOGGER.trace("adding " + selectedCluster.size()
                        + " instances to cluster " + clusterCount);
            }
            // and add that cluster to the final result
            clusters[clusterCount].addAll(selectedCluster);

            // remove instances in the cluster C so they are no longer considered
            instanceList.removeAll(selectedCluster);

            for (int i = 0; i < selectedCluster.size(); i++) {
                clusterProgressLogger.lightUpdate();
            }

            clusterCount++;
            // next iteration is over (G - C)
        }

        clusterProgressLogger.stop("Clustering completed.");
        return getClusters();
    }

    /**
     * Returns the list of clusters produced by clustering.
     *
     * @return A list of integer arrays, where each array represents a cluster
     *         and contains the index of the instance that belongs to a given
     *         cluster.
     */
    public List<int[]> getClusters() {
        final List<int[]> result = new ArrayList<int[]>(clusterCount);
        for (int i = 0; i < clusterCount; i++) {
            result.add(clusters[i].toIntArray());
        }
        return result;
    }
}
