/*
 * Copyright © 2005-2007 Institute for Computational Biomedicine,
 *                       Weill Medical College of Cornell University
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

import it.unimi.dsi.fastutil.ints.Int2BooleanAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2BooleanMap;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements the QT_Cluster algorithm. (QT stands for Quality
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
 * Date: Oct 2, 2005
 * Time: 6:23:51 PM
 */
public final class QTClusterer {
    private final int clusterCount;
    private final int[][] clusters;
    private final int[] clusterSizes;
    private final int instanceCount;
    private boolean clustersCannotOverlap;

    private final Int2BooleanMap jVisited;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param instanceCount The number of instances to cluster.
     */
    public QTClusterer(final int instanceCount) {
        clusterCount = instanceCount;
        // Data elements are available only if k<clusterSizes[l]
        clusters = new int[clusterCount][instanceCount];
        // number of instances in each cluster.
        clusterSizes = new int[clusterCount];
        this.instanceCount = instanceCount;

        resetTmpClusters();
        jVisited = new Int2BooleanAVLTreeMap();

        setClustersCannotOverlap(true);
    }

    private void resetTmpClusters() {
        // initialize instance index in clusters to value that is not supported.
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
     * Indicate that clusters cannot overlap. If clustersCannotOverlap is true,
     * then clustering will produce clusters that do not overlap. If
     * clustersCannotOverlap is false, overlapping is allowed, and some
     * instances will be part of several clusters.
     *
     * @param clustersCannotOverlap
     */
    public void setClustersCannotOverlap(final boolean clustersCannotOverlap) {
        this.clustersCannotOverlap = clustersCannotOverlap;
    }

    public List<int[]> cluster(final SimilarityDistanceCalculator calculator,
                               final float qualityThreshold) {
        final List<int[]> result = new ArrayList<int[]>();
        final Int2BooleanAVLTreeMap ignoreList = new Int2BooleanAVLTreeMap();      // set of instances to ignore. Map returns
        // true if instance must be ignored.
        cluster(result, calculator, qualityThreshold, ignoreList, instanceCount);
        return result;
    }

    private void cluster(final List<int[]> result,
                         final SimilarityDistanceCalculator calculator,
                         final float qualityThreshold,
                         final Int2BooleanMap ignoreList, int instancesLeft) {
        resetTmpClusters();
        if (instancesLeft <= 1) { // one instance -> one cluster
            if (instancesLeft <= 0) {
                return;
            }
            final int[] singletonCluster = new int[instancesLeft];

            // find the remaining instance:
            for (int i = 0; i < instanceCount; ++i) {
                if (!ignoreList.get(i)) {
                    singletonCluster[0] = i;
                }
            }

            result.add(singletonCluster);
            return;
        }

        for (int i = 0; i < instanceCount; ++i) { // i : cluster index
            if (clustersCannotOverlap && ignoreList.get(i)) {
                continue; // ignore instance i if part of previously selected clusters.
            }
            boolean done = false;
            addToCluster(i, i);
            jVisited.clear();
            while (!done && instancesLeft > 1) {
                double distance_i_j = Double.MAX_VALUE;
                int minDistanceInstanceIndex = -1;

                for (int j = 0; j < instanceCount; ++j) { // find instance j such that distance i,j minimum
                    if (ignoreList.get(j)) {
                        continue;
                    }
                    if (i != j) {
                        if (jVisited.get(j)) {
                            continue;
                        }

                        final double newDistance = calculator.distance(clusters[i],
                                clusterSizes[i], j);

                        if (newDistance < distance_i_j) {
                            distance_i_j = newDistance;
                            minDistanceInstanceIndex = j;
                            jVisited.put(j, true);
                        }
                    }
                }
                // grow clusters until min distance between new instance and cluster reaches quality threshold:
                if (distance_i_j > qualityThreshold) {
                    done = true;
                } else {
                    if (clustersCannotOverlap && ignoreList.get(minDistanceInstanceIndex)) {
                        done = true;
                    } else {
                        final boolean added;
                        added = addToCluster(minDistanceInstanceIndex, i);
                        if (!added && jVisited.get(minDistanceInstanceIndex)) {
                            done = true;
                        }
                    }
                }
            }
        }
        // identify cluster with maximum cardinality:
        int maxCardinality = 0;
        int selectedClusterIndex = -1;
        for (int l = 0; l < clusterSizes.length; ++l) {
            if (clusterSizes[l] > maxCardinality) {
                maxCardinality = clusterSizes[l];
                selectedClusterIndex = l;
            }
        }
        final int[] selectedCluster = getClusters().get(selectedClusterIndex).clone();
        result.add(selectedCluster);
        for (final int ignoreInstance : selectedCluster) {
            // mark instances of the selected cluster so that they are ignored in subsequent passes.
            ignoreList.put(ignoreInstance, true);
            --instancesLeft;
        }

        // recurse.
        cluster(result, calculator, qualityThreshold, ignoreList, instancesLeft);
    }

    /**
     * Add an instance to a cluster.
     *
     * @param instanceIndex
     * @param clusterIndex
     * @return True if instance appended to cluster, false otherwise
     */
    public boolean addToCluster(final int instanceIndex, final int clusterIndex) {
        assert instanceIndex != Integer.MAX_VALUE : "instance Index ==" + Integer.MAX_VALUE + " is not supported";
        // return immediately if instance already in cluster;
        for (int i = 0; i < clusterSizes[clusterIndex]; ++i) {
            if (clusters[clusterIndex][i] == instanceIndex) {
                return false;
            }
        }
        final int lastInstanceIndex = clusterSizes[clusterIndex];
        assert clusterIndex < clusters.length;
        assert lastInstanceIndex < clusters[clusterIndex].length;
        clusters[clusterIndex][lastInstanceIndex] = instanceIndex;
        clusterSizes[clusterIndex] += 1;
        return true;
    }

    /**
     * Returns the list of clusters produced by clustering.
     *
     * @return A list of integer arrays, where each array represents a cluster and contains the index of the
     *         instance that belongs to a given cluster.
     */
    public List<int[]> getClusters() {
        final List<int[]> result = new ArrayList<int[]>();
        final int resultClusterCount = clusterCount;

        for (int l = 0; l < resultClusterCount; ++l) {
            final int[] trimmedCluster = new int[clusterSizes[l]];
            System.arraycopy(clusters[l], 0, trimmedCluster, 0, clusterSizes[l]);
            result.add(trimmedCluster);
        }
        return result;
    }
}
