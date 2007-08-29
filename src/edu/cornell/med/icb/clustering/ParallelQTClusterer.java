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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedInteger;
import it.unimi.dsi.fastutil.ints.Int2BooleanMap;
import it.unimi.dsi.fastutil.ints.Int2BooleanAVLTreeMap;
import it.unimi.dsi.mg4j.util.ProgressLogger;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.ArrayList;

/**
 * Implements the QT Clustering algorithm. (QT stands for Quality
 * Threshold/diameter of the cluster)
 * See  http://en.wikipedia.org/wiki/Data_clustering#Types_of_clustering
 * and http://www.genome.org/cgi/content/full/9/11/1106 Heyer LJ et al 1999.
 * This implementation This implementation will run the clustering main loop
 * in multiple parallel threads if possible.  NOTE: The class itself is NOT
 * thread-safe.  Mulitple cluster calcuations cannot be executed with a single
 * instance of this class.  The implementation makes it
 * easy to plugin a new distance similarity measure or a new linkage method
 * (just implement the interface {@link SimilarityDistanceCalculator}).
 */
public final class ParallelQTClusterer extends AbstractQTClusterer {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER =
            Logger.getLogger(ParallelQTClusterer.class);

    /**
     * Number of threads to use when clustering.
     * @see ParallelTeam#getDefaultThreadCount()
     */
    private int threadCount;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     */
    public ParallelQTClusterer(final int numberOfInstances) {
        this(numberOfInstances, ParallelTeam.getDefaultThreadCount());
    }

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     * @param numberOfThreads Number of threads to use when clustering.
     */
    public ParallelQTClusterer(final int numberOfInstances,
                               final int numberOfThreads) {
        super(numberOfInstances);
        setThreadCount(numberOfThreads);
    }

    /**
     * Determine the number of threads to use when clustering.
     * @return Number of threads.
     */
    public int getThreadCount() {
        return threadCount;
    }

    /**
     * Specify the number of threads to use when clustering.
     * @param numberOfThreads Number of threads to use when clustering.
     * @throws IllegalArgumentException if the thread count is less than 1
     */
    public void setThreadCount(final int numberOfThreads) {
        if (numberOfThreads < 1) {
            throw new IllegalArgumentException("Thread count must be > 0");
        }
        threadCount = numberOfThreads;
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

        // don't create more threads than there are instances
        final int numberOfThreads =
                Math.max(1, Math.min(threadCount, instanceCount));
        LOGGER.debug("starting " + numberOfThreads + " threads");
        final ParallelTeam team = new ParallelTeam(numberOfThreads);

        try {
            cluster(team, result, calculator, qualityThreshold, ignoreList,
                    instanceCount, clusterProgressLogger);
        } catch (Exception e) {
            if (e instanceof RuntimeException) {
                throw (RuntimeException) e;
            } else {
                throw new ClusteringException(e);
            }
        }

        clusterProgressLogger.done();
        return result;
    }

    /**
     * Performs the actual clustering.
     *
     * @param team The team of threads for executing clustering
     * @param result A list that should be used to store the results.
     * @param calculator The {@link SimilarityDistanceCalculator}
     * that should be used when clustering
     * @param qualityThreshold
     * @param ignoreList
     * @param instances
     * @param progressLogger A {@link ProgressLogger}
     * that should used to update clustering progress.
     */
    protected void cluster(final ParallelTeam team,
                           final List<int[]> result,
                           final SimilarityDistanceCalculator calculator,
                           final float qualityThreshold,
                           final Int2BooleanMap ignoreList,
                           final int instances,
                           final ProgressLogger progressLogger) throws Exception {
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

        final ProgressLogger loopProgressLogger =
                new ProgressLogger(LOGGER, logInterval, "iterations");
        loopProgressLogger.displayFreeMemory = true;

        for (int clusterIndex = 0; clusterIndex < instanceCount; ++clusterIndex) { // i : cluster index
            final int i = clusterIndex;
            // ignore instance i if part of previously selected clusters.
            if (!clustersCannotOverlap || !ignoreList.get(i)) {
                boolean done = false;
                addToCluster(i, i);
                jVisited.clear();
                while (!done && instances > 1) {
                    final SharedDouble distance_i_j = new SharedDouble(Double.MAX_VALUE);
                    final SharedInteger minDistanceInstanceIndex = new SharedInteger(-1);
                    loopProgressLogger.expectedUpdates = instanceCount;
                    loopProgressLogger.start();

                    // find instance j such that distance i,j minimum
                    team.execute(new ParallelRegion() {
                        @Override
                        public void run() throws Exception {
                            execute(0, instanceCount - 1, new IntegerForLoop() {
                                @Override
                                public void run(final int first, final int last) throws Exception {
                                    for (int j = first; j <= last; ++j) {
                                        if (!ignoreList.get(j)) {
                                            if (i != j) {
                                                if (!jVisited.get(j)) {
                                                    final double newDistance =
                                                            calculator.distance(clusters[i],
                                                                    clusterSizes[i], j);
                                                    final int instance = j;
                                                    critical(new ParallelSection() {
                                                        @Override
                                                        public void run() throws Exception {
                                                            if (newDistance < distance_i_j.get()) {
                                                                distance_i_j.set(newDistance);
                                                                minDistanceInstanceIndex.set(instance);
                                                                jVisited.put(instance, true);
                                                            }
                                                        }
                                                    });

                                                }
                                            }
                                        }
                                        loopProgressLogger.update();
                                    }
                                }
                            });
                        }
                    });


                    // grow clusters until min distance between new instance
                    // and cluster reaches quality threshold:
                    if (distance_i_j.get() > qualityThreshold) {
                        done = true;
                    } else {
                        if (clustersCannotOverlap
                                && ignoreList.get(minDistanceInstanceIndex.get())) {
                            done = true;
                        } else {
                            final boolean added =
                                    addToCluster(minDistanceInstanceIndex.get(), i);
                            if (!added && jVisited.get(minDistanceInstanceIndex.get())) {
                                done = true;
                                LOGGER.info(String.format("Could not add instance minDistanceInstanceIndex=%d to cluster %d, distance was %f\n", minDistanceInstanceIndex.get(), i, distance_i_j.get()));

                            } else {
                                loopProgressLogger.expectedUpdates = instanceCount;
                            }
                        }
                    }
                }
            }

            loopProgressLogger.update();
        }

        // identify cluster with maximum cardinality
        int maxCardinality = 0;
        int selectedClusterIndex = -1;
        for (int l = 0; l < clusterSizes.length; ++l) {
            if (clusterSizes[l] > maxCardinality) {
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
        cluster(team, result, calculator, qualityThreshold, ignoreList,
                instancesLeft, progressLogger);
    }
}

