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

import it.unimi.dsi.mg4j.util.ProgressLogger;
import org.apache.log4j.Logger;

/**
 * Maximum distance linkage calculator that uses a local "cache" of distances between each
 * pair of instances being considered for clustering. Calculate the distance of a point to
 * a cluster as the maximum distance between the point and each point of the cluster.
 */
public abstract class CachingMaxLinkageDistanceCalculator
        extends MaxLinkageDistanceCalculator {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER =
            Logger.getLogger(CachingMaxLinkageDistanceCalculator.class);

    /**
     * Representation of the specified floating-point value according to the
     * IEEE 754 floating-point "double format" bit layout.
     */
    private static final long negativeZeroDoubleBits =
            Double.doubleToLongBits(-0.0d);

    /**
     * A matrix of the distances between every pair of instances.
     */
    private double[][] distanceCache;

    /**
     * Create a new {@link SimilarityDistanceCalculator}.
     */
    public CachingMaxLinkageDistanceCalculator() {
        super();
    }

    /**
     * Initializes the calculator with with the number of instances that should be considered
     * when clustering.
     *
     * @param numberOfInstances The largest number of instances that may be considered
     * for clustering.
     */
    @Override
    public void initialize(final int numberOfInstances) {
        super.initialize(numberOfInstances);
        distanceCache = new double[numberOfInstances][numberOfInstances];

        final ProgressLogger initializationProgressLogger =
                new ProgressLogger(LOGGER, ProgressLogger.DEFAULT_LOG_INTERVAL,
                        "instances initialized");
        initializationProgressLogger.displayFreeMemory = true;
        initializationProgressLogger.expectedUpdates = numberOfInstances;
        initializationProgressLogger.start("Initializing "
                + numberOfInstances + " instances.");


        // We don't assume the value distance(i,j) is the same as distance(j,i)
        // although we really could if we need to save a little bit of memory
        for (int i = 0; i < numberOfInstances; i++) {
            for (int j = 0; j < numberOfInstances; j++) {
                distanceCache[i][j] = distance(i, j);
            }
            initializationProgressLogger.update();
        }

        initializationProgressLogger.done();
    }

    /**
     * Returns the distance between an instance and the instances in a cluster.
     * The default implementation calculates maximum linkage (max of the
     * distances between instances in the cluster and instanceIndex).
     * @param cluster A list of indices that represent instances in a cluster
     * @param clusterSize The number of elements in the cluster
     * @param instanceIndex Index of the instance that is compared to the cluster.
     * @return the distance between an instance and the instances in a cluster.
     */
    @Override
    public final double distance(final int[] cluster, final int clusterSize,
                                 final int instanceIndex) {
        double maxDistance = ignoreDistance;

        for (int i = 0; i < clusterSize; i++) {
            final double a = distanceCache[cluster[i]][instanceIndex];
            final double b = maxDistance;

            // This code is inlined from java.lang.Math.max(a, b)
            if (a != a) {             // a is NaN
                maxDistance = a;
            } else if (a == 0.0d && b == 0.0d
                    && Double.doubleToLongBits(a) == negativeZeroDoubleBits) {
                maxDistance = b;
            } else if (a >= b) {
                maxDistance = a;
            } else {
                maxDistance = b;
            }
        }

        return maxDistance;
    }
}