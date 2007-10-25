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

/**
 * Maximum distance linkage calculator. Calculate the distance of a point to
 * a cluster as the maximum distance between the point and each point of the
 * cluster.
 *
 * @author Fabien Campagne
 * Date: Oct 4, 2005
 * Time: 5:35:52 PM
 */
public abstract class MaxLinkageDistanceCalculator
        implements SimilarityDistanceCalculator {
    /**
     * Representation of the specified floating-point value according to the
     * IEEE 754 floating-point "double format" bit layout.
     */
    private static final long negativeZeroDoubleBits =
            Double.doubleToLongBits(-0.0d);

    /**
     * Default value to use when pairs are missing or unknown.
     */
    protected final double ignoreDistance = Double.NEGATIVE_INFINITY;

    /**
     * Create a new {@link edu.cornell.med.icb.clustering.SimilarityDistanceCalculator}.
     */
    public MaxLinkageDistanceCalculator() {
        super();
    }

    /**
     * Initializes the calculator with with the number of instances that should be considered
     * when clustering.
     *
     * @param numberOfInstances The largest number of instances that may be considered
     * for clustering.
     */
    public void initialize(final int numberOfInstances) {
        // this implementation does not require any specific initialization
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
    public double distance(final int[] cluster, final int clusterSize,
                           final int instanceIndex) {
        double maxDistance = ignoreDistance;

        for (int i = 0; i < clusterSize; i++) {
            final double a = distance(cluster[i], instanceIndex);
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

    /**
     * When some distances between instance pairs are missing/unknown the
     * ignoreDistance is returned. The clustering algorithm uses ignoreDistance
     * to recognize cases when the distance is unknown.
     *
     * @return The minimum value, so that max(min, a) = a;
     */
    public final double getIgnoreDistance() {
        return ignoreDistance;
    }
}
