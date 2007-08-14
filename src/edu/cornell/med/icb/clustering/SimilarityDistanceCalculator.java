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
 * Calculates a similarity measure between two instances available to the
 * clustering algorithm, or between one instance and a cluster.
 *
 * @author Fabien Campagne
 * Date: Oct 2, 2005
 * Time: 6:32:55 PM
 */
public interface SimilarityDistanceCalculator {
    /**
     * Returns the distance between two instances.
     *
     * @param instanceIndex Index of the first instance.
     * @param otherInstanceIndex  Index of the second instance.
     * @return distance measure between the two instances.
     */
    double distance(int instanceIndex, int otherInstanceIndex);

    /**
     * Returns the distance between an instance and the instances in a cluster.
     *
     * @param cluster Cluster array
     * @param clusterSize Number of the cluster array that contain instances.
     * Other elements must not be accessed.
     * @param instanceIndex Index of the instance that is compared to the
     * cluster.
     * @return the distance between an instance and the instances in a cluster.
     */
    double distance(int[] cluster, int clusterSize, int instanceIndex);

    /**
     * When some distances between instance pairs are missing/unknown the
     * ignoreDistance is returned. The clustering algorithm uses ignoreDistance
     * to recognize cases when the distance is unknown.
     *
     * @return The distance value that the linkage method will ignore.
     */
    double getIgnoreDistance();
}
