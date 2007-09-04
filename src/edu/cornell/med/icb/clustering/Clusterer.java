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

import java.util.List;

public interface Clusterer {
    /**
     * Groups instances into clusters. Returns the indices of the instances
     * that belong to a cluster as an int array in the list result.
     *
     * @param calculator The distance calculator to
     * @param qualityThreshold The clustering algorithm quality threshold.
     * @return The list of clusters.
     */
    List<int[]> cluster(final SimilarityDistanceCalculator calculator,
                        final float qualityThreshold);

    /**
     * Returns the list of clusters produced by clustering.
     *
     * @return A list of integer arrays, where each array represents a cluster
     * and contains the index of the instance that belongs to a given cluster.
     */
    List<int[]> getClusters();

    /**
     * Indicate that clusters cannot overlap. If clustersCannotOverlap is true,
     * then clustering will produce clusters that do not overlap. If
     * clustersCannotOverlap is false, overlapping is allowed, and some
     * instances will be part of several clusters.
     *
     * @param cannotOverlap Indicates whether or not clusters can overlap
     */
    void setClustersCannotOverlap(final boolean cannotOverlap);

    /**
     * If clustersCannotOverlap is true, then clustering will produce
     * clusters that do not overlap. If clustersCannotOverlap is false,
     * overlapping is allowed, and some instances will be part of
     * several clusters.
     *
     * @return whether or not clusters can overlap
     */
    boolean isClustersCannotOverlap();
}
