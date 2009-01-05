/*
 * Copyright (C) 2007-2009 Institute for Computational Biomedicine,
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

import java.util.List;

/**
 * Common interface definition for clustering algorithms.  Clustering is the
 * classification of objects into different groups, or more precisely,
 * the partitioning of a data set into subsets (clusters), so that the
 * data in each subset (ideally) share some common trait - often proximity
 * according to some defined distance measure.
 */
public interface Clusterer {
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
    List<int[]> cluster(final SimilarityDistanceCalculator calculator,
                        final double qualityThreshold);

    /**
     * Returns the list of clusters produced by clustering.
     *
     * @return A list of integer arrays, where each array represents a cluster
     * and contains the index of the instance that belongs to a given cluster.
     */
    List<int[]> getClusters();
}
