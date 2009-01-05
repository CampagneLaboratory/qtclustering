/*
 * Copyright (C) 2005-2009 Institute for Computational Biomedicine,
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
