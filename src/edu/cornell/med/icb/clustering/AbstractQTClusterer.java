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

import it.unimi.dsi.mg4j.util.ProgressLogger;

/**
 * Base class for clustering implementations.
 */
public abstract class AbstractQTClusterer implements Clusterer {
    /**
     * The time interval for a new log in milliseconds.
     *
     * @see it.unimi.dsi.mg4j.util.ProgressLogger#DEFAULT_LOG_INTERVAL
     */
    protected long logInterval = ProgressLogger.DEFAULT_LOG_INTERVAL;

    /**
     * The total number of instances to cluster.
     */
    protected final int instanceCount;

    /**
     * The total number of clusters created by clustering the instances.
     */
    protected int clusterCount;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     */
    public AbstractQTClusterer(final int numberOfInstances) {
        super();
        if (numberOfInstances < 0) {
            throw new IllegalArgumentException("Number of instances ("
                    + numberOfInstances + ") must not be negative");
        }
        clusterCount = numberOfInstances;
        instanceCount = numberOfInstances;
    }

    /**
     * Get the the progress logging interval.
     *
     * @return the logging interval in milliseconds.
     */
    public final long getLogInterval() {
        return logInterval;
    }

    /**
     * Set the the progress logging interval.
     *
     * @param interval the logging interval in milliseconds.
     */
    public final void setLogInterval(final long interval) {
        this.logInterval = interval;
    }
}
