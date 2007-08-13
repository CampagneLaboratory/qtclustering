package edu.cornell.med.icb.clustering;

/**
 * Calculates a similarity measure between two instances available to the
 * clustering algorithm, or between one instance and a cluster.
 *
 * @author Fabien Campagne
 *
 * Date: Oct 2, 2005
 * Time: 6:32:55 PM
 */
public interface SimilarityDistanceCalculator {
    /**
     * Returns the distance between two instances. *
     *
     * @param instanceIndex Index of the first instance.
     * @param otherInstanceIndex  Index of the second instance.
     * @return distance measure between the two instances.
     */
    double distance(final int instanceIndex, final int otherInstanceIndex);

    /**
     * Returns the distance between an instance and the instances in a cluster.
     */
    double distance(int[] cluster, int clusterSize, int instanceIndex);

    /**
     * When some distances between instance pairs are missing/unknown the ignoreDistance is returned.
     * The clustering algorithm uses ignoreDistance to recognize cases when the distance is unknown.
     *
     * @return The distance value that the linkage method will ignore.
     */
    double getIgnoreDistance();
}
