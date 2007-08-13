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
public abstract class MaxLinkageDistanceCalculator implements SimilarityDistanceCalculator {
    /**
     * Returns the distance between an instance and the instances in a cluster.
     * The default implementation calculates maximum linkeage (max of the
     * distances between instances in the cluster and instanceIndex).
     * @param cluster Cluster array
     * @param clusterSize Number of the cluster array that contain instances.
     * Other elements must not be accessed.
     * @param instanceIndex Index of the instance that is compared to the
     * cluster.
     */

    public final double distance(final int[] cluster, final int clusterSize,
                                 final int instanceIndex) {
        double maxDistance = 0;

        for (int i = 0; i < clusterSize; ++i) {
            final int anInstance = cluster[i];

            maxDistance = Math.max(distance(anInstance, instanceIndex), maxDistance);
        }

        return maxDistance;
    }

    public final double getIgnoreDistance() {
        return Integer.MIN_VALUE; // Return the minimum integer value, so that max(min, a)=a;
    }
}
