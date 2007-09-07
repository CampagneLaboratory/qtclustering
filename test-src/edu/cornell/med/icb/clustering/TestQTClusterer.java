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

import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * User: Fabien Campagne
 * Date: Oct 2, 2005
 * Time: 6:59:35 PM
 * To change this template use File | Settings | File Templates.
 */
public final class TestQTClusterer {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER =
            Logger.getLogger(TestQTClusterer.class);

    /**
     * Default delta to use when comparing floating point values.
     */
    private static final double DELTA = 0.00001;

    /**
     * This test validates that each instance will be placed into it's own
     * cluster when there is no overlap between them.
     */
    @Test
    public void oneInstancePerCluster() {
        // put one instance in each cluster, total two instances
        final QTClusterer clusterer = new QTClusterer(2);
        final SimilarityDistanceCalculator distanceCalculator =
                new MaxLinkageDistanceCalculator() {
            @Override
            public double distance(final int instanceIndex,
                                   final int otherInstanceIndex) {
                if (instanceIndex != otherInstanceIndex) {
                    return 100;
                } else {
                    return 0;
                }
            }
        };

        assertEquals(100d, distanceCalculator.distance(0, 1), DELTA);
        assertEquals(100d, distanceCalculator.distance(1, 0), DELTA);
        assertEquals(0d, distanceCalculator.distance(0, 0) , DELTA);
        assertEquals(0d, distanceCalculator.distance(1, 1), DELTA);

        final List<int[]> clusters = clusterer.cluster(distanceCalculator, 2);
        assertNotNull(clusters);
        assertEquals(2, clusters.size());
        assertEquals(1, clusters.get(0).length);
        assertEquals(1, clusters.get(1).length);
        assertEquals(0, clusters.get(0)[0]);
        assertEquals(1, clusters.get(1)[0]);

    }

    @Test
    public void fourInstanceClusteringInOneCluster() {
        // put one instance in each cluster, total two instances
        final QTClusterer clusterer = new QTClusterer(4);
        final SimilarityDistanceCalculator distanceCalculator = new MaxLinkageDistanceCalculator() {
            @Override
            public double distance(final int instanceIndex, final int otherInstanceIndex) {
                if (instanceIndex == 0 && otherInstanceIndex == 1 ||
                        instanceIndex == 1 && otherInstanceIndex == 0) {
                    return 0;    // instances 0 and 1 belong to same cluster
                } else {
                    return 10;
                }
            }
        };
        // instances 0,1,2,3 go to cluster 1 (distance(0,1)=0; distance(2,0)=10<=threshold)

        assertEquals(0d, distanceCalculator.distance(0, 1), DELTA);
        assertEquals(0d, distanceCalculator.distance(1, 0), DELTA);
        assertEquals(10d, distanceCalculator.distance(0, 0), DELTA);
        assertEquals(10d, distanceCalculator.distance(1, 1), DELTA);
        assertEquals(10d, distanceCalculator.distance(0, 2), DELTA);
        assertEquals(10d, distanceCalculator.distance(2, 0), DELTA);
        assertEquals(10d, distanceCalculator.distance(2, 3), DELTA);
        final List<int[]> clusters = clusterer.cluster(distanceCalculator, 10);
        assertNotNull(clusters);
        assertEquals("Expected one cluster", 1, clusters.size());
        final int[] cluster = clusters.get(0);
        assertEquals("First cluster must have size 4", 4, cluster.length);
        assertTrue("Instance 0 in cluster 0", ArrayUtils.contains(cluster, 0));
        assertTrue("Instance 1 in cluster 0", ArrayUtils.contains(cluster, 1));
        assertTrue("Instance 2 in cluster 0", ArrayUtils.contains(cluster, 2));
        assertTrue("Instance 3 in cluster 0", ArrayUtils.contains(cluster, 3));
    }

    @Test
    public void fourInstanceClusteringInThreeClusters() {
        // put one instance in each cluster, total two instances
        final QTClusterer clusterer = new QTClusterer(4);
        final SimilarityDistanceCalculator distanceCalculator = new MaxLinkageDistanceCalculator() {

            @Override
            public double distance(final int instanceIndex, final int otherInstanceIndex) {
                if (instanceIndex == 0 && otherInstanceIndex == 1 ||
                        instanceIndex == 1 && otherInstanceIndex == 0) {
                    return 0;    // instances 0 and 1 belong to same cluster
                } else {
                    return 11;
                }
            }
        };
        assertEquals(0d, distanceCalculator.distance(0, 1), DELTA);
        assertEquals(0d, distanceCalculator.distance(1, 0), DELTA);
        assertEquals(11d, distanceCalculator.distance(0, 2), DELTA);
        assertEquals(11d, distanceCalculator.distance(2, 0), DELTA);
        assertEquals(11d, distanceCalculator.distance(2, 3), DELTA);
        final List<int[]> clusters = clusterer.cluster(distanceCalculator, 10);
        assertNotNull(clusters);
        assertEquals("Incorrect number of clusters", 3, clusters.size());
        assertEquals("First cluster must have size 2", 2, clusters.get(0).length);
        assertEquals("Second cluster must have size 1", 1, clusters.get(1).length);
        assertEquals("Third cluster must have size 1", 1, clusters.get(2).length);
        assertEquals("Instance 0 in cluster 0", 0, clusters.get(0)[0]);
        assertEquals("Instance 1 in cluster 0", 1, clusters.get(0)[1]);
        assertEquals("Instance 2 in cluster 1", 2, clusters.get(1)[0]);
        assertEquals("Instance 3 in cluster 2", 3, clusters.get(2)[0]);
    }

    @Test
    public void fourInstanceClusteringInFourClusters() {
        // put one instance in each cluster, total two instances
        final QTClusterer clusterer = new QTClusterer(4);
        final SimilarityDistanceCalculator distanceCalculator = new MaxLinkageDistanceCalculator() {

            @Override
            public double distance(final int instanceIndex, final int otherInstanceIndex) {
                if (instanceIndex == 0 && otherInstanceIndex == 1 ||
                        instanceIndex == 1 && otherInstanceIndex == 0) {
                    return 0;    // instances 0 and 1 belong to same cluster
                } else {
                    return 10;
                }
            }
        };

        final List<int[]> clusters = clusterer.cluster(distanceCalculator, 2);
        assertNotNull(clusters);
        assertEquals("Incorrect number of clusters", 3, clusters.size());
        assertEquals("First cluster must have size 2", 2, clusters.get(0).length);
        assertEquals("Second cluster must have size 1", 1, clusters.get(1).length);
        assertEquals("Third cluster must have size 1", 1, clusters.get(2).length);
        assertEquals("Instance 0 in cluster 0", 0, clusters.get(0)[0]);
        assertEquals("Instance 1 in cluster 0", 1, clusters.get(0)[1]);
        assertEquals("Instance 2 in cluster 2", 2, clusters.get(1)[0]);
        assertEquals("Instance 3 in cluster 3", 3, clusters.get(2)[0]);
    }

    @Test
    public void zeroDistanceCalculator() {
        final QTClusterer clusterer = new QTClusterer(4);
        final SimilarityDistanceCalculator distanceCalculator = new MaxLinkageDistanceCalculator() {

            @Override
            public double distance(final int instanceIndex, final int otherInstanceIndex) {
                return 0;    // instances 0-3 belong to the same cluster
            }
        };

        final List<int[]> clusters = clusterer.cluster(distanceCalculator, 2);
        assertNotNull(clusters);
        assertEquals("Expected one cluster", 1, clusters.size());
        final int[] cluster = clusters.get(0);
        assertEquals("First cluster must have size 4", 4, cluster.length);
        assertTrue("Instance 0 in cluster 0", ArrayUtils.contains(cluster, 0));
        assertTrue("Instance 1 in cluster 0", ArrayUtils.contains(cluster, 1));
        assertTrue("Instance 2 in cluster 0", ArrayUtils.contains(cluster, 2));
        assertTrue("Instance 3 in cluster 0", ArrayUtils.contains(cluster, 3));
    }

    /**
     * This test validates that the clusterer will not throw any
     * errors when passed zero instances.
     */
    @Test
    public void zeroInstances() {
        final Clusterer iterative = new QTClusterer(0);
        final Clusterer recursive = new RecursiveQTClusterer(0);
        final SimilarityDistanceCalculator distanceCalculator =
                new MaxLinkageDistanceCalculator() {
                    @Override
                    public double distance(final int instanceIndex,
                                           final int otherInstanceIndex) {
                        return Math.abs(instanceIndex - otherInstanceIndex);
                    }
                };

        List<int[]> result = iterative.cluster(distanceCalculator, 0);
        assertNotNull(result);
        assertEquals(0, result.size());

        result = recursive.cluster(distanceCalculator, 0);
        assertNotNull(result);
        assertEquals(0, result.size());
    }

    /**
     * This test validates that the clusterer will not not allow a negative
     * instance count.
     */
    @Test (expected=IllegalArgumentException.class)
    public void illegalInstanceCount() {
        new QTClusterer(-1);
    }

    /**
     * This test validates that a dataset is clustered correctly using
     * various different values of thresholds.
     */
    @Test
    public void multipleThresholds() {
        // raw data to test
        final int[] data = {
            1, 2, 3, 3, 2, 1, 42, 43, 4, 6
        };

        // list of expected results per threshold tested
        final List[] expectedResults = new List[6];
        // threshold = 0 ( each instance in it's own cluster )
        expectedResults[0] = new ArrayList<int[]>();
        for (final int i : data) {
            expectedResults[0].add(new int[] { i });
        }

        // threshold = 1
        expectedResults[1] = new ArrayList();
        expectedResults[1].add(new int[] { 1, 1, 2, 2 });
        expectedResults[1].add(new int[] { 3, 3, 4 });
        expectedResults[1].add(new int[] { 42, 43 });
        expectedResults[1].add(new int[] { 6 });

        // threshold = 2
        expectedResults[2] = new ArrayList();
        expectedResults[2].add(new int[] { 1, 1, 2, 2, 3, 3 });
        expectedResults[2].add(new int[] { 42, 43 });
        expectedResults[2].add(new int[] { 4, 6 });

        // threshold = 3
        expectedResults[3] = new ArrayList();
        expectedResults[3].add(new int[] { 1, 1, 2, 2, 3, 3, 4 });
        expectedResults[3].add(new int[] { 42, 43 });
        expectedResults[3].add(new int[] { 6 });

        // threshold = 4 (same as 3)
        expectedResults[4] = new ArrayList(expectedResults[3]);

        // threshold = 5
        expectedResults[5] = new ArrayList();
        expectedResults[5].add(new int[] { 1, 1, 2, 2, 3, 3, 4, 6 });
        expectedResults[5].add(new int[] { 42, 43 });


        final Clusterer clusterer = new QTClusterer(data.length);
        // Distance function that deturns the difference between instances
        final SimilarityDistanceCalculator distanceCalculator =
                new MaxLinkageDistanceCalculator() {
                    @Override
                    public double distance(final int i, final int j) {
                        return Math.abs(data[i] - data[j]);
                    }
                };

        for (int i = 0; i <= 5; i++) {
            final List<int[]> clusters = clusterer.cluster(distanceCalculator, i);
            assertNotNull("Cluster at threshold " + i, clusters);

            LOGGER.debug("Iterative clusters - threshold = " + i);
            final List<int[]> expectedCluster = expectedResults[i];

            int j = 0;
            for (final int[] cluster : clusters) {
                // convert instance indexes from the cluster to data
                final int[] result = new int[cluster.length];
                for (int k = 0; k < result.length; k++) {
                    result[k] = data[cluster[k]];
                }
                LOGGER.debug(j + ":" + ArrayUtils.toString(result));
                final int[] expectedResult = expectedCluster.get(j);
                assertTrue("Cluster " + j + "with threshold " + i
                        + "does not match expected",
                        ArrayUtils.isEquals(expectedResult, result));
                j++;
            }
        }
    }
}
