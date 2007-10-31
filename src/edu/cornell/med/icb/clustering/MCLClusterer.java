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

import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

/**
 * http://micans.org/mcl/
 */
public final class MCLClusterer implements Clusterer {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOGGER = Logger.getLogger(MCLClusterer.class);

    /**
     * The total number of instances to cluster.
     */
    private final int instanceCount;

    /**
     * The total number of clusters created by clustering the instances.
     */
    private int clusterCount;

    /**
     * The list of clusters.
     */
    private final IntArrayList[] clusters;

    /**
     * Default executable name for mcl.  The default assumes that the command is already
     * on the path.
     */
    private static final String DEFAULT_MCL_COMMAND = System.getProperty("MCL_COMMAND", "mcl");

    /**
     * Name of mcl executable.  This should be a fully qualified path unless it
     * is located on the execution path.
     */
    private String mclCommand = DEFAULT_MCL_COMMAND;

    /**
     * Construct a new quality threshold clusterer.
     *
     * @param numberOfInstances The number of instances to cluster.
     *                          i.e., |G| where G is the set of instances
     */
    public MCLClusterer(final int numberOfInstances) {
        this(numberOfInstances, DEFAULT_MCL_COMMAND);
    }

    public MCLClusterer(final int numberOfInstances, final String mclExecutable) {
        super();
        if (numberOfInstances < 0) {
            throw new IllegalArgumentException("Number of instances ("
                    + numberOfInstances + ") must not be negative");
        }

        instanceCount = numberOfInstances;
        clusters = new IntArrayList[numberOfInstances];
        for (int i = 0; i < numberOfInstances; i++) {
            clusters[i] = new IntArrayList();                        // NOPMD
        }

        mclCommand = mclExecutable;
    }


    /**
     * Creates clusters that are read directly from an MCL output file.  This constructor is
     * simply intended to take MCL output and conform it to the
     * {@link edu.cornell.med.icb.clustering.Clusterer} interface.  When you already have
     * an MCL output file, you would skip the call to
     * {@link #cluster(SimilarityDistanceCalculator, double)} and jump directly to
     * {@link #getClusters()} after calling this constructor.
     * @param reader The reader/file to read from
     * @throws IOException if there is a problem reading from the Reader
     */
    public MCLClusterer(final Reader reader) throws IOException {
        super();

        int numberOfInstances = 0;
        final List<IntArrayList> clusterArray = new ArrayList<IntArrayList>();

        TSVReader tsvReader = null;
        try {
            tsvReader = new TSVReader(reader);
            clusterCount = 0;
            while (tsvReader.hasNext()) {
                if (!tsvReader.isEmptyLine() && !tsvReader.isCommentLine()) {
                    tsvReader.next();
                    final IntArrayList cluster = new IntArrayList();
                    for (int i = 0; i < tsvReader.numTokens(); i++) {
                        cluster.add(tsvReader.getInt());
                        numberOfInstances++;
                    }
                    clusterArray.add(cluster);
                    clusterCount++;
                } else {
                    tsvReader.skip();
                }
            }
        } finally {
            if (tsvReader != null) {
                try {
                    tsvReader.close();
                } catch (IOException ioe) {
                    // ignore
                }
            }
        }

        instanceCount = numberOfInstances;
        clusters = new IntArrayList[clusterCount];
        int i = 0;
        for (final IntArrayList cluster : clusterArray) {
            clusters[i] = cluster;
            i++;
        }
    }

    /**
     * Groups instances into clusters. Returns the indices of the instances that belong to a cluster
     * as an int array in the list result.
     *
     * @param calculator The {@link SimilarityDistanceCalculator}
     * that should be used when clustering
     * @param qualityThreshold The QT clustering algorithm quality threshold (d)
     * @return The list of clusters.
     */
    public List<int[]> cluster(final SimilarityDistanceCalculator calculator,
                               final double qualityThreshold) {
        if (mclCommand == null) {
            throw new IllegalStateException("mcl command not set!");
        }

        // reset cluster results
        clusterCount = 0;
        for (int i = 0; i < instanceCount; i++) {
            clusters[i].clear();
        }

        BufferedReader br = null;
        try {
            final File mclInputFile = File.createTempFile("mcl-input", ".txt");
            writeMCLInputFile(mclInputFile, calculator, qualityThreshold);

            final File mclOutputFile = File.createTempFile("mcl-output", ".txt");
            final String[] command = {
                    mclCommand, mclInputFile.getAbsolutePath(), "--abc",
                    "-o", mclOutputFile.getAbsolutePath()
            };

            LOGGER.info("Executing: " + ArrayUtils.toString(command));

            final ProcessBuilder builder = new ProcessBuilder(command);
            builder.redirectErrorStream(true);
            final Process process = builder.start();
            final InputStream is = process.getInputStream();
            final InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr);
            String line;
            while ((line = br.readLine()) != null) {
                LOGGER.info(line);
            }
            process.waitFor();

            LOGGER.info("Program terminated!");

            readMCLOutputFile(mclOutputFile);
        } catch (IOException e) {
            LOGGER.error("Counldn't create MCL file", e);
            throw new ClusteringException("Counldn't create MCL file", e);
        } catch (InterruptedException e) {
            LOGGER.error("Interrupted!", e);
            Thread.currentThread().interrupt();
        } finally {
            IOUtils.closeQuietly(br);
        }

        return getClusters();
    }

    /**
     * Returns the list of clusters produced by clustering.
     *
     * @return A list of integer arrays, where each array represents a cluster and contains the
     *         index of the instance that belongs to a given cluster.
     */
    public List<int[]> getClusters() {
        final List<int[]> result = new ArrayList<int[]>(clusterCount);
        for (int i = 0; i < clusterCount; i++) {
            result.add(clusters[i].toIntArray());
        }
        return result;
    }

    /**
     * Writes instance pairs that are within the distance threshold to a file
     * in the proper format for mcl to read.
     *
     * @param file The file to write the instances to be clustered
     * @param calculator The {@link SimilarityDistanceCalculator}
     * that should be used when clustering
     * @param qualityThreshold The QT clustering algorithm quality threshold (d)
     * @throws java.io.IOException if the file cannot be written to
     */
    private void writeMCLInputFile(final File file,
                                   final SimilarityDistanceCalculator calculator,
                                   final double qualityThreshold) throws IOException {
        LOGGER.debug("Writing mcl file to: " + file.getName());
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(file);
            for (int i = 0; i < instanceCount; i++) {
                // note that we don't need to start from 0 here otherwise there will be duplicates
                for (int j = i; j < instanceCount; j++) {
                    final double distance = calculator.distance(i, j);
                    // if the distance is within the threshold, add this pair
                    // we also want to explicitly include each instance with itself
                    if (distance <= qualityThreshold || i == j) {
                        writer.printf("%d\t%d\t%f" + IOUtils.LINE_SEPARATOR, i, j, 1.0);
                    }
                }
            }
        } finally {
            IOUtils.closeQuietly(writer);
        }
    }

    private void readMCLOutputFile(final File file) throws IOException {
        TSVReader tsvReader = null;
        FileReader fileReader = null;
        try {
            fileReader = new FileReader(file);
            tsvReader = new TSVReader(fileReader);
            while (tsvReader.hasNext()) {
                if (!tsvReader.isEmptyLine() && !tsvReader.isCommentLine()) {
                    tsvReader.next();
                    final IntArrayList cluster = clusters[clusterCount];
                    for (int i = 0; i < tsvReader.numTokens(); i++) {
                        cluster.add(tsvReader.getInt());
                    }
                    clusterCount++;
                } else {
                    tsvReader.skip();
                }
            }
        } finally {
            IOUtils.closeQuietly(fileReader);
            if (tsvReader != null) {
                tsvReader.close();
            }
        }
    }
}
