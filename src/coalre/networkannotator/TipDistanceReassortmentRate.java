/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package coalre.networkannotator;

import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Computes the reassortment rates based on the distance to the nearest tip 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class TipDistanceReassortmentRate extends ReassortmentAnnotator {

    

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("reassortment_distances.txt");
        double burninPercentage = 10.0;
        double maxTipDistance = 1.0;
        int[] removeSegments = new int[0];


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
            		"maximal distance to a tip to be considered in rate calculation:" + 
                    + maxTipDistance;
        }
    }

    public TipDistanceReassortmentRate(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader
        ReassortmentLogReader logReader = new ReassortmentLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

	      System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
        
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
	        for (Network network : logReader){	        	
	        	pruneNetwork(network, options.removeSegments);
	        	computeReassortmentLeaveDist(network, ps, options.maxTipDistance);	        	
	        	ps.print("\n");
	        }
	        System.exit(0);
	        ps.close();
        }
        System.out.println("\nDone!");
    }
    

        
    /**
     * gets how many reticulation events happen on the trunk vs. not on the trunk
     * The trunk is define as any edge on the network that has descendents that are more than minTipDistance 
     * away from that node
     * @param network
     * @param ps
     * @param minTipDistance
     */
    private void computeReassortmentLeaveDist(Network network, PrintStream ps, double minTipDistance){    	
        
        // get the length of the network        
        List<NetworkEdge> allEdges = network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .collect(Collectors.toList());

        
        double totalLength = 0;
        int totalReassortments = 0;
        
        // for each node, get the distance to the nearest tip (inefficient due to visiting nodes multiple times)
		for (NetworkNode n : network.getNodes()){
			// get the distance to the nearest tip of the node below the edge
			double minTipDist = Double.POSITIVE_INFINITY;
			if (n.getParentEdges().get(0).isRootEdge())
				continue;
			
			if (!n.isLeaf()) {
				for (NetworkEdge children : n.getChildEdges()) {
					minTipDist = Math.min(minTipDist, getMinTipDistance(children));
				}
			}else {
				minTipDist = 0;
			}
				

			
			if (minTipDist >= minTipDistance)
				continue;

			if (n.isReassortment())
				totalReassortments++;

			
			// add the parents to the total length weighted by the observation probability given the segments
			for (NetworkEdge parent : n.getParentEdges()) {
				if (parent.getLength()+minTipDist >= minTipDistance) {
					// only count the part of the edge that is above the minTipDistance
					double correctedLength = parent.getLength() - (minTipDistance - minTipDist);
					totalLength += correctedLength*(1.0-Math.pow(0.5, parent.hasSegments.cardinality()-1));
				}else {					
					totalLength += parent.getLength()*(1.0-Math.pow(0.5, parent.hasSegments.cardinality()-1));
				}
			}
        }
        				

        ps.print(totalReassortments + "\t" + totalLength + "\t" + totalReassortments/totalLength);

    }
    
    private double getMinTipDistance(NetworkEdge edge) {
    	
    	double minTipDist = Double.POSITIVE_INFINITY;
		if (edge.childNode.isLeaf())
			return edge.getLength();

		for (NetworkEdge children : edge.childNode.getChildEdges()) {
            minTipDist = Math.min(minTipDist, getMinTipDistance(children));
    	}
		
    	return minTipDist + edge.getLength();
    }
        
    public static String helpMessage =
            "TrunkReassortment - counts how many reassortment events happened on trunk and non-trunk nodes.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-trunkDefinition {MostRecentSample, TipDistance} Choose trunk definition method.\n"
                    + "                         (default MostRecentSample)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-minTipDistance     		minimum distance between internal network node\n"
                    + "                         and tip node such that the internal node is considered trunk.\n"
                    + "                         If not  specified, the trunk is any node between samples\n"
                    + "                         height=0 and the root.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'reassortment_distances.txt'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve TrunkReassortment options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, NetworkAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;


                case "-maxTipDistance":
                    if (args.length<=i+1) {
                        printUsageAndError("-minTipDistance must be followed by a number.");
                    }

                    try {
                        options.maxTipDistance =
                                Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("minTipDistance must be a positive number. ");
                     }

                    i += 1;
                    break;
                    
                case "-removeSegments":
                    if (args.length<=i+1) {
                        printUsageAndError("-removeSegments must be followed by at least one number.");
                    }

                    try {
                    	String[] argarray = args[i + 1].split(",");
                    	options.removeSegments = new int[argarray.length];
                    	for (int j = 0; j < argarray.length; j++)
                    		options.removeSegments[j] = Integer.parseInt(argarray[j]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
                     }

                    i += 1;
                    break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	NetworkAnnotatorOptions options = new NetworkAnnotatorOptions();


        getCLIOptions(args, options);
        

        // Run ACGAnnotator
        try {
            new TipDistanceReassortmentRate(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}