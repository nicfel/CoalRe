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

import coalre.distribution.NetworkEvent;
import coalre.distribution.NetworkIntervals;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import javax.swing.*;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import cern.colt.Arrays;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * marks the clades in a reassortment network by using a table with clades and sepecifying a segment to follow for clade labelling 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class getCladeHeights extends ReassortmentAnnotator {

    
    double totalLength = 0;
    int totalReassortments = 0;
    NetworkEdge mrcaRootEdge;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile;
        File cladeFile;
        double burninPercentage = 10.0;


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Clade file: " + cladeFile + "\n" +
                    "Burnin percentage: " + burninPercentage + "\n" +
                    "-----------------------------------------\n";
                    
                    
        }
    }

    public getCladeHeights(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");

        // read in the clade table as a csv file

        // Initialise reader
        ReassortmentLogReader logReader = new ReassortmentLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

        System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
	      
	    List<String> lines = null;
    	try {
    		lines = java.nio.file.Files.readAllLines(options.cladeFile.toPath());
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    	
        String[] clades = new String[lines.size()-1];
        String[] tipNames = new String[lines.size()-1];
        
		for (int i = 1; i < lines.size(); i++) {
			// split line on , or \t
			String[] split = lines.get(i).split("[,\t]");
			clades[i-1] = split[1];
			tipNames[i-1] = split[0];
		}	
		
		List<String> uniqueClades = new ArrayList<>();
		for (String clade : clades) {
			if (!uniqueClades.contains(clade)) {
				uniqueClades.add(clade);
			}
		}
		// add the sorted combination of unique clades seperated by + to unique clades
		Collections.sort(uniqueClades);
		
		
		
		
             
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	int count = 0;
        	
        	// start printing the header
        	ps.print("Sample");
        	
    		// for each clade in uniqueClades, make a list with all the tips in it
    		List<List<String>> cladeTips = new ArrayList<>();
    		for (String clade : uniqueClades) {
    			List<String> tips = new ArrayList<>();
                for (int i = 0; i < clades.length; i++) {
                    if (clades[i].contentEquals(clade)) {
                        tips.add(tipNames[i]);
                    }
                }
                cladeTips.add(tips);       
    		}


    		int sample = 0;
	        for (Network network : logReader){	  
	        	// for each unique clade, get the edge above the
	        	if (sample==0) {
					for (String clade : uniqueClades) {
		                // loop over all segments and print the clade name
						for (int i = 0; i < network.getSegmentCount(); i++) {
		                    ps.print("\tmin:" + clade + ":" + i);
		                    ps.print("\tmax:" + clade + ":" + i);
		                }
					}
					ps.println();
	        	}
	        	
	        	
	        	ps.print(sample);
				for (int c=0; c < uniqueClades.size(); c++) {
					// loop over all segments
					for (int i = 0; i < network.getSegmentCount(); i++) {
						// label the clades in the network
						double[] minMaxHeights = getMinMaxClade(network, cladeTips.get(c).toArray(new String[0]), i);
						ps.print("\t" + minMaxHeights[0]);
						ps.print("\t" + minMaxHeights[1]);
					}
				}
				ps.println();
				sample++;
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }


	private double[] getMinMaxClade(Network network, String[] tipNames, int followSegment) {
		// get the most recent edge that has all tips below it when following followSegment		
		hasClades(network.getRootEdge(), tipNames, followSegment);	
		return new double[] {mrcaRootEdge.childNode.getHeight(), mrcaRootEdge.parentNode.getHeight()};
	}

	private int hasClades(NetworkEdge edge, String[] tipNames, int followSegment) {
		// get all leaves below rootEdge after following follow segment
		int noTips = 0;
		
		if (edge.childNode.isLeaf()) {
			// check if the leaf is in the tipNames
			for (String tipName : tipNames) {
				if (edge.childNode.getTaxonLabel().equals(tipName)) {
					return 1; // found a leaf with a clade
				}
			}
			return 0;
		}else {
			for (NetworkEdge e : edge.childNode.getChildEdges()) {
				if (e.hasSegments.get(followSegment)) {
					noTips += hasClades(e, tipNames, followSegment);
				}
			}
		}
		if (noTips == tipNames.length) {
			mrcaRootEdge = edge;
			return 0;
		}
		return noTips;
	}


	public static String helpMessage =
	        "getCladeHeights - Annotates reassortment networks with min/max clade heights\n"
	      + "                   for specified clades across segments.\n"
	      + "\n"
	      + "Usage: java getCladeHeights [options]\n"
	      + "\n"
	      + "Required options:\n"
	      + "--------------------------------------------------------------\n"
	      + "-tree <file>             Input file containing sampled reassortment networks.\n"
	      + "-clade <file>            TSV or CSV file listing tip-to-clade mappings.\n"
	      + "                         Format: one line per tip with columns [tipName,cladeName]\n"
	      + "-out <file>              Output file to write clade height summaries.\n"
	      + "\n"
	      + "Optional arguments:\n"
	      + "--------------------------------------------------------------\n"
	      + "-burnin <percentage>     Percentage of networks to discard as burn-in.\n"
	      + "                         (Default: 10.0)\n"
	      + "-help                    Display this help message and exit.\n"
	      + "\n"
	      + "Output format:\n"
	      + "--------------------------------------------------------------\n"
	      + "Tab-delimited file with one row per network sample and columns:\n"
	      + "  Sample, followed by min and max height for each clade and segment:\n"
	      + "  e.g. Sample\tmin:CladeA:0\tmax:CladeA:0\tmin:CladeB:1\t...\n"
	      + "\n"
	      + "Example:\n"
	      + "  java getCladeHeights -tree networks.log -clade clades.csv -out clade_heights.tsv\n";

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
                    
				case "-tree":
					if (args.length <= i + 1)
						printUsageAndError("-treeFile must be followed by a file name.");

					options.inFile = new File(args[i + 1]);
					i += 1;
					break;
					
				case "-clade":
					if (args.length <= i + 1)
						printUsageAndError("-cladeFile must be followed by a file name.");
					options.cladeFile = new File(args[i + 1]);
					i += 1;
					break;
				
				case "-out":
					if (args.length <= i + 1)
						printUsageAndError("-out must be followed by a file name.");
					options.outFile = new File(args[i + 1]);
					i += 1;
					break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
            
			if (i >= args.length)
				break;
        }

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
            new getCladeHeights(options);

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