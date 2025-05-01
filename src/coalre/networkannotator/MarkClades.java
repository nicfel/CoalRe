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
public class MarkClades extends ReassortmentAnnotator {

    
    double totalLength = 0;
    int totalReassortments = 0;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile;
        File cladeFile;
        double burninPercentage = 10.0;
        int followSegment = 0;
        int[] removeSegments = new int[0];


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Clade file: " + cladeFile + "\n" +
                    "Burnin percentage: " + burninPercentage + "\n" +
                    "Follow segment: " + followSegment + "\n" +
                    "Remove segments: " + (removeSegments.length > 0 ? removeSegments[0] : "") + "\n" +
                    "-----------------------------------------\n";
                    
                    
        }
    }

    public MarkClades(NetworkAnnotatorOptions options) throws IOException {

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
             
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	int count = 0;
	        for (Network network : logReader){	  
//	        	pruneNetwork(network, options.removeSegments);
	        	labelClades(network, clades, tipNames, options.followSegment);
	        	// mark the remaining edges based on the parental edges
	        	markRemainingEdges(network, options.followSegment);
	        	markReassortmentEvents(network);
	        	ps.println(network);
	        	
	        	count++;
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }
    

        
    private void markReassortmentEvents(Network network) {
		List<NetworkNode> reassortmentNodes = network.getNodes().stream().filter(e -> e.isReassortment())
				.collect(Collectors.toList());
		
		// make a map with all reassortment nodes and keep track of the new type label w/o changing it yet
		Map<NetworkNode, String> reassortmentMap = new HashMap<>();

		for (NetworkNode r : reassortmentNodes) {
			// get the type label of the parent nodes
			List<String> typeLabels = new ArrayList<>();
			for (NetworkEdge e : r.getParentEdges()) {
				typeLabels.add(e.parentNode.getTypeLabel());				
			}
			// remove all type labels that say multiple
			
			
			// remove duplicates
			Collections.sort(typeLabels);
			for (int i = typeLabels.size() - 1; i>0; i--) {
				if (typeLabels.get(i-1).equals(typeLabels.get(i))) {
					typeLabels.remove(i);
					i--;
				}
			}
			
			// if type label size is >1 combined them using a plus sign
			String typeLabel = "";
			if (typeLabels.size() == 2) {
				for (int i = 0; i < typeLabels.size(); i++) {
					typeLabel += typeLabels.get(i);
					if (i < typeLabels.size() - 1) {
						typeLabel += "+";
					}
				}
			} else if (typeLabels.size() == 1) {
				typeLabel = typeLabels.get(0);
			} else {
				throw new IllegalArgumentException("Type label size is 0 or >2: " + typeLabels);
			}
			
			reassortmentMap.put(r, typeLabel);
		}
		
		// now, assign the type labels to the reassortment nodes
		for (NetworkNode r : reassortmentMap.keySet()) {
			r.setTypeLabel(reassortmentMap.get(r));
		}
		
		
	}



	private void markRemainingEdges(Network n, int followSegment) {
		
		// get all edges that do not have the followSegment
		List<NetworkEdge> edges = n.getEdges().stream()
				.filter(e -> e.hasSegments.get(followSegment) == false)
				.collect(Collectors.toList());
    	do {
    		// remove all edges for which childNode has a type label
			for (int i = edges.size()-1; i >= 0; i--) {
				if (edges.get(i).childNode.getTypeLabel() != null) {
					// remove edge from list
					edges.remove(i);					
				}
			}
			
			// for all remaining edges, assign the type label of the parentNODE to childNode
			for (NetworkEdge edge : edges) {
				if (edge.parentNode.getTypeLabel() != null) {
                    edge.childNode.setTypeLabel(edge.parentNode.getTypeLabel());
                }
			}    		
    	} while (edges.size() > 0);
		
	}



	private void labelClades(Network network, String[] clades, String[] tipNames, int followSegment) {
		labelNodes(network.getRootEdge(), clades, tipNames, followSegment);		
	}



	private List<String> labelNodes(NetworkEdge edge, String[] clades, String[] tipNames, int followSegment) {
		// get all leaves below rootEdge after following follow segment
		List<String> hasClades = new ArrayList<>();
		
		if (edge.childNode.isLeaf()) {
			// find the clades for this leaf
			String leaf = edge.childNode.getTaxonLabel();
			for (int i = 0; i < clades.length; i++) {
				if (leaf.contentEquals(tipNames[i])) {
					hasClades.add(clades[i]);
					edge.childNode.setTypeLabel(hasClades.get(0));

					return hasClades;
				}
			}			
		}else {
			for (NetworkEdge e : edge.childNode.getChildEdges()) {
				if (e.hasSegments.get(followSegment)) {
					hasClades.addAll(labelNodes(e, clades, tipNames, followSegment));
				}
			}
		}
		// remove duplicates in hasClades
		Collections.sort(hasClades);
		for (int i = hasClades.size() - 1; i>0; i--) {
			if (hasClades.get(i-1).equals(hasClades.get(i))) {
				hasClades.remove(i);
				i--;
			}
		}
		
		// add a new metadata entry to the edge
		if (hasClades.size()==1)
			edge.childNode.setTypeLabel(hasClades.get(0));
		else if (hasClades.size() > 1) {
			edge.childNode.setTypeLabel("multiple");
		}
		
		return hasClades;
		
		
	}



	/**
     * gets how many reticulation events happen on the trunk vs. not on the trunk
     * The trunk is define as any edge on the network that has descendents that are more than minTipDistance 
     * away from that node
     * @param network
     * @param ps
     * @param minTipDistance
     */
    private void getEvents(Network network, PrintStream ps, int count){    	
   	
		NetworkEvent prevEvent = null;

		double waitTime = 0;
		List<NetworkNode> reassortmentNodes = network.getNodes().stream()
                .filter(e -> e.isReassortment())
                .collect(Collectors.toList());
		
        for (NetworkNode r : reassortmentNodes) {   
	        ps.print(count + "\t" + r.getHeight()+  "\t" + r.getParentEdges().get(0).hasSegments +"\t" + 
	        		r.getParentEdges().get(1).hasSegments + "\n");	        
        }
    	
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
                    
                case "-followSegment":
                    if (args.length<=i+1) {
                        printUsageAndError("-removeSegments must be followed by at least one number.");
                    }

                    try {
                		options.followSegment = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("followSegment must be an integer");
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
            new MarkClades(options);

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