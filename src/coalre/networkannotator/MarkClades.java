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
        int printRelativeToSegment = -1;


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Clade file: " + cladeFile + "\n" +
                    "Burnin percentage: " + burninPercentage + "\n" +
                    "Follow segment: " + followSegment + "\n" +
                    "Remove segments: " + (removeSegments.length > 0 ? removeSegments[0] : "") + "\n" +
                    "Print relative to segment: " + printRelativeToSegment + " (if -1 the whole network is printed to file)\n" +
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
		
		List<String> uniqueClades = new ArrayList<>();
		for (String clade : clades) {
			if (!uniqueClades.contains(clade)) {
				uniqueClades.add(clade);
			}
		}
		// add the sorted combination of unique clades seperated by + to unique clades
		Collections.sort(uniqueClades);
		
		List<String> combinedClades = new ArrayList<>();
		for (int i = 0; i < uniqueClades.size(); i++) {
			for (int j = i + 1; j < uniqueClades.size(); j++) {
				combinedClades.add(uniqueClades.get(i) + "+" + uniqueClades.get(j));
			}
		}
		for (int i = 0; i < uniqueClades.size(); i++) {
			combinedClades.add(uniqueClades.get(i) + "+other");
			
		}

		
		
		uniqueClades.addAll(combinedClades);
		uniqueClades.add("other");
             
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	int count = 0;
            List<String> leafNodes = new ArrayList<>();

	        for (Network network : logReader){	  
//	        	pruneNetwork(network, options.removeSegments);
	        	labelClades(network, clades, tipNames, options.followSegment);
	        	// mark the remaining edges based on the parental edges
	        	markRemainingEdges(network, options.followSegment);
	        	markReassortmentEvents(network);
	        	
				if (options.printRelativeToSegment >= 0) {
					if (count==0) {
		            	for (NetworkNode networkNode : network.getNodes()){
		            		if (networkNode.isLeaf()){
		            			leafNodes.add(networkNode.getTaxonLabel());
		            		}
		        		}

					}
					// build a tree by following segment options.printRelativeToSegment. Use the marked network to highlight reassortment events along branches
					Tree tree = getSingleChildTree(network, options.printRelativeToSegment, leafNodes);
					
					Node noSingleRoot = convertToNonSingleChildTree(tree.getRoot(), uniqueClades);	
					
					Tree newTree = new Tree(noSingleRoot);
					
					ps.println(noSingleRoot +";");
					
				}else {
					ps.println(network);
				}     	
	        	
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
			edge.childNode.setTypeLabel("other");
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
    
    public static final String helpMessage =
    	    "MarkClades - Annotates reassortment networks with clade labels based on segment paths.\n" +
    	    "\n" +
    	    "Usage:\n" +
    	    "  java MarkClades [options]\n" +
    	    "\n" +
    	    "Required options:\n" +
    	    "  -tree <file>              Input BEAST2 log file with sampled reassortment networks.\n" +
    	    "  -clade <file>             TSV or CSV file with columns: <tipName>,<cladeLabel>.\n" +
    	    "  -out <file>               Output file to write annotated networks or trees.\n" +
    	    "\n" +
    	    "Optional options:\n" +
    	    "  -burnin <percent>         Percentage of samples to discard as burn-in (default: 10.0).\n" +
    	    "  -followSegment <index>    Segment index to follow for clade inference (default: 0).\n" +
    	    "  -removeSegments <list>    Comma-separated list of segments to remove (e.g., 1,2).\n" +
    	    "  -printSegment <index>     Print segment-specific tree instead of full network.\n" +
    	    "  -help                     Print this message and exit.\n" +
    	    "\n" +
    	    "Output:\n" +
    	    "  If -printSegment is provided, a tree with clade and reassortment annotations is output.\n" +
    	    "  Otherwise, the full annotated network is written to the output file.";

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
                case "-printSegment":
					if (args.length <= i + 1) {
						printUsageAndError("-printRelativeToSegment must be followed by a number.");
					}

					try {
						options.printRelativeToSegment = Integer.parseInt(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing printRelativeToSegment.");
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
    
    
    private Tree getSingleChildTree(Network network, int segment, List<String> leafNodes){
    	// get teh root of this segment tree
        List<NetworkNode> rootEdge = network.getNodes().stream()
                .filter(e -> e.isCoalescence())
                .filter(e -> e.getChildEdges().get(0).hasSegments.get(segment))
                .filter(e -> e.getChildEdges().get(1).hasSegments.get(segment))
                .collect(Collectors.toList());
        
        double height = 0;
		for (NetworkNode node : rootEdge) {
			if (node.getHeight() > height) {
				height = node.getHeight();
			}
		}
		NetworkNode rootNode = rootEdge.get(0);
		for (NetworkNode node : rootEdge) {
			if (node.getHeight() == height) {
				rootNode = node;
            }
		}
        
        
    	Node root = getNextNode(rootNode, network.getSegmentCount(), segment, leafNodes);
    	   	
    	Tree tree = new Tree(root);
    	return tree;
    }
    
    private Node getNextNode(NetworkNode networkNode, int nrSegments, int segment, List<String> leafNodes){
    	
    	Node node = new Node();
    	
    	if (networkNode.isLeaf()){    		
    		node.setHeight(networkNode.getHeight());
    		node.setID(networkNode.getTaxonLabel());
    		node.setNr(leafNodes.indexOf(networkNode.getTaxonLabel()));
    		
    		return node;
    	}else{    	
    		
        	node.setHeight(networkNode.getHeight());
	    	List<Node> newNodes = new ArrayList<>();
	    	for (NetworkEdge childEdge : networkNode.getChildEdges()){
	    		if (childEdge.hasSegments.get(segment)){
	    			Node newNode = getNextNode(childEdge.childNode, nrSegments, segment, leafNodes);
	    			
	    			newNode.metaDataString = "isReassortment=" + (childEdge.childNode.isReassortment() ? 1:0);
	        		
	            	for (int i = 0; i < nrSegments; i++)
	            		newNode.metaDataString = newNode.metaDataString + ",seg" + i + "=" + (childEdge.hasSegments.get(i)? 1 : 0);

	            	newNode.metaDataString = newNode.metaDataString + ",segsCarried=" + childEdge.hasSegments.cardinality();
	            	newNode.metaDataString = newNode.metaDataString + ",type" + "=" + childEdge.childNode.getTypeLabel();
	            	

	    			newNodes.add(newNode);
	    		}    		
	    	}
	    	if (networkNode.isCoalescence()){
	    		node.metaDataString = "isReassortment=0";
	    		if (newNodes.size()==2){
	    			node.setLeft(newNodes.get(0));
	    			node.setRight(newNodes.get(1));
	    			
		    		newNodes.get(0).setParent(node);
		    		newNodes.get(1).setParent(node);

	    		}else{
	    			node.setLeft(newNodes.get(0));
		    		newNodes.get(0).setParent(node);
	    		}
	    	}else if (networkNode.isReassortment()){
	    		node.setLeft(newNodes.get(0));
	    		newNodes.get(0).setParent(node);
	    		node.metaDataString = "isReassortment=0";
	    	}
	    	return node;
    	}   	
    }
    
    private Node convertToNonSingleChildTree(Node n, List<String> uniqueClades) {
    	Node node = new Node();
    	
    	if (n.isLeaf()){    		
    		node.setHeight(n.getHeight());
    		node.setID(n.getID());
    		node.setNr(n.getNr());
    		
			String[] event = n.metaDataString.split(",");
			String type="";
			for (String e : event) {
				if (e.startsWith("type")) {
					type=e.split("=")[1];
				}
			}
			

    		node.metaDataString = "type=" +type;
    		
    		return node;
    	}else if (n.getChildCount() == 2){    	
        	node.setHeight(n.getHeight());
        	List<Node> newNodes = new ArrayList<>();
        	for (Node child : n.getChildren()){
        		// check if the node is a single child node
        		int reassortmentEvents = 0;
        		List<String> events = new ArrayList<>();
				while (child.getChildCount() == 1) {
					String[] event = child.metaDataString.split(",");
					for (String e : event) {
						if (e.startsWith("type")) {
							events.add(e.split("=")[1]);
						}
					}
					
					reassortmentEvents++;
					child = child.getChild(0);

				}
        		Node newNode = convertToNonSingleChildTree(child, uniqueClades);
				// count each unique event
				if (!newNode.isLeaf()) {
					newNode.metaDataString = "Events=" + reassortmentEvents;
				}else {
					newNode.metaDataString = newNode.metaDataString+",Events=" + reassortmentEvents;                
				}
				
				int[] eventCount = new int[uniqueClades.size()];
				for (int i = 0; i < events.size(); i++) {
					eventCount[uniqueClades.indexOf(events.get(i))]++;						
				}	

				for (int i = 0; i < eventCount.length; i++) {
					newNode.metaDataString += ", " + uniqueClades.get(i) + "=" + eventCount[i];
				}

        		newNodes.add(newNode);
        	}
        	
//    		node.metaDataString = "isReassortment=0";
    		if (newNodes.size()==2){
    			node.setLeft(newNodes.get(0));
    			node.setRight(newNodes.get(1));
    			
	    		newNodes.get(0).setParent(node);
	    		newNodes.get(1).setParent(node);
    		}
	    	return node;
    	}
    	
    	return null;
	}


    
    private void printTree(Tree tree, PrintStream out, int sample){
	    out.print("tree STATE_" + sample + " = ");
	    // Don't sort, this can confuse CalculationNodes relying on the tree
	    //tree.getRoot().sort();
	    final int[] dummy = new int[1];
	    final String newick = tree.getRoot().toSortedNewick(dummy, true);
	    out.print(newick);
	    out.print(";");
    }
    
    
}