/*
 * Copyright (C) 2024 Nicola Müller <nicola.felix.mueller@gmail.com>
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

import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.networkannotator.ReassortmentAnnotator;
import coalre.networkannotator.ReassortmentLogReader;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

/**
 * Turn the network into a single child by folling a user specified segment.
 * It then looks for all coalescent events where one child has a reassortment event and the other does not.
 * It then computes the cluster size for each child and prints the results to a file.
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class ClusterSizeComparison extends ReassortmentAnnotator {
    List<NetworkNode> allTrunkNodes;
    List<Double> leaveDistance;
    List<Boolean> isTrunkNode;
    int stateCount;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("sizecomp.tsv");
        double burninPercentage = 10.0;
        int segment = 0;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "segment/chromsome or plasmid to use as base tree " + segment + "\n"+
                    "-----------------------------------------\n";
        }
    }

    public ClusterSizeComparison(NetworkAnnotatorOptions options) throws IOException {

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
	      
    	int segmentCount=-1;

        int counter=1;
        
        // keeps track of the leave nodes
        List<String> leafNodes = new ArrayList<>();
        
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	ps.print("iteration\ttime\tleafsWith\tleafsWithout\n");
          	
        	boolean first = true;
	        for (Network network : logReader){	  
	        	
	        	segmentCount = network.getSegmentCount();
	        	
	        	if (first){        		
	            	for (NetworkNode networkNode : network.getNodes()){
	            		if (networkNode.isLeaf()){
	            			leafNodes.add(networkNode.getTaxonLabel());
	            		}
	        		}
	        		first = false;
	        	}
	        	
	        	
	        	// print the header for the file
	        	// calculate the cluster size for each node that has two children with different segments
                // loop over all internal nodes and check if exactly one of the children has a reasosrtment event (i.e. single child node below)
	        	for (NetworkNode n : network.getNodes()){
	        		if (n.isCoalescence()) {        			
	        			if (n.getChildEdges().get(0).hasSegments.get(options.segment) && n.getChildEdges().get(1).hasSegments.get(options.segment)) {
                            boolean leftIsSingleChild = n.getChildEdges().get(0).childNode.isReassortment();
                            boolean rightIsSingleChild = n.getChildEdges().get(1).childNode.isReassortment();
                            boolean bothAreSingleChildren = leftIsSingleChild && rightIsSingleChild;
                            if ((leftIsSingleChild || rightIsSingleChild) && !bothAreSingleChildren) {
                                int countwith = 0;
                                int countwithout = 0;
                                if (leftIsSingleChild) {
                                    countwith = calculateClusterSize(n.getChildEdges().get(0), options.segment);
                                    countwithout = calculateClusterSize(n.getChildEdges().get(1), options.segment);                                
                                }else {
                                	countwith = calculateClusterSize(n.getChildEdges().get(1), options.segment);
                                	countwithout = calculateClusterSize(n.getChildEdges().get(0), options.segment);
                                }
                                ps.print(counter + "\t" + n.getHeight() + "\t" + countwith + "\t" + countwithout + "\n");
                        }
	        		}
	        		
	        		
            }   }
	        	counter=counter+1;
	        }
	        ps.close();
        }
        
        System.out.println("\nDone!");
    }
    

	
	private int calculateClusterSize(NetworkEdge e, int segment) {
		int clusterSize = 0;
		if (e.childNode.isLeaf()) {
			clusterSize++;
		}else {
			for (NetworkEdge childEdge : e.childNode.getChildEdges()) {
                if (childEdge.hasSegments.get(segment)) {
					clusterSize += calculateClusterSize(childEdge, segment);
				}
			}
		}		
		return clusterSize;
	}
	
	
    private Tree getSingleChildTree(Network network, int segment, List<String> leafNodes){
    	// get the root of this segment tree by starting at the network root
    	// and following the segment down until we reach the first coalescent event
    	// where both children have the segment
        NetworkNode segmentRoot = findSegmentRoot(network.getRootEdge().childNode, segment);
        
        if (segmentRoot == null)
        	throw new IllegalArgumentException("root of segment tree not found");

    	Node root = getNextNode(segmentRoot, network.getSegmentCount(), segment, leafNodes);
    	   	
    	Tree tree = new Tree(root);
    	return tree;
    }
    
    private NetworkNode findSegmentRoot(NetworkNode currentNode, int segment) {
    	// Check if this is a coalescent node where both children have the segment
    	if (currentNode.isCoalescence()) {
    		List<NetworkEdge> childEdges = currentNode.getChildEdges();
    		if (childEdges.size() == 2 &&
    			childEdges.get(0).hasSegments.get(segment) &&
    			childEdges.get(1).hasSegments.get(segment)) {
    			// Found the root of the segment tree
    			return currentNode;
    		}
    	}
    	
    	// Continue following the segment down
    	for (NetworkEdge childEdge : currentNode.getChildEdges()) {
    		if (childEdge.hasSegments.get(segment)) {
    			NetworkNode result = findSegmentRoot(childEdge.childNode, segment);
    			if (result != null) {
    				return result;
    			}
    		}
    	}
    	
    	// Segment root not found in this branch
    	return null;
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
        
	    			newNode.metaDataString = "isReassortment=" + (childEdge.childNode.isReassortment()?1:0);
	        		
	            	for (int i = 0; i < nrSegments; i++)
	            		newNode.metaDataString = newNode.metaDataString + ",seg" + i + "=" + (childEdge.hasSegments.get(i)? 1 : 0);

	            	newNode.metaDataString = newNode.metaDataString + ",segsCarried=" + childEdge.hasSegments.cardinality();

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
	
    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Plasmid Tree Mapper");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel chromosomeIndexLabel = new JLabel("Index of the chromosome (or plasmid)\nto output, starts counting at 0:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        JTextField chromosomeIndex = new JTextField(20);
        chromosomeIndex.setText(Integer.toString(options.segment));
        chromosomeIndex.setEditable(true);
//        minTipDistance.setEnabled(false);        

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);

        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(chromosomeIndexLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
                        .addComponent(chromosomeIndex))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton))
                );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(chromosomeIndexLabel)
                        .addComponent(chromosomeIndex))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.segment = Integer.parseInt(chromosomeIndex.getText());
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select Reassortment Network log file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Reassortment Event Locator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "ClusterSizeComparison - Analyzes reassortment networks to compute cluster sizes.\n"
                    + "\n"
                    + "Usage: ClusterSizeComparison [options] <inputLogFile> [outputFile]\n"
                    + "\n"
                    + "Options:\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display this usage information.\n"
                    + "-burnin <percentage>     Specify the percentage of the log file to discard\n"
                    + "                         as burn-in. Default is 10%.\n"
                    + "-segment <index>         Specify the index of the segment to analyze,\n"
                    + "                         starting from 0. Default is 0.\n"
                    + "\n"
                    + "If no output file is specified, the default output file name is 'sizecomp.tsv'.\n"
                    + "The inputLogFile is mandatory and must be a valid reassortment network log file.";

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
                case "-segment":
                    if (args.length<=i+1) {
                        printUsageAndError("-segment must be one of the segments in the network.");
                    }

                    try {
                        options.segment = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
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

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new ClusterSizeComparison(options);

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