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

import beast.app.treeannotator.CladeSystem;
//import bacter.Conversion;
//import bacter.ConversionGraph;
//import bacter.Locus;
//import bacter.util.BacterACGLogReader;
import beast.app.util.Utils;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.math.statistic.DiscreteStatistics;
import coalre.distribution.NetworkIntervals;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.network.SegmentTreeInitializer;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A rewrite of TreeAnnotator that outputs reassortment distances 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class ReassortmentDistance {

    private enum SummaryStrategy { MEAN, MEDIAN }

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("reassortment_distances.txt");
        double burninPercentage = 10.0;
        double convSupportThresh = 50.0;
        SummaryStrategy summaryStrategy = SummaryStrategy.MEAN;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "Node height and conv. site summary: " + summaryStrategy;
        }
    }

    public ReassortmentDistance(NetworkAnnotatorOptions options) throws IOException {

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


        
        boolean first = true;
        int segments = -1;
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
	        for (Network network : logReader){
	        	computeReassortmenDistance(network, ps);
	        	ps.print("\n");
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }
    
    private void computeReassortmenDistance(Network network, PrintStream ps){
    	// get all reassortment nodes
    	
//    	List<NetworkEdge> allEdges = network.getEdges().stream()
//    			.filter(e -> !e.isRootEdge())
//    			.collect(Collectors.toList());
//    	double networkLength = 0.0;
//    	
//    	for( NetworkEdge edge : allEdges)
//    		networkLength += edge.getLength();
    	
        List<NetworkEdge> sourceEdges = network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.isReassortment())
                .collect(Collectors.toList());
        // get the for each segment the MRCA between the reassortment node and the attachment point
        List<NetworkEdge> parentEdges = new ArrayList<>();
//        ps.print("\t");
        for (NetworkEdge edge : sourceEdges){
        	// get the reassortment node
        	NetworkNode node = edge.parentNode;
        	BitSet segmentsLeft = node.getParentEdges().get(0).hasSegments;
        	BitSet segmentsRight = node.getParentEdges().get(1).hasSegments;
        	

        	// Get the common ancestor between each pair of segments that was split in the reassortment event
        	ps.print("\t"); 
        	boolean isFirst = true; 
        	for (int i = 0; i < segmentsLeft.size(); i++){
        		for (int j = 0; j < segmentsRight.size();j++){
        			if (segmentsLeft.get(i) && segmentsRight.get(j)){
        				NetworkNode left = node.getParentEdges().get(0).parentNode;
        				NetworkNode right =  node.getParentEdges().get(1).parentNode;
        				// track both segments up the network to see when they share a common ancestor
        				while (!left.equals(right)){
        					if (left.getHeight()>right.getHeight()){
        						right = getParentNode(right, j);
        					}else{
        						left = getParentNode(left,i);
        					}        						
        				}
        				if (isFirst)
        					isFirst = false;
        				else
        					ps.print(",");
        					
        				ps.print(i + "-" + j + ":" + (right.getHeight()-node.getHeight()));
    				}
    			}
        	}	
        }
       

    }
    
    private NetworkNode getParentNode(NetworkNode node, int segment){
    	for (NetworkEdge edge : node.getParentEdges()){
    		if (edge.hasSegments.get(segment)){
    			node = edge.parentNode;
    			return node;
    		}
    	}
    	return null;
    }
    
    private int getReticulationCount(Network network){
    	Set<NetworkNode> networkNodes  = network.getNodes();
    	int reticulationCount = 0;
    	for (NetworkNode networkNode : networkNodes)
    		if (networkNode.isReassortment())
    			reticulationCount++;
    	
    	return reticulationCount;
    }
    
    private boolean[] summaryReticulationCountNetworks(SummaryStrategy summaryStrategy, int[] reticulationCount){
    	boolean[] summaryNetwork = new boolean[reticulationCount.length];
		double val = 0.0;
    	if (summaryStrategy == SummaryStrategy.MEAN){
    		for (int i = 0; i<reticulationCount.length; i++) val += reticulationCount[i];
    		val /= reticulationCount.length;
    	}else{
    		System.err.println("median not coded yet");
        	System.exit(0); 		
    	}
    	
    	// get the networks with the closest number of reticulation counts
    	double min_dist = Double.POSITIVE_INFINITY;
    	for (int i = 0; i < summaryNetwork.length; i++){
    		double dist = Math.abs(reticulationCount[i] - val);
			if (dist < min_dist)
				min_dist = dist;
    	}
    	// get all networks that are min_dist away
    	for (int i = 0; i < summaryNetwork.length; i++){
    		if (Math.abs(reticulationCount[i] - val)<(min_dist + 1e-6))
    			summaryNetwork[i] = true;
    	}   	
    	
    	return summaryNetwork;	
    }
    
    private void removeNonGeneticSegments2(Network network){
    	// remove segments from edges if they are "above" the segment tree root
    	Set<NetworkNode> networkNodes  = network.getNodes();
    	
    	double[] maxLeafHeight = new double[network.getSegmentCount()];

    	double[] segmentTreeRootHeight = new double[network.getSegmentCount()];
    	
    	
    	// get the leave with the largest height
    	for (NetworkNode networkNode : networkNodes){
    		if (networkNode.isLeaf()){
    			for (int i = 0; i < network.getSegmentCount(); i++){
    				if (networkNode.getParentEdges().get(0).hasSegments.get(i)){
    					if (networkNode.getHeight() > maxLeafHeight[i])
    						maxLeafHeight[i] = networkNode.getHeight();
    				}
    			}
    		}
    	}
    	
   	
    	
    	System.out.println(network.getExtendedNewick());
    	System.out.println(Arrays.toString(maxLeafHeight));

    	System.exit(0);
    	
    	
//    	network.get
//    	network.getRootEdge()
    	
    }
    
    private void collapseReticulations(Network network){
    	// remove all loops
    	Set<NetworkNode> networkNodes  = network.getNodes();
    	for (NetworkNode networkNode : networkNodes){
    		if (networkNode.isReassortment()){
    			// check if the two parent edges have the same parent
    			if (networkNode.getParentEdges().get(0).parentNode.equals(networkNode.getParentEdges().get(1).parentNode)){
    				// remove the reassortment node
    				NetworkEdge networkEdge = networkNode.getChildEdges().get(0);
    				NetworkNode newParent = networkNode.getParentEdges().get(0).parentNode.
    						getParentEdges().get(0).parentNode;

    				networkEdge.parentNode = newParent;
    				newParent.removeChildEdge(networkNode.getParentEdges().get(0).parentNode.
    						getParentEdges().get(0));
    				newParent.addChildEdge(networkEdge);
    			}
    				
    		}
    	}
    	
		// merge pairs of reassortment nodes if the have the same parent node
		boolean cont = false;
		do{
			cont  = collapseEdges(network);
		}while(cont);
    }
    
    private boolean collapseEdges(Network network){
    	Set<NetworkNode> networkNodes = network.getNodes();
    	
		// collapse reassortment events that end up on the same branch
    	for (NetworkNode networkNode : networkNodes){
    		// checks if the node is the root node
    		if (!networkNode.getParentEdges().get(0).isRootEdge()){
	    		if (networkNode.isReassortment()){
	    			// check if there is a group of reassortment events that are not separated by coalescent events
	    			List<NetworkNode> reassortmentGroup = new ArrayList<NetworkNode>();
	    			reassortmentGroup.add(networkNode);
	    			addNeighbouringReassortmenEvents(reassortmentGroup, networkNode);
	    			if (reassortmentGroup.size()>1){   			
	    				for (int a = 0; a < reassortmentGroup.size(); a++){
	    					for (int b = a+1; b < reassortmentGroup.size(); b++){
	    						if (reassortmentGroup.get(a).getParentEdges().get(0).parentNode.
	    								equals(reassortmentGroup.get(b).getParentEdges().get(0).parentNode)){
	    							// remove the coalescent event
	    							removeCoalescentReassortmentEvent(reassortmentGroup.get(a).getParentEdges().get(0).parentNode,
	    									reassortmentGroup.get(a).getParentEdges().get(0), reassortmentGroup.get(b).getParentEdges().get(0));
	    							return true;
	    						}else if (reassortmentGroup.get(a).getParentEdges().get(0).parentNode.
	    								equals(reassortmentGroup.get(b).getParentEdges().get(1).parentNode)){
	    							removeCoalescentReassortmentEvent(reassortmentGroup.get(a).getParentEdges().get(0).parentNode,
	    									reassortmentGroup.get(a).getParentEdges().get(0), reassortmentGroup.get(b).getParentEdges().get(1));
	    							return true;
	    						}else if (reassortmentGroup.get(a).getParentEdges().get(1).parentNode.
	    								equals(reassortmentGroup.get(b).getParentEdges().get(0).parentNode)){
	    							removeCoalescentReassortmentEvent(reassortmentGroup.get(a).getParentEdges().get(1).parentNode,
	    									reassortmentGroup.get(a).getParentEdges().get(1), reassortmentGroup.get(b).getParentEdges().get(0));
	    							return true;
	    						}else if (reassortmentGroup.get(a).getParentEdges().get(1).parentNode.
	    								equals(reassortmentGroup.get(b).getParentEdges().get(1).parentNode)){
	    							removeCoalescentReassortmentEvent(reassortmentGroup.get(a).getParentEdges().get(1).parentNode,
	    									reassortmentGroup.get(a).getParentEdges().get(1), reassortmentGroup.get(b).getParentEdges().get(1));
	    							return true;
	    						}
	    					}
	    				}		
	    			}		    			
	    		}
    		}
    	}
		return false;
    }


	private void addNeighbouringReassortmenEvents(List<NetworkNode> reassortmentGroup, NetworkNode networkNode){
    	if (!networkNode.isReassortment())
    		return;
    	
    	if (networkNode.getParentEdges().get(0).parentNode.isReassortment()){
    		reassortmentGroup.add(networkNode.getParentEdges().get(0).parentNode);
    		addNeighbouringReassortmenEvents(reassortmentGroup, networkNode.getParentEdges().get(0).parentNode);
    	}
    	
    	if (networkNode.getParentEdges().get(1).parentNode.isReassortment()){
    		reassortmentGroup.add(networkNode.getParentEdges().get(1).parentNode);
    		addNeighbouringReassortmenEvents(reassortmentGroup, networkNode.getParentEdges().get(1).parentNode);
    	}
    }
	
    private void removeCoalescentReassortmentEvent(NetworkNode parentNode, NetworkEdge networkEdge1,
			NetworkEdge networkEdge2) {
    	
    	if (!parentNode.isCoalescence())
    		throw new IllegalArgumentException("should be coal event");
    	
		// TODO make the new reassortment event the average height
		if (networkEdge2.childNode.getHeight()>networkEdge1.childNode.getHeight()){
			// remove the coalescent event
			networkEdge1.childNode.addParentEdge(parentNode.getParentEdges().get(0));
			parentNode.getParentEdges().get(0).childNode = networkEdge1.childNode;
			networkEdge1.childNode.removeParentEdge(networkEdge1);
			
			
			NetworkNode childNode = networkEdge2.childNode;
			
			// remove the other reassortment event
			childNode.removeParentEdge(networkEdge2);			
			childNode.getChildEdges().get(0).childNode.addParentEdge(childNode.getParentEdges().get(0));
			childNode.getChildEdges().get(0).childNode.removeParentEdge(childNode.getChildEdges().get(0));
			
			
		}else{
			// remove the coalescent event
			networkEdge2.childNode.addParentEdge(parentNode.getParentEdges().get(0));
			parentNode.getParentEdges().get(0).childNode = networkEdge2.childNode;
			networkEdge2.childNode.removeParentEdge(networkEdge2);
			
			NetworkNode childNode = networkEdge1.childNode;
			
			// remove the other reassortment event
			childNode.removeParentEdge(networkEdge1);			
			childNode.getChildEdges().get(0).childNode.addParentEdge(childNode.getParentEdges().get(0));
			childNode.getChildEdges().get(0).childNode.removeParentEdge(childNode.getChildEdges().get(0));			
		}

		
	}
    
    private void removeNonGeneticSegments(Network network){
    	for (int i = 0; i < network.getSegmentCount(); i++){
    		network.getRootEdge().hasSegments.set(i, false);
    		removeSegmentsFromEdge(i, network.getRootEdge().childNode);
//    		System.exit(0);
    	}
    }
    private void removeSegmentsFromEdge(int segment, NetworkNode networkNode){
    	List<NetworkEdge> childEdge = networkNode.getChildEdges();
    	
    	if (networkNode.isCoalescence()){
    		int c = 0;
    		for (NetworkEdge edge : childEdge){
    			if (edge.hasSegments.get(segment))
    				c++;
    		}
    		if (c==1){
	    		for (NetworkEdge edge : childEdge){
	    			if (edge.hasSegments.get(segment)){
	    				edge.hasSegments.set(segment, false);
	    				removeSegmentsFromEdge(segment, edge.childNode);
	    			}
	    		}
    		}
    	}else if (networkNode.isReassortment()){
    		for (NetworkEdge edge : childEdge){
    			if (edge.hasSegments.get(segment)){
    				edge.hasSegments.set(segment, false);
    				removeSegmentsFromEdge(segment, edge.childNode);
    			}
    		}   		
    	}else{
    	}
    		
    }

    private void printBitSetHeader(PrintStream ps, BitSet bitSet) {
        boolean isFirst = true;
        for (int b = bitSet.nextSetBit(0); b >= 0; b = bitSet.nextSetBit(b+1)) {

            if (isFirst)
                isFirst = false;
            else
                ps.print("_");

            ps.print(b/2);

            if (b == Integer.MAX_VALUE) {
                break; // or (i+1) would overflow
            }
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
        dialog.setTitle("Reassortment Network Annotator");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel summaryMethodLabel = new JLabel("Position summary method:");
        JLabel thresholdLabel = new JLabel("Posterior conversion support threshold:");
        JCheckBox geneFlowCheckBox = new JCheckBox("Record gene flow");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        JTextField gfOutFilename = new JTextField(20);
        gfOutFilename.setEditable(false);
        gfOutFilename.setEnabled(false);
        JButton gfOutFileButton = new JButton("Choose File");
        gfOutFileButton.setEnabled(false);

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);

        JComboBox<SummaryStrategy> heightMethodCombo = new JComboBox<>(SummaryStrategy.values());

        JSlider thresholdSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.convSupportThresh));
        thresholdSlider.setMajorTickSpacing(50);
        thresholdSlider.setMinorTickSpacing(10);
        thresholdSlider.setPaintTicks(true);
        thresholdSlider.setPaintLabels(true);
        thresholdSlider.setSnapToTicks(true);

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
                        .addComponent(summaryMethodLabel)
                        .addComponent(thresholdLabel)
                        .addComponent(geneFlowCheckBox))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
                        .addComponent(heightMethodCombo)
                        .addComponent(thresholdSlider)
                        .addComponent(gfOutFilename))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton)
                        .addComponent(gfOutFileButton)));

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
                        .addComponent(summaryMethodLabel)
                        .addComponent(heightMethodCombo,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(thresholdLabel)
                        .addComponent(thresholdSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(geneFlowCheckBox)
                        .addComponent(gfOutFilename)
                        .addComponent(gfOutFileButton)));

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.convSupportThresh = thresholdSlider.getValue();
            options.summaryStrategy = (SummaryStrategy)heightMethodCombo.getSelectedItem();
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
            inFileChooser.setDialogTitle("Select ACG log file to summarize");
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

        geneFlowCheckBox.addActionListener(e -> {
            boolean newValue = geneFlowCheckBox.isSelected();
            gfOutFilename.setEnabled(newValue);
            gfOutFileButton.setEnabled(newValue);
        });

        JFileChooser gfOutFileChooser = new JFileChooser();
        gfOutFileButton.addActionListener(e -> {
            gfOutFileChooser.setDialogTitle("Select gene flow output file name.");
            if (options.inFile != null)
                gfOutFileChooser.setCurrentDirectory(options.inFile);
            else
                gfOutFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            int returnVal = gfOutFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                gfOutFilename.setText(gfOutFileChooser.getSelectedFile().getName());
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
        frame.setTitle("ACGAnnotator");
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
            "ACGAnnotator - produces summaries of Bacter ACG log files.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-positions {mean,median} Choose position summary method.\n"
                    + "                         (default mean)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-threshold percentage    Choose minimum posterior probability\n"
                    + "                         for including conversion in summary.\n"
                    + "                         (Default 50%)\n"
                    + "-recordGeneFlow gfFile   Record posterior distribution of gene\n"
                    + "                         flow in given file.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'summary.tree'.";

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
     * Retrieve ACGAnnotator options from command line.
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

                case "-positions":
                    if (args.length<=i+1) {
                        printUsageAndError("-positions must be followed by either 'MEAN' or 'MEDIAN'.");
                    }

                    if (args[i+1].toLowerCase().equals("mean")) {
                        options.summaryStrategy = SummaryStrategy.MEAN;

                        i += 1;
                        break;
                    }

                    if (args[i+1].toLowerCase().equals("median")) {
                        options.summaryStrategy = SummaryStrategy.MEDIAN;

                        i += 1;
                        break;
                    }

                    printUsageAndError("-positions must be followed by either 'MEAN' or 'MEDIAN'.");

                case "-threshold":
                    if (args.length<=i+1) {
                        printUsageAndError("-threshold must be followed by a number.");
                    }

                    try {
                        options.convSupportThresh =
                                Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Threshold must be a number between 0 and 100.");
                    }

                    i += 1;
                    break;

                case "-recordGeneFlow":
                    if (args.length<=i+1) {
                        printUsageAndError("-recordGeneFlow must be followed by a file name.");
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
            new ReassortmentDistance(options);

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