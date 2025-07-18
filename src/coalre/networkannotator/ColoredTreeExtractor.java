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

import beast.base.core.BEASTObject;
import beast.base.core.Function;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Maps reassortment events onto segments
 * 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class ColoredTreeExtractor extends ReassortmentAnnotator {

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("summary.tree");
        double burninPercentage = 0;
        int[] removeSegments = new int[0];
        int outputSegment = -1;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "Outputs segment: " + outputSegment;
        }
    }

    public ColoredTreeExtractor(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader

        ReassortmentLogReader logReader = new ReassortmentLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

        // keeps track of the leave nodes
        List<String> leafNodes = new ArrayList<>();
               
        PrintStream ps = new PrintStream(options.outFile);
        
        
        // build the clades
        boolean first = true;
        int c= 1;
        for (Network network : logReader){
        	if (first){        		
            	for (NetworkNode networkNode : network.getNodes()){
            		if (networkNode.isLeaf()){
            			leafNodes.add(networkNode.getTaxonLabel());
            		}
        		}
        	}
        	pruneNetwork(network, options.removeSegments);

        	// get the tree with single child nodes back
        	Tree tree = getSingleChildTree(network, options.outputSegment, leafNodes);
        	
        	if (first){
        		tree.init(ps);
        		ps.print("\n");
        		first = false;
        	}
        	
        	
        	printTree(tree, ps, c);
        	c++;
        	
        	ps.print("\n");
        }
   		ps.println("End;");
   		ps.close();
        System.out.println("\nDone!");

    }    
    
    static Tree getSingleChildTree(Network network, int segment, List<String> leafNodes){
    	// get teh root of this segment tree
        List<NetworkNode> rootEdge = network.getNodes().stream()
                .filter(e -> e.isCoalescence())
                .filter(e -> !e.getParentEdges().get(0).hasSegments.get(segment))
                .filter(e -> e.getChildEdges().get(0).hasSegments.get(segment))
                .filter(e -> e.getChildEdges().get(1).hasSegments.get(segment))
                .collect(Collectors.toList());
        
        if (rootEdge.size()!=1)
        	throw new IllegalArgumentException("root of segment tree not found");

    	Node root = getNextNode(rootEdge.get(0), network.getSegmentCount(), segment, leafNodes);
    	   	
    	Tree tree = new Tree(root);
    	return tree;
    }
    
    private static Node getNextNode(NetworkNode networkNode, int nrSegments, int segment, List<String> leafNodes){
    	
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
    
    private void printTree(Tree tree, PrintStream out, int sample){
	    out.print("tree STATE_" + sample + " = ");
	    // Don't sort, this can confuse CalculationNodes relying on the tree
	    //tree.getRoot().sort();
	    final int[] dummy = new int[1];
	    final String newick = tree.getRoot().toSortedNewick(dummy, true);
	    out.print(newick);
	    out.print(";");
    }
    
    
    
    

//    /**
//     * Use a GUI to retrieve ACGAnnotator options.
//     *
//     * @param options options object to populate using GUI
//     * @return true if options successfully collected, false otherwise
//     */
//    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {
//
//        boolean[] canceled = {false};
//
//        JDialog dialog = new JDialog((JDialog)null, true);
//        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
//        dialog.setLocationRelativeTo(null);
//        dialog.setTitle("Reassortment Network Annotator");
//
//        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
//        JLabel outFileLabel = new JLabel("Output file:");
//        JLabel burninLabel = new JLabel("Burn-in percentage:");
//        JLabel summaryMethodLabel = new JLabel("Position summary method:");
//        JLabel thresholdLabel = new JLabel("Posterior conversion support threshold:");
//        JCheckBox geneFlowCheckBox = new JCheckBox("Record gene flow");
//
//        JTextField inFilename = new JTextField(20);
//        inFilename.setEditable(false);
//        JButton inFileButton = new JButton("Choose File");
//
//        JTextField outFilename = new JTextField(20);
//        outFilename.setText(options.outFile.getName());
//        outFilename.setEditable(false);
//        JButton outFileButton = new JButton("Choose File");
//
//        JTextField gfOutFilename = new JTextField(20);
//        gfOutFilename.setEditable(false);
//        gfOutFilename.setEnabled(false);
//        JButton gfOutFileButton = new JButton("Choose File");
//        gfOutFileButton.setEnabled(false);
//
//        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
//                0, 100, (int)(options.burninPercentage));
//        burninSlider.setMajorTickSpacing(50);
//        burninSlider.setMinorTickSpacing(10);
//        burninSlider.setPaintTicks(true);
//        burninSlider.setPaintLabels(true);
//        burninSlider.setSnapToTicks(true);
//
//        JComboBox<SummaryStrategy> heightMethodCombo = new JComboBox<>(SummaryStrategy.values());
//
//        JSlider thresholdSlider = new JSlider(JSlider.HORIZONTAL,
//                0, 100, (int)(options.convSupportThresh));
//        thresholdSlider.setMajorTickSpacing(50);
//        thresholdSlider.setMinorTickSpacing(10);
//        thresholdSlider.setPaintTicks(true);
//        thresholdSlider.setPaintLabels(true);
//        thresholdSlider.setSnapToTicks(true);
//
//        Container cp = dialog.getContentPane();
//        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
//        cp.setLayout(boxLayout);
//
//        JPanel mainPanel = new JPanel();
//
//        GroupLayout layout = new GroupLayout(mainPanel);
//        mainPanel.setLayout(layout);
//        layout.setAutoCreateGaps(true);
//        layout.setAutoCreateContainerGaps(true);
//
//        layout.setHorizontalGroup(layout.createSequentialGroup()
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(logFileLabel)
//                        .addComponent(outFileLabel)
//                        .addComponent(burninLabel)
//                        .addComponent(summaryMethodLabel)
//                        .addComponent(thresholdLabel)
//                        .addComponent(geneFlowCheckBox))
//                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
//                        .addComponent(inFilename)
//                        .addComponent(outFilename)
//                        .addComponent(burninSlider)
//                        .addComponent(heightMethodCombo)
//                        .addComponent(thresholdSlider)
//                        .addComponent(gfOutFilename))
//                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
//                        .addComponent(inFileButton)
//                        .addComponent(outFileButton)
//                        .addComponent(gfOutFileButton)));
//
//        layout.setVerticalGroup(layout.createSequentialGroup()
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(logFileLabel)
//                        .addComponent(inFilename,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE)
//                        .addComponent(inFileButton))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(outFileLabel)
//                        .addComponent(outFilename,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE)
//                        .addComponent(outFileButton))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(burninLabel)
//                        .addComponent(burninSlider,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(summaryMethodLabel)
//                        .addComponent(heightMethodCombo,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(thresholdLabel)
//                        .addComponent(thresholdSlider,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(geneFlowCheckBox)
//                        .addComponent(gfOutFilename)
//                        .addComponent(gfOutFileButton)));
//
//        mainPanel.setBorder(new EtchedBorder());
//        cp.add(mainPanel);
//
//        JPanel buttonPanel = new JPanel();
//
//        JButton runButton = new JButton("Analyze");
//        runButton.addActionListener((e) -> {
//            options.burninPercentage = burninSlider.getValue();
//            options.convSupportThresh = thresholdSlider.getValue();
//            options.summaryStrategy = (SummaryStrategy)heightMethodCombo.getSelectedItem();
//            dialog.setVisible(false);
//        });
//        runButton.setEnabled(false);
//        buttonPanel.add(runButton);
//
//        JButton cancelButton = new JButton("Quit");
//        cancelButton.addActionListener((e) -> {
//            dialog.setVisible(false);
//            canceled[0] = true;
//        });
//        buttonPanel.add(cancelButton);
//
//        JFileChooser inFileChooser = new JFileChooser();
//        inFileButton.addActionListener(e -> {
//            inFileChooser.setDialogTitle("Select ACG log file to summarize");
//            if (options.inFile == null)
//                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//            int returnVal = inFileChooser.showOpenDialog(dialog);
//
//            if (returnVal == JFileChooser.APPROVE_OPTION) {
//                options.inFile = inFileChooser.getSelectedFile();
//                inFilename.setText(inFileChooser.getSelectedFile().getName());
//                runButton.setEnabled(true);
//            }
//        });
//
//        JFileChooser outFileChooser = new JFileChooser();
//        outFileButton.addActionListener(e -> {
//            outFileChooser.setDialogTitle("Select output file name.");
//            if (options.inFile != null)
//                outFileChooser.setCurrentDirectory(options.inFile);
//            else
//                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//
//            outFileChooser.setSelectedFile(options.outFile);
//            int returnVal = outFileChooser.showOpenDialog(dialog);
//
//            if (returnVal == JFileChooser.APPROVE_OPTION) {
//                options.outFile = outFileChooser.getSelectedFile();
//                outFilename.setText(outFileChooser.getSelectedFile().getName());
//            }
//        });
//
//        geneFlowCheckBox.addActionListener(e -> {
//            boolean newValue = geneFlowCheckBox.isSelected();
//            gfOutFilename.setEnabled(newValue);
//            gfOutFileButton.setEnabled(newValue);
//        });
//
//        JFileChooser gfOutFileChooser = new JFileChooser();
//        gfOutFileButton.addActionListener(e -> {
//            gfOutFileChooser.setDialogTitle("Select gene flow output file name.");
//            if (options.inFile != null)
//                gfOutFileChooser.setCurrentDirectory(options.inFile);
//            else
//                gfOutFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//
//            int returnVal = gfOutFileChooser.showOpenDialog(dialog);
//
//            if (returnVal == JFileChooser.APPROVE_OPTION) {
//                gfOutFilename.setText(gfOutFileChooser.getSelectedFile().getName());
//            }
//        });
//
//        cp.add(buttonPanel);
//
//        dialog.pack();
//        dialog.setResizable(false);
//        dialog.setVisible(true);
//
//        return !canceled[0];
//    }
//
//    /**
//     * Prepare JFrame to which ACGAnnotator output streams will be
//     * directed.
//     */
//    private static void setupGUIOutput() {
//
//        JFrame frame = new JFrame();
//        frame.setTitle("ACGAnnotator");
//        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
//
//        JTextArea textArea = new JTextArea(25, 80);
//        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
//        textArea.setEditable(false);
//        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);
//
//        JButton closeButton = new JButton("Close");
//        closeButton.addActionListener(e -> System.exit(0));
//        JPanel buttonPanel = new JPanel();
//        buttonPanel.add(closeButton);
//        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);
//
//        // Redirect streams to output window:
//        OutputStream out = new OutputStream() {
//            @Override
//            public void write(int b) throws IOException {
//                SwingUtilities.invokeLater(() -> {
//                    if ((char)b == '\r') {
//                        int from = textArea.getText().lastIndexOf("\n") + 1;
//                        int to = textArea.getText().length();
//                        textArea.replaceRange(null, from, to);
//                    } else
//                        textArea.append(String.valueOf((char) b));
//                });
//            }
//        };
//
//        System.setOut(new PrintStream(out, true));
//        System.setErr(new PrintStream(out, true));
//
//        frame.pack();
//        frame.setVisible(true);
//    }

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
                    
                case "-outputSegment":
                    if (args.length<=i+1) {
                        printUsageAndError("-outputSegment must be followed by exactly one number.");
                    }

                    try {
                		options.outputSegment = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("outputSegment must be an integer");
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

//        if (args.length == 0) {
//            // Retrieve options from GUI:
//
//            try {
//                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
//            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
//                Log.warning.println("Error setting cross-platform look and feel.");
//            }
//
//            try {
//                SwingUtilities.invokeAndWait(() -> {
//                    if (!getOptionsGUI(options))
//                        System.exit(0);
//
//                    setupGUIOutput();
//                });
//            } catch (InterruptedException | InvocationTargetException e) {
//                e.printStackTrace();
//            }
//
//
//        } else {
            getCLIOptions(args, options);
//        }

        // Run ACGAnnotator
        try {
            new ColoredTreeExtractor(options);

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