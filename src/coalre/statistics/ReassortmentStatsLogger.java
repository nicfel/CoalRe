package coalre.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class ReassortmentStatsLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree associated with network.",
            new ArrayList<>());


    Network network;
    int segCount;
    boolean logObservable = false;


    public ReassortmentStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        segCount = network.getSegmentCount();
        
        if(segmentTreesInput.get().size()>0)
        	logObservable = true;

    }

    @Override
    public void init(PrintStream out) {
        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

		if (logObservable){
	        for (int i = 0; i < segCount; i++){
	        	for (int j = i+1; j < segCount; j++){
	                out.print(prefix + "jointReassortment." + i + "_" + j + "\t");        		
	                out.print(prefix + "splitReassortment." + i + "_" + j + "\t");        		
	        	}
	        }    
		}else{
	        for (int i = 0; i < segCount; i++){
	        	for (int j = i+1; j < segCount; j++){
	                out.print(prefix + "jointObsReassortment." + i + "_" + j + "\t");        		
	                out.print(prefix + "splitObsReassortment." + i + "_" + j + "\t");        		
	        	}
	        }    

		}
		
    }

    @Override
    public void log(long sample, PrintStream out) {
    	int[][] jointReassortmentCount = new int[segCount][segCount];
    	int[][] splitReassortmentCount = new int[segCount][segCount];
    	
		double[] rootHeights = new double[segCount];
		if (logObservable)
			for (int i = 0; i < segCount; i++)
				rootHeights[i] = segmentTreesInput.get().get(i).getRoot().getHeight();
		else
			for (int i = 0; i < segCount; i++)
				rootHeights[i] = Double.POSITIVE_INFINITY;


    	
		List<NetworkNode> networkNodes = new ArrayList<>(network.getNodes());

    	final List<NetworkNode> reassortmentNodes = networkNodes.stream()
    			.filter(NetworkNode::isReassortment)
    			.collect(Collectors.toList());
    	
    	for (int k = 0; k < reassortmentNodes.size(); k++){
    		List<NetworkEdge> parentEdges = reassortmentNodes.get(k).getParentEdges();
    		List<NetworkEdge> childEdge = reassortmentNodes.get(k).getChildEdges();
    		
    		final BitSet childSegment = childEdge.get(0).hasSegments;
    		final BitSet parent1Segments = parentEdges.get(0).hasSegments;
    		final BitSet parent2Segments = parentEdges.get(1).hasSegments;
    		
    		// count which edges go where
            for (int i = 0; i < segCount; i++){
            	if (childSegment.get(i) && reassortmentNodes.get(k).getHeight()<rootHeights[i]){
	            	for (int j = i+1; j < segCount; j++){
	            		if(childSegment.get(j) && reassortmentNodes.get(k).getHeight()<rootHeights[j]){
	            			if (parent1Segments.get(i) && parent1Segments.get(j)){
	            				jointReassortmentCount[i][j]++;
	            			}else if(parent2Segments.get(i) && parent2Segments.get(j)){
	            				jointReassortmentCount[i][j]++;
	            			}else{
	            				splitReassortmentCount[i][j]++;
	            			}	            			
	            		}
	            	}
            	}
            } 
    		
    	}
    	
        for (int i = 0; i < segCount; i++){
        	for (int j = i+1; j < segCount; j++){
                out.print(jointReassortmentCount[i][j]+ "\t" +
                		splitReassortmentCount[i][j]+ "\t");
        	}
        } 
	}

    @Override
    public void close(PrintStream out) {

    }
}
