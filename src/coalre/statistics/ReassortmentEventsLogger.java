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

public class ReassortmentEventsLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree associated with network.",
            new ArrayList<>());


    Network network;
    int segCount;
    boolean logObservable = false;


    public ReassortmentEventsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        segCount = network.getSegmentCount();
        
        if(segmentTreesInput.get().size()>0)
        	logObservable = true;

    }

    @Override
    public void init(PrintStream out) {	
        out.print("beforeEvent:parentLineage1:parentlineage2\t");        		

    }

    @Override
    public void log(long sample, PrintStream out) {
    	
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
    		
//    		BitSet childSegment = new BitSet();
//    		BitSet parent1Segments = new BitSet();
//    		BitSet parent2Segments = new BitSet();
//    		childSegment = (BitSet) childEdge.get(0).hasSegments.clone();
//    		childSegment = (BitSet) childEdge.get(0).hasSegments.clone();
//    		childSegment = (BitSet) childEdge.get(0).hasSegments.clone();
//    		
    		final BitSet childSegment = (BitSet) childEdge.get(0).hasSegments.clone();
    		final BitSet parent1Segments = (BitSet) parentEdges.get(0).hasSegments.clone();
    		final BitSet parent2Segments = (BitSet) parentEdges.get(1).hasSegments.clone();
    		
    		int beforeEvent = 0, parentLineage1 = 0, parentLineage2 = 0;
    		
    		// count which edges go where
            for (int i = 0; i < segCount; i++){            	
            	if (reassortmentNodes.get(k).getHeight()>rootHeights[i]){
            		childSegment.set(i, false);
            		parent1Segments.set(i, false);
            		parent2Segments.set(i, false);
            	}
            } 
            if (childSegment.cardinality()>1)
            	out.print(childSegment + ":" + parent1Segments + ":" + parent2Segments + ":" + reassortmentNodes.get(k).getHeight()  + "\t");
    	}
    	
	}

    @Override
    public void close(PrintStream out) {

    }
}
