package coalre.operators;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

public class NetworkExchange extends DivertSegmentOperator {
	final public Input<Boolean> isNarrowInput = new Input<>("isNarrow", "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);
	
	//TODO: What is this doing??
//    @Override
//    public void initAndValidate() {
//    }

    @Override
    public double networkProposal() {
        double logHR = 0.0;
        
        final Network network = networkInput.get(this);
        
        if (isNarrowInput.get()) {
        	logHR = narrow(network);
        } else {
            logHR = wide(network);
        }

        return logHR;
    }
    
    
    public double narrow(final Network network) {
    	
//    	List<NetworkNode> internalNodes = new ArrayList<>(network.getInternalNodes());
//    	final int internalNodesCount = internalNodes.size();
//        if (internalNodesCount <= 1) {
//            return Double.NEGATIVE_INFINITY;
//        }
        
    	System.out.println("yo");
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
        
        List<NetworkEdge> possibleGrandParentEdges = networkEdges.stream()
        		.filter(e -> !e.isLeafEdge())
        		.filter(e -> !e.isRootEdge())
        		.filter(e -> e.childNode.isCoalescence()) // good idea?
        		.collect(Collectors.toList());
        
        NetworkEdge grandParentEdge = possibleGrandParentEdges.
        		get(Randomizer.nextInt(possibleGrandParentEdges.size()));
        
        NetworkNode grandParent = grandParentEdge.childNode;
        
        
        
        List<NetworkEdge> possibleParentEdges = grandParent.getChildEdges(); // 2
//        		.stream()
//        		.filter(e -> !e.isLeafEdge()).collect(Collectors.toList());
        
        //TODO: Would this work?
        final int parentId = Randomizer.nextInt(possibleParentEdges.size());
        NetworkEdge parentEdge = possibleParentEdges.get(parentId);
        NetworkEdge auntEdge = super.getSisterEdge(parentEdge);
        
        NetworkNode parent = parentEdge.childNode;
        NetworkNode aunt = auntEdge.childNode;
        
//        System.out.println(possibleParentEdges.size());
//        final int[] ids = new Random().ints(0, possibleParentEdges.size())
//        		.distinct().limit(2).toArray();
        
        
//        NetworkEdge parentEdge = possibleParentEdges.get(ids[0]);
//        NetworkEdge auntEdge = possibleParentEdges.get(ids[1]);
//        
//        
//        NetworkNode parent = parentEdge.childNode;
//        NetworkNode aunt = auntEdge.childNode;
        
        if (parent.getHeight() < aunt.getHeight()) {
        	auntEdge = possibleParentEdges.get(parentId);
        	parentEdge = super.getSisterEdge(auntEdge);
        	
        	parent = parentEdge.childNode;
        	aunt = auntEdge.childNode;
        } 
        
        if( parent.isLeaf() ) {
        	System.out.println("Doooesn't count");
            // tree with dated tips
            return Double.NEGATIVE_INFINITY;
        }
        
        // Hastings ratio not updated for now

        List<NetworkEdge> possibleChildEdges = parent.getChildEdges();
        final NetworkEdge childEdge = possibleChildEdges.get(Randomizer.nextInt(possibleChildEdges.size()));
        exchangeEdges(childEdge, auntEdge, parent, grandParent);
        
        BitSet childSegs = childEdge.hasSegments;
        BitSet auntSegs = auntEdge.hasSegments;
        
        removeSegmentsFromAncestors(grandParentEdge, auntSegs);
        removeSegmentsFromAncestors(parentEdge, childSegs);
        
        addSegmentsToAncestors(grandParentEdge, childSegs);
        addSegmentsToAncestors(parentEdge, auntSegs);
        

        
        
 
        
    	
    	return 1;
    }
    
    public double wide(final Network network) {
    	
    	//TODO: Everything
    	
    	return 0.0;
    }
    
    
    /* exchange sub-nets whose root are i and j */

    protected void exchangeEdges(NetworkEdge i, NetworkEdge j,
    		NetworkNode p, NetworkNode jP) {
        // precondition p -> i & jP -> j
//        replace(p, i, j);
    	System.out.println("Exchange happening");
        p.removeChildEdge(i);
        jP.removeChildEdge(j);
        p.addChildEdge(j);
        jP.addChildEdge(i);
//        replace(jP, j, i);
        // postcondition p -> j & p -> i
    }
    
    
    
    
    
       /**
     * Check that each edge is ancestral to at least one segment.
     *
     * @return true if all edges are ancestral.
     */
    public boolean allEdgesAncestral() {
        Set<NetworkNode> nodeList = networkInput.get().getNodes();
        for (NetworkNode node : nodeList) {
            for (NetworkEdge parentEdge : node.getParentEdges()) {
                if (parentEdge.hasSegments.isEmpty())
                    return false;
            }
        }

        return true;
    }

}
