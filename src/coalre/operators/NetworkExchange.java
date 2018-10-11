package coalre.operators;

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

public class NetworkExchange extends NetworkOperator {

    @Override
    public double networkProposal() {
        double logHR = 0.0;

        List<NetworkEdge> sourceEdges = network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>1)
                .collect(Collectors.toList());

        if (sourceEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceEdges.size());

        NetworkEdge sourceEdge = sourceEdges.get(Randomizer.nextInt(sourceEdges.size()));
        NetworkEdge destEdge = getSpouseEdge(sourceEdge);

        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
        logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);

        network.startEditing(this);

        logHR -= addSegmentsToAncestors(destEdge, segsToDivert);
        logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert);

        if (!allEdgesAncestral())
            return Double.NEGATIVE_INFINITY;

        logHR += getLogConditionedSubsetProb(destEdge.hasSegments);

        int reverseSourceEdgeCount = (int)(network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>1)
                .count());

        logHR += Math.log(1.0/reverseSourceEdgeCount);

        return logHR;
    }
    
    
    public double narrow() {
    	
//    	List<NetworkNode> internalNodes = new ArrayList<>(network.getInternalNodes());
//    	final int internalNodesCount = internalNodes.size();
//        if (internalNodesCount <= 1) {
//            return Double.NEGATIVE_INFINITY;
//        }
        
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
        
        List<NetworkEdge> possibleGrandParentEdges = networkEdges.stream()
        		.filter(e -> !e.isLeafEdge()).collect(Collectors.toList());
        
        NetworkNode grandParent = possibleGrandParentEdges.
        		get(Randomizer.nextInt(possibleGrandParentEdges.size())).parentNode;
        
        
        
        List<NetworkEdge> possibleParentEdges = grandParent.getChildEdges();
        
        //TODO: Would this work?
//        NetworkEdge parentEdge = possibleParentEdges.get(Randomizer.nextInt(possibleParentEdges.size()));
//        NetworkEdge auntEdge = super.getSisterEdge(parentEdge);
//        
//        NetworkNode parent = parentEdge.childNode;
//        NetworkNode aunt = auntEdge.childNode;
        
        final int[] ids = new Random().ints(possibleParentEdges.size())
        		.distinct().limit(2).toArray();
        
        NetworkNode parent = possibleParentEdges.get(ids[0]).childNode;
        NetworkNode uncle = possibleParentEdges.get(ids[1]).childNode;
        
        if (parent.getHeight() < uncle.getHeight()) {
        	parent = possibleParentEdges.get(ids[1]).childNode;
        	uncle = possibleParentEdges.get(ids[0]).childNode;
        }
        
        // Hastings ratio not updated for now
 
        
    	
    	return 0;
    }
    
    
    /* exchange sub-nets whose root are i and j */

    protected void exchangeNodes(NetworkEdge i, NetworkEdge j,
    		NetworkNode p, NetworkNode jP) {
        // precondition p -> i & jP -> j
        replace(p, i, j);
        replace(jP, j, i);
        // postcondition p -> j & p -> i
    }
    
    
    
    
    
    


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
        double logP = 0.0;

        segsToRemove = (BitSet)segsToRemove.clone();
        segsToRemove.and(edge.hasSegments);

        if (segsToRemove.isEmpty())
            return logP;

        edge.hasSegments.andNot(segsToRemove);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isReassortment()) {

            logP += Math.log(0.5)*segsToRemove.cardinality();

            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), segsToRemove);

        } else {

            segsToRemove.andNot(getSisterEdge(edge).hasSegments);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);

        }

        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd) {
        double logP = 0.0;

        segsToAdd = (BitSet)segsToAdd.clone();
        segsToAdd.andNot(edge.hasSegments);

        if (segsToAdd.isEmpty())
            return logP;

        edge.hasSegments.or(segsToAdd);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isReassortment()) {

            BitSet segsToAddLeft = new BitSet();
            BitSet segsToAddRight = new BitSet();

            for (int segIdx=segsToAdd.nextSetBit(0); segIdx != -1;
                    segIdx=segsToAdd.nextSetBit(segIdx+1)) {
                if (Randomizer.nextBoolean())
                    segsToAddLeft.set(segIdx);
                else
                    segsToAddRight.set(segIdx);

                logP += Math.log(0.5);
            }

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);

        } else {

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
        }

        return logP;
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
