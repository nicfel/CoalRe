package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public abstract class NetworkOperator extends Operator {

    public Input<Network> networkInput = new Input<>("network",
            "Network on which to operate",
            Input.Validate.REQUIRED);

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree associated with network.",
            new ArrayList<>());

    protected Network network;
    protected List<Tree> segmentTrees;
    int[] segmentTreeMap;
    
	protected BitSet segmentsChanged = new BitSet();
	
	


    @Override
    public void initAndValidate() {
        network = networkInput.get();
        segmentTrees = segmentTreesInput.get();

        if (segmentTrees.isEmpty()) {
            Log.warning.println("Warning: no segment trees provided to " +
                    "network operator " + getID() + ".");
        } else {
            if (segmentTrees.size() != network.getSegmentCount())
                throw new IllegalArgumentException("Number of segment trees provided " +
                        "to operator must match number of segments in network.");
            // ensure the order of the segment trees is correct
            segmentTreeMap = new int[network.getSegmentCount()];
            
            if (network.segmentNames==null) {
            	network.segmentNames=new String[segmentTreesInput.get().size()];
            	for (int j = 0; j < segmentTreesInput.get().size(); j++) {
            		network.segmentNames[j] = segmentTreesInput.get().get(j).getID();
            	}
            }
            
            int c=0;
            for (int i = 0; i < network.getSegmentCount(); i++) {
            	for (int j = 0; j < segmentTreesInput.get().size(); j++) {
            		if (segmentTreesInput.get().get(j).getID().contentEquals(network.segmentNames[i])) {
            			segmentTreeMap[i] = j;
            			c++;
            		}
            	}
            }
            if (c!=network.getSegmentCount()) {
            	throw new IllegalArgumentException("not all segments found");
            }         
            

        }

    }

//    static int count = 0;

    public double proposal() {
		if (network.segmentTreesNeedInit) {			
			System.out.println("restore segment trees in first move");
			// update all segment trees			
            for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
          		network.updateSegmentTree(segmentTrees.get(segmentTreeMap[segIdx]), segmentTreeMap[segIdx]); 
            
            network.segmentTreesNeedInit = false;
			return Double.POSITIVE_INFINITY;
		}


//        count += 1;
//
//        System.out.println("count = " + count);
    	
//        System.out.println(network);
        
//        for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
//            System.out.println(segmentTrees.get(segIdx) +";");
//        System.out.println("=========");

//        String[] trees = new String[segmentTrees.size()];
//        if (this instanceof TipReheight) {
//            if (((TipReheight) this).taxonsetInput.get().getTaxonId(0).contentEquals("RVA/1070/Thailand|2013-01-01")) {
//                System.out.println(network);
//                for (int i =0; i < segmentTrees.size(); i++) {
//                    trees[i] = segmentTrees.get(i) + ";";
//                }
//            }
//        }
    	// set segmentsChanged to true
//        segmentsChanged.set(0, network.getSegmentCount(), true);

        double logHR = 0.0;
//        try {
        	logHR += networkProposal();
//        } catch (Exception e) {
//            Log.warning.println("re-initialize segment tree map");
//            e.printStackTrace();
//            // update segment trees and return Double.NEGATIVE_INFINITY
//            for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
//          		network.updateSegmentTree(segmentTrees.get(segmentTreeMap[segIdx]), segmentTreeMap[segIdx]);            
//        }

//        if (logHR>Double.NEGATIVE_INFINITY) {
//            for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
//            	if (segmentsChanged.get(segIdx))
//            		network.updateSegmentTree(segmentTrees.get(segmentTreeMap[segIdx]), segmentTreeMap[segIdx]);
//        }



//        if (this instanceof TipReheight) {
//            if (((TipReheight) this).taxonsetInput.get().getTaxonId(0).contentEquals("RVA/1070/Thailand|2013-01-01")) {
//                System.out.println("\n\n\n\n\n\n\n\n\n\n\n\n\n");
//                for (int i =0; i < segmentTrees.size(); i++) {
//                    System.out.println(segmentTrees.get(i).getID());
//                    System.out.println(trees[i]);
//                    System.out.println(segmentTrees.get(i) + ";");
//                }
//                System.out.println(network);
//
//                System.out.println("\n=========");
//            }
//        }

        
//        System.out.println(getID());
//        System.out.println(network);
        
//        System.out.println(network);
//        System.out.println(".............");
//        
//        for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
//            System.out.println(segmentTrees.get(segIdx) +";");


        return logHR;
    }

    /**
     * Propose a new network state.  The network state will be
     * used to update the segment trees once the network proposal
     * is complete.
     *
     * @return log of HR for proposal.
     */
    protected abstract double networkProposal();

    /**
     * Retrieve sister of given edge
     * @param childEdge child edge
     * @return sister of given child edge
     */
    protected NetworkEdge getSisterEdge(NetworkEdge childEdge) {
        int idx = childEdge.parentNode.getChildEdges().indexOf(childEdge);
        int otherIdx = (idx + 1) % 2;

        return childEdge.parentNode.getChildEdges().get(otherIdx);
    }

    /**
     * Retrieve spouse of given edge
     * @param parentEdge parent edge
     * @return spouse of given parent edge
     */
    protected NetworkEdge getSpouseEdge(NetworkEdge parentEdge) {
        int idx = parentEdge.childNode.getParentEdges().indexOf(parentEdge);
        int otherIdx = (idx + 1) % 2;

        return parentEdge.childNode.getParentEdges().get(otherIdx);
    }

    /**
     * Chooses a number of segments from a subset of the provided source segments.
     * A first segment is chosen uniformly at random will definitely be included in the subset.
     * Another segment chosen uniformly at random will not be included in the subset.
     * Each of the remaining segments are included in the subset with probability 0.5.
     *
     * @param sourceSegments set of segments from which to choose
     * @return the segment subset or null if conditional subsetting is impossible.
     */
    protected BitSet getRandomConditionedSubset(BitSet sourceSegments) {

        if (sourceSegments.cardinality()<2)
            return null;

        BitSet destSegments = new BitSet();

        do {
            destSegments.clear();

            for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1;
                 segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

                if (Randomizer.nextBoolean())
                    destSegments.set(segIdx);
            }

        } while (destSegments.cardinality() == 0
                || destSegments.cardinality() == sourceSegments.cardinality());

        return destSegments;
    }

    /**
     * Compute the probability of choosing a particular subset of source segments
     * using the getRandomConditionedSubset method.
     *
     * @param sourceSegments set of segments used as the argument to getRandomConditionedSubset
     * @return log probability of subset
     */
    protected double getLogConditionedSubsetProb(BitSet sourceSegments) {

        if (sourceSegments.cardinality()<2)
            return Double.NEGATIVE_INFINITY;

        return sourceSegments.cardinality()*Math.log(0.5)
                - Math.log(1.0 - 2.0*Math.pow(0.5, sourceSegments.cardinality()));
    }
    
    public void replace(final NetworkNode node, final NetworkEdge child, final NetworkEdge replacement) {
    	node.removeChildEdge(child);
    	node.addChildEdge(replacement);
    }

    /**
     * Will be used in the next version of BEAST to prevent date trait cloning
     * from breaking the BEAuti model.
     */
    public boolean notCloneable() {
        return true;
    }
}
