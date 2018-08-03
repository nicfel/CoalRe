package coalre.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;

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

    final public double proposal() {
        double logHR = networkProposal();

        if (logHR>Double.NEGATIVE_INFINITY) {
            for (Tree segmentTree : segmentTreesInput.get())
                updateSegmentTree(segmentTree);
        }

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
     * Chooses a number of segments form a subset of the provided source segments.
     * A first segment is chosen uniformly at random will definitely be included in the subset.
     * Another segment chosen uniformly at random will not be included in the subset.
     * Each of the remaining segments are included in the subset with probability 0.5.
     *
     * @param sourceSegments set of segments from which to choose
     * @return the segment subset or null if conditional subsetting is impossible.
     */
    protected BitSet getRandomConditionedSubset(BitSet sourceSegments) {

        if (sourceSegments.cardinality()<2) {
            return null;
        }

        BitSet destSegments = new BitSet();

        do {

            destSegments.clear();

            for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1;
                 segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

                if (Randomizer.nextBoolean()) {
                    destSegments.set(segIdx);
                }

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

    /**
     * Modify segmentTree so that it matches the
     * @param segmentTree
     */
    protected void updateSegmentTree(Tree segmentTree) {

    }
}
