package coalre.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.math.Binomial;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;

import java.util.BitSet;

public abstract class NetworkOperator extends Operator {

    public Input<Network> networkInput = new Input<>("network",
            "Network on which to operate",
            Input.Validate.REQUIRED);

    public static final double LOG2 = Math.log(2.0);

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
     * @return the segment subset.
     */
    protected BitSet getRandomConditionedSubset(BitSet sourceSegments) {

        if (sourceSegments.cardinality()<2) {
            throw new RuntimeException("Too few segments to draw conditioned subset.");
        }

        BitSet destSegments = new BitSet();

        int firstSegPos = Randomizer.nextInt(sourceSegments.cardinality());
        int lastSegPos;
        do {
            lastSegPos = Randomizer.nextInt(sourceSegments.cardinality());
        } while (lastSegPos == firstSegPos);

        int segPos = 0;
        for (int segIdx=sourceSegments.nextSetBit(0); segIdx >=0;
             segIdx = sourceSegments.nextSetBit(segIdx+1), segPos += 1) {

            if (segPos == firstSegPos) {
                destSegments.set(segIdx);
                continue;
            }

            if (segPos == lastSegPos)
                continue;

            if (Randomizer.nextBoolean()) {
                destSegments.set(segIdx);
            }

        }

        return destSegments;
    }

    /**
     * Compute the probability of choosing a particular subset of source segments
     * using the getRandomConditionedSubset method.
     *
     * @param sourceSegments set of segments used as the argument to getRandomConditionedSubset
     * @return log probability of subset
     */
    protected double getConditionedSubsetProb(BitSet sourceSegments) {
        return Math.log(1.0/sourceSegments.cardinality()) +
                Math.log(1.0/(sourceSegments.cardinality()-1)) +
                (sourceSegments.cardinality()-2)*Math.log(0.5);
    }
}
