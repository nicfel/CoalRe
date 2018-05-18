package coalre.operators;

import beast.core.Input;
import beast.core.Operator;
import coalre.network.Network;
import coalre.network.NetworkEdge;

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
}
