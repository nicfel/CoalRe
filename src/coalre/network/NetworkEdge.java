package coalre.network;

import java.util.*;

public class NetworkEdge {

    public NetworkNode parentNode, childNode;
    public BitSet hasSegments;

    public NetworkEdge(NetworkNode parentNode, NetworkNode childNode,
                       BitSet hasSegments) {
        this.parentNode = parentNode;
        this.childNode = childNode;

        this.hasSegments = hasSegments;
    }

    public NetworkNode getParentNode() {
        return parentNode;
    }

    public NetworkNode getChildNode() {
        return childNode;
    }

    public void setParentNode(NetworkNode newParentNode) {
        parentNode = newParentNode;
    }

    public void setChildNode(NetworkNode newChildNode) {
        childNode = newChildNode;
    }

    public BitSet getHasSegments() {
        return hasSegments;
    }

    public double getReassortmentObsProb() {
        return 1.0 - Math.pow(0.5, hasSegments.cardinality()-1);
    }

}
