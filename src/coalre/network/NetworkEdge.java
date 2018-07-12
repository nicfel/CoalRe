package coalre.network;

import java.util.*;

public class NetworkEdge {

    public NetworkNode parentNode, childNode;
    public BitSet hasSegments;

    public NetworkEdge() { }

    public NetworkEdge(NetworkNode parentNode, NetworkNode childNode,
                       BitSet hasSegments) {
        this.parentNode = parentNode;
        this.childNode = childNode;

        this.hasSegments = hasSegments;
    }

    public double getReassortmentObsProb() {
        // There are always two reassortment configurations that
        // produce an unobserved reassortment: 1111 and 0000
        // (assuming 4 segs on lineage)
        return 1.0 - 2.0*Math.pow(0.5, hasSegments.cardinality());
    }

    public double getLength() {
        return parentNode.getHeight() - childNode.getHeight();
    }

    public boolean isRootEdge() {
        return parentNode == null;
    }

}
