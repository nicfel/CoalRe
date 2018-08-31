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

    public NetworkEdge getCopy() {
        return getCopy(new HashMap<>());
    }

    public NetworkEdge getCopy(Map<NetworkNode,NetworkNode> seenNodes) {

        NetworkEdge edgeCopy = new NetworkEdge(null, null, (BitSet)hasSegments.clone());
        NetworkNode childNodeCopy;
        boolean traverse = true;
        if (seenNodes.containsKey(childNode)) {
            childNodeCopy = seenNodes.get(childNode);
            traverse = false;
        } else {
            childNodeCopy = new NetworkNode();
            childNodeCopy.setHeight(childNode.getHeight());
            childNodeCopy.setTaxonLabel(childNode.getTaxonLabel());
            childNodeCopy.setTaxonIndex(childNode.getTaxonIndex());
            seenNodes.put(childNode, childNodeCopy);
        }

        childNodeCopy.addParentEdge(edgeCopy);

        if (traverse) {
            for (NetworkEdge childEdge : childNode.getChildEdges()) {
                NetworkEdge childEdgeCopy = childEdge.getCopy(seenNodes);
                childNodeCopy.addChildEdge(childEdgeCopy);
            }
        }

        return edgeCopy;
    }

}
