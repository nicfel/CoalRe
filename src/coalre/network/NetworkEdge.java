package coalre.network;

import java.util.*;

public class NetworkEdge {

    public NetworkNode parentNode, childNode;
    public BitSet hasSegments;

    public NetworkEdge(NetworkNode parentNode, NetworkNode childNode,
                       BitSet hasSegments) {
        this.parentNode = parentNode;
        this.childNode = childNode;

        this.hasSegments = new BitSet();
        this.hasSegments.or(hasSegments);
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

    public double getReassortmentObsProb() {
        return 1.0 - Math.pow(0.5, hasSegments.cardinality()-1);
    }

    public String getExtendedNewick() {
        return getExtendedNewick(new ArrayList<NetworkNode>()) + ";";
    }

    private String getExtendedNewick(List<NetworkNode> seenReassortmentNodes) {
        StringBuilder result = new StringBuilder();

        boolean traverse = true;
        int hybridID = -1;
        if (childNode.isReassortment()) {
            hybridID = seenReassortmentNodes.indexOf(childNode);

            if (hybridID<0) {
                traverse = false;
                seenReassortmentNodes.add(childNode);
                hybridID = seenReassortmentNodes.size()-1;
            }
        }

        if (traverse && !childNode.isLeaf()) {
            result.append("(");

            boolean isFirst = true;
            for (NetworkEdge childEdge : childNode.getChildEdges()) {
                if (isFirst)
                    isFirst = false;
                else
                    result.append(",");

                result.append(childEdge.getExtendedNewick(seenReassortmentNodes));
            }

            result.append(")");
        }

        if (childNode.getTaxonLabel() != null)
            result.append(childNode.getTaxonLabel());

        if (hybridID>=0) {
            result.append("#H").append(hybridID);
        }

        result.append("[&segments=").append(hasSegments).append("]");

        if (parentNode != null)
            result.append(":").append(parentNode.getHeight() - childNode.getHeight());
        else
            result.append(":0.0");

        return result.toString();
    }
}
