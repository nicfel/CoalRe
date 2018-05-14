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

    public String getExtendedNewick() {
        return getExtendedNewick(new ArrayList<NetworkNode>());
    }

    private String getExtendedNewick(List<NetworkNode> seenReassortmentNodes) {
        String result = "";

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

        if (traverse) {
            if (!childNode.children.isEmpty()) {
                result += "(";

                for (NetworkEdge childEdge : childNode.children) {
                    result += childEdge.getExtendedNewick(seenReassortmentNodes);
                }

                result += ")";
            }
        }

        if (hybridID>=0) {
            result += "#H" + hybridID;
        }

        if (parentNode != null)
            result += ":" + (parentNode.getHeight() - childNode.getHeight());
        else
            result += ":0.0";

        result += "[&segments=" + hasSegments.toString() + "]";

        return result;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        NetworkEdge that = (NetworkEdge) o;
        return Objects.equals(parentNode, that.parentNode) &&
                Objects.equals(childNode, that.childNode) &&
                Objects.equals(hasSegments, that.hasSegments);
    }

    @Override
    public int hashCode() {

        return Objects.hash(parentNode, childNode, hasSegments);
    }
}
