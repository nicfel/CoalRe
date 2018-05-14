package coalre.network;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Objects;

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
