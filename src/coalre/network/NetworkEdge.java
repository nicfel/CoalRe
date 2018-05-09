package coalre.network;

import java.util.Arrays;
import java.util.Objects;

public class NetworkEdge {
    public NetworkNode parentNode, childNode;

    public Boolean[] hasSegments;

    public NetworkEdge(NetworkNode parentNode, NetworkNode childNode,
                       Boolean[] hasSegments) {
        this.parentNode = parentNode;
        this.childNode = childNode;

        this.hasSegments = new Boolean[hasSegments.length];
        //TODO there is probably a better way of doing this
        for (int i = 0; i < hasSegments.length; i++){
            if (hasSegments[i]==null)
                this.hasSegments[i] = false;
            else
                this.hasSegments[i] = hasSegments[i];
        }
    }

    public NetworkNode getParentNode() {
        return parentNode;
    }

    public NetworkNode getChildNode() {
        return childNode;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        NetworkEdge that = (NetworkEdge) o;
        return Objects.equals(parentNode, that.parentNode) &&
                Objects.equals(childNode, that.childNode) &&
                Arrays.equals(hasSegments, that.hasSegments);
    }

    @Override
    public int hashCode() {

        int result = Objects.hash(parentNode, childNode);
        result = 31 * result + Arrays.hashCode(hasSegments);
        return result;
    }
}
