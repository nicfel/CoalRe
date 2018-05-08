package coalre.distribution;

import coalre.network.NetworkNode;

import java.util.Objects;

public class NetworkLineage {
        NetworkNode node;
        boolean isLeft;

        public NetworkLineage(NetworkNode node, boolean isLeft) {
            this.node = node;
            this.isLeft = isLeft;
        }

        public NetworkNode getParentNode() {
            if (isLeft)
                return node.getParent();
            else
                return node.getSecondParent();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            NetworkLineage that = (NetworkLineage) o;
            return isLeft == that.isLeft &&
                    Objects.equals(node, that.node);
        }

        @Override
        public int hashCode() {

            return Objects.hash(node, isLeft);
        }
}
