package coalre.statistics;

import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public class NetworkStats {

    public static int getReassortmentCount(Network network) {
        return (int)network.getNodes().stream().filter(NetworkNode::isReassortment).count();
    }

    public static double getTotalEdgeLength(Network network) {
        return network.getEdges().stream().filter(e -> !e.isRootEdge()).
                map(NetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
    }

    public static double getTotalHeight(Network network) {
        return network.getRootEdge().getChildNode().getHeight();
    }
}
