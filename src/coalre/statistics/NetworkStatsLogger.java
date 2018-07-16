package coalre.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.PrintStream;

public class NetworkStatsLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);

    Network network;

    public NetworkStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
    }

    @Override
    public void init(PrintStream out) {

        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

        out.print(prefix + "height\t" +
                prefix + "totalLength\t" +
                prefix + "reassortmentNodeCount\t");
    }

    @Override
    public void log(long sample, PrintStream out) {

        out.print(getTotalHeight(network) + "\t" +
                getTotalEdgeLength(network) + "\t" +
                getReassortmentCount(network) + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }

    public static int getReassortmentCount(Network network) {
        return (int)network.getNodes().stream().filter(NetworkNode::isReassortment).count();
    }

    public static double getTotalEdgeLength(Network network) {
        return network.getEdges().stream().filter(e -> !e.isRootEdge()).
                map(NetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
    }

    public static double getTotalHeight(Network network) {
        return network.getRootEdge().childNode.getHeight();
    }
}
