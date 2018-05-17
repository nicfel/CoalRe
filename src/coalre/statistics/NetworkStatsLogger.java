package coalre.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import coalre.network.Network;

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

        out.print(NetworkStats.getTotalHeight(network) + "\t" +
                NetworkStats.getTotalEdgeLength(network) + "\t" +
                NetworkStats.getReassortmentCount(network) + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
