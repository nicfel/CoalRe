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

    public NetworkStatsLogger() { }

    @Override
    public void initAndValidate() {
    }

    @Override
    public void init(PrintStream out) {

        String prefix = getID() == null ? "networkStat." : getID() + ".";

        out.print(prefix + "height\t" +
                prefix + "totalLength\t" +
                prefix + "reassortmentNodeCount\t");
    }

    @Override
    public void log(long sample, PrintStream out) {

        out.print(NetworkStats.getTotalHeight(networkInput.get()) + "\t" +
                NetworkStats.getTotalEdgeLength(networkInput.get()) + "\t" +
                NetworkStats.getReassortmentCount(networkInput.get()) + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
