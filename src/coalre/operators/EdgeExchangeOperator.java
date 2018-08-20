package coalre.operators;

import beast.util.Randomizer;
import coalre.network.NetworkEdge;

import java.util.ArrayList;
import java.util.List;

public class EdgeExchangeOperator extends DivertSegmentOperator {

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public double networkProposal() {
        double logHR = 0.0;

        // Select two edges at random

        List<NetworkEdge> edges = new ArrayList<>(network.getEdges());

        NetworkEdge edge1, edge2;
        do {
            edge1 = edges.get(Randomizer.nextInt(edges.size()));
        } while (edge1.isRootEdge());

        do {
            edge2 = edges.get(Randomizer.nextInt(edges.size()));
        } while(edge2.isRootEdge() || edge2 == edge1);

        // Check for possibility of exchange

        // Update segment mapping (only redraw segments that are not shared)

        // Perform exchange

        return logHR;
    }
}
