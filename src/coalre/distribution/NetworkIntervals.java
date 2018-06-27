package coalre.distribution;



import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import coalre.network.Network;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Nicola Felix Mueller
 */
public class NetworkIntervals extends CalculationNode {
    public Input<Network> networkInput = new Input<>("network",
            "network for which to calculate the intervals", Validate.REQUIRED);

    private Network network;

    private boolean isDirty;

    private List<NetworkEvent> networkEventList;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        isDirty = true;
    }

    List<NetworkEvent> getNetworkEventList() {
        update();
        return networkEventList;
    }

    void update() {
        if (!isDirty)
            return;

        System.out.println("Updating networkIntervals");

        networkEventList = network.getNodes().stream().map(n -> {
            NetworkEvent event = new NetworkEvent();
            event.time = n.getHeight();
            event.node = n;
            switch(n.getChildCount()) {
                case 0:
                    event.type = NetworkEvent.NetworkEventType.SAMPLE;
                    break;

                case 1:
                    event.type = NetworkEvent.NetworkEventType.REASSORTMENT;
                    break;

                case 2:
                    event.type = NetworkEvent.NetworkEventType.COALESCENCE;
                    break;

                default:
                    throw new RuntimeException("Network node has illegal number of children.");
            }
            return event;
        }).sorted(Comparator.comparingDouble(e -> e.time)).collect(Collectors.toList());

        int lineages = 0;
        double reassortmentObsProb = 0;

        for (NetworkEvent event : networkEventList) {
            switch(event.type) {
                case SAMPLE:
                    lineages += 1;
                    reassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb();
                    break;

                case REASSORTMENT:
                    lineages += 1;
                    reassortmentObsProb -= event.node.getChildEdges().get(0).getReassortmentObsProb();
                    reassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb();
                    reassortmentObsProb += event.node.getParentEdges().get(1).getReassortmentObsProb();
                    break;

                case COALESCENCE:
                    lineages -= 1;
                    reassortmentObsProb -= event.node.getChildEdges().get(0).getReassortmentObsProb();
                    reassortmentObsProb -= event.node.getChildEdges().get(1).getReassortmentObsProb();
                    reassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb();
                    break;
            }

            event.lineages = lineages;
            event.logReassortmentObsProb = reassortmentObsProb;
        }

        isDirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        isDirty = true;

        return true;
    }

    @Override
    protected void restore() {
        isDirty = true;
    }
}