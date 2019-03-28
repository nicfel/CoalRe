package coalre.distribution;



import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import coalre.network.Network;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan and Nicola Felix Mueller
 */
public class NetworkIntervals extends CalculationNode {
    public Input<Network> networkInput = new Input<>("network",
            "network for which to calculate the intervals", Validate.REQUIRED);

    public Input<RealParameter> binomialProbInput = new Input<>("binomialProb",
            "Probability of a given segment choosing a particular parent.");

    private Network network;

    private List<NetworkEvent> networkEventList, storedNetworkEventList;

    public boolean eventListDirty = true;

    @Override
    public void initAndValidate() {
        network = networkInput.get();

        storedNetworkEventList = new ArrayList<>();
    }

    List<NetworkEvent> getNetworkEventList() {
        update();

        return networkEventList;
    }

    public double getBinomialProb() {
        return binomialProbInput.get() != null
                ? binomialProbInput.get().getArrayValue()
                : 0.5;
    }

    void update() {
        if (!eventListDirty)
            return;

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
        double totalReassortmentObsProb = 0;

        for (NetworkEvent event : networkEventList) {
            switch(event.type) {
                case SAMPLE:
                    lineages += 1;
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    break;

                case REASSORTMENT:
                    lineages += 1;
                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(1).getReassortmentObsProb(getBinomialProb());

                    event.segsToSort = event.node.getChildEdges().get(0).hasSegments.cardinality();
                    event.segsSortedLeft = event.node.getParentEdges().get(0).hasSegments.cardinality();
                    break;

                case COALESCENCE:
                    lineages -= 1;
                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb -= event.node.getChildEdges().get(1).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    break;
            }

            event.lineages = lineages;
            event.totalReassortmentObsProb = totalReassortmentObsProb;
        }

        eventListDirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        eventListDirty = true;

        return true;
    }

    @Override
    protected void restore() {
        List<NetworkEvent> tmp = networkEventList;
        networkEventList = storedNetworkEventList;
        storedNetworkEventList = tmp;

        super.restore();
    }

    @Override
    protected void store() {
        storedNetworkEventList.clear();
        storedNetworkEventList.addAll(networkEventList);

        super.store();
    }
}