package coalre.distribution;



import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

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

    public enum NetworkEventType { SAMPLE, COALESCENCE, REASSORTMENT }

    List<NetworkEvent> networkEventList;

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

        // TODO: update event list

        networkEventList = network.getNodes().stream().map(n -> {
            NetworkEvent event = new NetworkEvent();
            event.time = n.getHeight();
            event.node = n;
            switch(n.getChildCount()) {
                case 0:
                    event.type = NetworkEventType.SAMPLE;
                    break;

                case 1:
                    event.type = NetworkEventType.REASSORTMENT;
                    break;

                case 2:
                    event.type = NetworkEventType.COALESCENCE;
                    break;

                default:
                    throw new RuntimeException("Network node has illegal number of children.");
            }
            return event;
        }).sorted(Comparator.comparingDouble(e -> e.time)).collect(Collectors.toList());

        int lineages = 0;

        for (NetworkEvent event : networkEventList) {
            switch(event.type) {
                case SAMPLE:
                case REASSORTMENT:
                    lineages += 1;
                    break;
                case COALESCENCE:
                    lineages -= 1;
                    break;
            }

            event.lineages = lineages;
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