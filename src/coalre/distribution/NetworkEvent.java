package coalre.distribution;

import coalre.network.NetworkEdge;

import java.util.ArrayList;
import java.util.List;

public class NetworkEvent {
    NetworkIntervals.NetworkEventType type;
    double time;
    int lineages;

    List<NetworkEdge> lineagesAdded = new ArrayList<>();
    List<NetworkEdge> lineagesRemoved = new ArrayList<>();
}
