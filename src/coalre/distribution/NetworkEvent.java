package coalre.distribution;

import java.util.ArrayList;
import java.util.List;

public class NetworkEvent {
    NetworkIntervals.NetworkEventType type;
    double time;
    int lineages;

    List<NetworkLineage> lineagesAdded = new ArrayList<>();
    List<NetworkLineage> lineagesRemoved = new ArrayList<>();
}
