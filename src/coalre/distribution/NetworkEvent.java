package coalre.distribution;

import coalre.network.NetworkNode;

public class NetworkEvent {
    NetworkIntervals.NetworkEventType type;
    double time;
    NetworkNode node;

    int lineages;
    double logReassortmentObsProb;
}
