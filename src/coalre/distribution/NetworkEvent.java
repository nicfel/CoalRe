package coalre.distribution;

import coalre.network.NetworkNode;

public class NetworkEvent {
    public enum NetworkEventType { SAMPLE, COALESCENCE, REASSORTMENT }

    public NetworkEventType type;
    public double time;
    public NetworkNode node;

    public int lineages;
    public double totalReassortmentObsProb;
}
