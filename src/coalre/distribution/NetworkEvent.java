package coalre.distribution;

import coalre.network.NetworkNode;

import java.util.BitSet;

public class NetworkEvent {
    public enum NetworkEventType { SAMPLE, COALESCENCE, REASSORTMENT }

    public NetworkEventType type;
    public double time;

    /**
     * Number of segments on a reassorting lineage.
     */
    public int segsToSort;
    
    /**
     * Which segments went left, needed for plasmids.
     */
    public BitSet segsLeft;
    
    /**
     * whcih segments went right, needed for plasmids.
     */
    public BitSet segsRight;
    
    public boolean empty;



    /**
     * Number of segments sent to the first parent.
     */
    public int segsSortedLeft;

    public int lineages;
    public double totalReassortmentObsProb;

    /**
     * Only used when setting up event list.
     * May not point to a compatible node at other times.
     */
    public NetworkNode node;
    
    @Override
    public String toString() {
		return type +":"+ lineages +":" +totalReassortmentObsProb;
    }
}
