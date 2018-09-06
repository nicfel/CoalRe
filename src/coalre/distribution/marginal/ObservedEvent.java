package coalre.distribution.marginal;

import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class ObservedEvent {
    public enum EventType {COALESCENCE, SAMPLE};
    public EventType type;
    public double time;
    public Node[] segTreeNodes;

    public String taxonLabel;
    public int taxonIndex;

    public ObservedEvent(int nSegments) {
        segTreeNodes = new Node[nSegments];
    }
}
