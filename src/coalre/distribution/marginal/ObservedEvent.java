package coalre.distribution.marginal;

import beast.evolution.tree.Node;

import java.util.HashSet;
import java.util.Set;

public class ObservedEvent {
    public enum EventType {COALESCENCE, SAMPLE};
    public EventType type;
    public double time;
    public Set<Node> nodes = new HashSet<>();

}
