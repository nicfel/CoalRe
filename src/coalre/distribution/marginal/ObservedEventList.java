package coalre.distribution.marginal;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class ObservedEventList extends CalculationNode {


    public Input<List<Tree>> segmentTreesInput = new Input<>(
            "segmentTree",
            "Segment tree to compute probability density of.",
            new ArrayList<>());

    private List<ObservedEvent> eventList;
    private List<Pair<Node,Integer>> nodeIdxPairs;

    private boolean listIsDirty;

    public ObservedEventList() {
        eventList = new ArrayList<>();
        nodeIdxPairs = new ArrayList<>();
    };

    @Override
    public void initAndValidate() {
        listIsDirty = true;
    }

    private void update() {
        if (!listIsDirty)
            return;

        nodeIdxPairs.clear();
        for (int segIdx=0; segIdx<segmentTreesInput.get().size(); segIdx++) {
            Node[] nodes = segmentTreesInput.get().get(segIdx).getNodesAsArray();
            for (int i=0; i<nodes.length; i++) {
                nodeIdxPairs.add(new Pair<>(nodes[i], segIdx));
            }
        }

        // Sort entries by their corresponding node ages and leaves by their node numbers.
        // (Allows us to distinguish coincident samples on different lineages, which
        // are not part of the same event - unlike coincident coalescent events.)
        nodeIdxPairs.sort(new Comparator<Pair<Node, Integer>>() {
            @Override
            public int compare(Pair<Node, Integer> p1, Pair<Node, Integer> p2) {

                if (p1.value1.getHeight() < p2.value1.getHeight()) return -1;
                if (p1.value1.getHeight() > p2.value1.getHeight()) return 1;

                if (p1.value1.isLeaf() && p2.value1.isLeaf()) {
                    if (p1.value1.getNr()<p2.value1.getNr()) return -1;
                    if (p1.value1.getNr()>p2.value1.getNr()) return 1;
                }

                return 0;
            }
        });

        ObservedEvent currentEvent = null;
        double currentTime = Double.NEGATIVE_INFINITY;
        int currentTaxonIdx = -1;

        for (Pair<Node,Integer> nodeIdxPair : nodeIdxPairs) {
            Node node = nodeIdxPair.value1;
            int segIdx = nodeIdxPair.value2;

            if (node.getHeight()>currentTime || (node.isLeaf() && node.getNr() != currentTaxonIdx)) {
                currentTime = node.getHeight();

                currentEvent = new ObservedEvent(getNSegments());
                currentEvent.time = currentTime;
                currentEvent.segTreeNodes[segIdx] = node;

                if (node.isLeaf()) {
                    currentEvent.type = ObservedEvent.EventType.SAMPLE;
                    currentEvent.taxonLabel = node.getID();
                    currentEvent.taxonIndex = node.getNr();
                    currentTaxonIdx = node.getNr();
                } else {
                    currentEvent.type = ObservedEvent.EventType.COALESCENCE;
                }

                eventList.add(currentEvent);
            } else {
                currentEvent.segTreeNodes[segIdx] = node;
            }
        }

        listIsDirty = false;
    }

    public List<ObservedEvent> getEventList() {
        update();
        return eventList;
    }

    public int getNSegments() {
        return segmentTreesInput.get().size();
    }

    @Override
    protected void restore() {
        super.restore();

        listIsDirty = true;
    }

    @Override
    protected boolean requiresRecalculation() {
        listIsDirty = true;
        return true;
    }
}
