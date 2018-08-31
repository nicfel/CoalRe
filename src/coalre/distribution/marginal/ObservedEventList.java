package coalre.distribution.marginal;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class ObservedEventList extends CalculationNode {


    public Input<List<Tree>> segmentTreesInput = new Input<>(
            "segmentTree",
            "Segment tree to compute probability density of.",
            new ArrayList<>());

    private List<ObservedEvent> eventList;
    private List<Node> nodeList;

    private boolean listIsDirty;

    public ObservedEventList() {
        eventList = new ArrayList<>();
        nodeList = new ArrayList<>();
    };

    @Override
    public void initAndValidate() {
        listIsDirty = true;
    }

    private void update() {
        if (!listIsDirty)
            return;

        nodeList.clear();
        for (Tree segmentTree : segmentTreesInput.get()) {
            nodeList.addAll(Arrays.asList(segmentTree.getNodesAsArray()));
        }
        nodeList.sort(Comparator.comparingDouble(Node::getHeight));

        ObservedEvent currentEvent = null;
        double currentTime = Double.NEGATIVE_INFINITY;

        for (Node node : nodeList) {
            if (node.getHeight()>currentTime) {
                currentTime = node.getHeight();

                currentEvent = new ObservedEvent();
                currentEvent.time = currentTime;
                currentEvent.nodes.add(node);

                eventList.add(currentEvent);
            } else {
                currentEvent.nodes.add(node);
            }
        }

        listIsDirty = false;
    }

    public List<ObservedEvent> getEventList() {
        update();
        return eventList;
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
