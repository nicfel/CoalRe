package coalre.network;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

public class SegmentTreeInitializer extends BEASTObject implements StateNodeInitialiser {

    public Input<Network> networkInput = new Input<>("network",
            "Network used to initialize segment trees.", Input.Validate.REQUIRED);

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree to initialize.", new ArrayList<>());

    int nSegments;
    Network network;
    List<Tree> segmentTrees;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        segmentTrees = segmentTreesInput.get();
        nSegments = network.getSegmentCount();

        if (segmentTrees.size() != nSegments)
            throw new IllegalArgumentException("Number of segment trees must match number of segments.");
    }

    @Override
    public void initStateNodes() {

        for (int segIdx=0; segIdx<nSegments; segIdx++) {
            network.updateSegmentTree(segmentTrees.get(segIdx), segIdx);
        }

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.addAll(segmentTreesInput.get());
    }
}
