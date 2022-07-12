package coalre.network;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Tree;
import cern.colt.Arrays;

import java.util.ArrayList;
import java.util.List;

public class SegmentTreeInitializer extends BEASTObject implements StateNodeInitialiser {

    public Input<Network> networkInput = new Input<>("network",
            "Network used to initialize segment trees.", Input.Validate.REQUIRED);

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree to initialize.", new ArrayList<>());

    public Input<Integer> segmentIndexInput = new Input<>("segmentIndex",
            "If only one segment tree is provided, this is its index.");

    int nSegments;
    Network network;
    List<Tree> segmentTrees;
    Integer segmentIndex;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        segmentTrees = segmentTreesInput.get();
        nSegments = network.getSegmentCount();
        segmentIndex = segmentIndexInput.get();

        if (segmentIndex == null) {
            if (segmentTrees.size() != nSegments)
                throw new IllegalArgumentException("Number of segment trees must match number of segments.");
        } else {
            if (segmentTrees.size() != 1)
                throw new IllegalArgumentException("Must only provide one segment tree when segmentIndex is specified.");
            if (segmentIndex >= nSegments)
                throw new IllegalArgumentException("Illegal segment index given.");
        }
    }

    @Override
    public void initStateNodes() {

        if (segmentIndex != null) {
            network.updateSegmentTree(segmentTrees.get(0), segmentIndexInput.get());
        } else {
            for (int segIdx = 0; segIdx < nSegments; segIdx++) {
                network.updateSegmentTree(segmentTrees.get(segIdx), segIdx);
            }
        }

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.addAll(segmentTreesInput.get());
    }
}
