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
        
        networkInput.get().baseName = segmentTreesInput.get().get(0).getID();
        for (int segIdx=1; segIdx<nSegments; segIdx++) {
    		String newBase = "";
        	for (int i = 0; i < segmentTreesInput.get().get(segIdx).getID().length(); i++) {
        		if (i>=networkInput.get().baseName.length()) {
        			networkInput.get().baseName = newBase;
        		}            			
				if (networkInput.get().baseName.substring(0,i+1).contentEquals(segmentTreesInput.get().get(segIdx).getID().substring(0, i+1))) {
					newBase = networkInput.get().baseName.substring(0,i+1);
				}
				
        	}
        	networkInput.get().baseName = newBase;
        }         

        networkInput.get().segmentNames = new String[nSegments];        
        for (int segIdx=0; segIdx<nSegments; segIdx++) {
        	networkInput.get().segmentNames[segIdx] = segmentTreesInput.get().get(segIdx).getID().replace(networkInput.get().baseName, "");
        }
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
