package coalre.util;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

public class DateOffsetInitializer extends BEASTObject implements StateNodeInitialiser {

    public Input<RealParameter> dateOffsetInput = new Input<>("dateOffset", 
    		"keeps track of how much the dates have change", Validate.REQUIRED);

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree to initialize.", new ArrayList<>());

    List<Tree> segmentTrees;

    @Override
    public void initAndValidate() {
    	
        segmentTrees = segmentTreesInput.get();
    }

    @Override
    public void initStateNodes() {
    	double time = 0.0;

        for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++) {
            time = Math.max(time, segmentTreesInput.get().get(0).getDate(0.0));
        }
        
        dateOffsetInput.get().setValue(time);

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(dateOffsetInput.get());
    }
}
