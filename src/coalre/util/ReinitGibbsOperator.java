package coalre.util;

import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;

public class ReinitGibbsOperator extends Operator {

    public Input<StateNode> stateNodeInput = new Input<>("stateNode",
            "StateNode object to Gibbs-sample via reinitialization.",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {

        stateNodeInput.get().startEditing(this);
        stateNodeInput.get().initAndValidate();

        return Double.POSITIVE_INFINITY;
    }

}
