package coalre.operators;

import beast.core.Input;
import beast.core.Operator;
import coalre.network.Network;

public abstract class NetworkOperator extends Operator {

    public Input<Network> networkInput = new Input<>("network",
            "Network on which to operate",
            Input.Validate.REQUIRED);
}
