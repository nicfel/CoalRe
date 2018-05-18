package coalre.operators;

import beast.core.Input;
import beast.core.Operator;
import coalre.network.Network;

public abstract class NetworkOperator extends Operator {

    public static final double LOG2 = Math.log(2.0);

    public Input<Network> networkInput = new Input<>("network",
            "Network on which to operate",
            Input.Validate.REQUIRED);
}
