package coalre.distribution;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;

import java.util.List;
import java.util.Random;

public class NetworkDistribution extends Distribution {
    public Input<NetworkIntervals> networkIntervalsInput = new Input<>("networkIntervals",
            "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	if (networkIntervalsInput.get().isDirty())
    		return true;
    	return false;
    }


}
