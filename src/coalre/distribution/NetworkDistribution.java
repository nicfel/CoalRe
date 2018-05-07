
package coalre.distribution;


import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;

// should be the same as the original tree distribution but allowing for Structured tree intervals 
// additionally allows recalculation of which daughter lineages were involved in a coalescent event
@Description("Distribution on a tree, typically a prior such as Coalescent or Yule")
public class NetworkDistribution extends Distribution {
    public Input<NetworkIntervals> networkIntervalsInput = new Input<NetworkIntervals>("networkIntervals",
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
        final NetworkIntervals ti = networkIntervalsInput.get();
        assert ti.isDirtyCalculation();
        return true;
    }
    
 	/** Indicate that the tree distribution can deal with dated tips in the tree
	 * Some tree distributions like the Yule prior cannot handle this.
	 * @return true by default
	 */
	public boolean canHandleTipDates() {
		return true;
	}
}
