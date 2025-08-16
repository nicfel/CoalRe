package coalre.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import coalre.network.Network;


@Description("calculates the differences between the entries of a vector")
public class emptyNetworkEdgesPrior extends CalculationNode implements Function {
    final public Input<Network> networkInput = new Input<>("network", "ignore difference after that index");


    boolean needsRecompute = true;
    double emptyEdgesCount;
    double storedEmptyEdgesCount;
    

    @Override
    public void initAndValidate() {

    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return emptyEdgesCount;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	emptyEdgesCount = (double) networkInput.get().getEdges().stream()
    			.filter(e -> e.hasSegments.cardinality() == 0)
    			.count();
    	
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return emptyEdgesCount;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
		storedEmptyEdgesCount = emptyEdgesCount;
        super.store();
    }

    @Override
    public void restore() {
    	double tmp = storedEmptyEdgesCount;
    	storedEmptyEdgesCount = emptyEdgesCount;
		emptyEdgesCount = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }
} // class Sum
