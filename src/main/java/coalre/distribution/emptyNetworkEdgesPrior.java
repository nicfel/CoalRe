package coalre.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import beast.base.inference.distribution.ParametricDistribution;
import coalre.network.NetworkEdge;

import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a reassortment network using under" +
        " the framework of Mueller (2018).")
public class emptyNetworkEdgesPrior extends NetworkDistribution {

	final public Input<ParametricDistribution> nrEventsDistributionInput = new Input<>("nrEventsDistribution", "distribution used to calculate prior, e.g. normal, beta, gamma.", Validate.REQUIRED);
	final public Input<ParametricDistribution> emptyLengthDistributionInput = new Input<>("emptyLengthDistribution", "distribution used to calculate prior, e.g. normal, beta, gamma.", Validate.REQUIRED);

	/**
     * shadows distInput *
     */
    protected ParametricDistribution nrEventsDistribution;
    protected ParametricDistribution emptyLengthDistribution;

    @Override
    public void initAndValidate(){
    	nrEventsDistribution = nrEventsDistributionInput.get();
    	emptyLengthDistribution = emptyLengthDistributionInput.get();
        calculateLogP();
    }

    public double calculateLogP() {
    	
    	logP = 0.0;
    	
    	// get how many reassortment edges are empty
        double nrEmptyReassortmentEdges = (double) networkIntervalsInput.get().networkInput.get().getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()==0)
                .filter(e -> e.childNode.isReassortment())
                .count();
        
        RealScalar<NonNegativeReal> nrEmptyEdges =
                new RealScalarParam<>(nrEmptyReassortmentEdges, NonNegativeReal.INSTANCE);

        
        // get how many reassortment events are empty
        List<NetworkEdge> emptyEdges = networkIntervalsInput.get().networkInput.get().getEdges().stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality()==0)
				.collect(Collectors.toList());
        
        double overalLength = 0.0;
        for (int i = 0; i < emptyEdges.size(); i++)
        	overalLength += emptyEdges.get(i).getLength();

        // NOTE: this deliberately reproduces the pre-migration behaviour, which looks like a
        // bug. The BEAST 2 code assigned nrEmptyReassortmentEdges (not overalLength) into
        // overalLengthForInit, and then built lengthEmptyEdges from arrayForInit anyway -- so
        // overalLength was computed and discarded, and emptyLengthDistribution scored the
        // empty-edge COUNT, not the total length. Preserved verbatim because changing it would
        // alter this model's likelihood; flagged in tmp/b3migration/NOTES.md for the author.
        RealScalar<NonNegativeReal> lengthEmptyEdges =
                new RealScalarParam<>(nrEmptyReassortmentEdges, NonNegativeReal.INSTANCE);
           
        
        logP += nrEventsDistribution.calcLogP(nrEmptyEdges);
        logP += emptyLengthDistribution.calcLogP(lengthEmptyEdges);
        
        return logP;
    }    

}
