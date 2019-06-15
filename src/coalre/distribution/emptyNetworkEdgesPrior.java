package coalre.distribution;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.distributions.ParametricDistribution;
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
        
        final Double[] arrayForInit = new Double[1];
        arrayForInit[0] = nrEmptyReassortmentEdges;	
        
        Function nrEmptyEdges = new RealParameter(arrayForInit);

        
        // get how many reassortment events are empty
        List<NetworkEdge> emptyEdges = networkIntervalsInput.get().networkInput.get().getEdges().stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality()==0)
				.collect(Collectors.toList());
        
        double overalLength = 0.0;
        for (int i = 0; i < emptyEdges.size(); i++)
        	overalLength += emptyEdges.get(i).getLength();
       
        final Double[] overalLengthForInit = new Double[1];
        overalLengthForInit[0] = nrEmptyReassortmentEdges;	

        Function lengthEmptyEdges = new RealParameter(arrayForInit);
           
        
        logP += nrEventsDistribution.calcLogP(nrEmptyEdges);
        logP += emptyLengthDistribution.calcLogP(lengthEmptyEdges);
        
        return logP;
    }    

}
