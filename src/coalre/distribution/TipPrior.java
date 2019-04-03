package coalre.distribution;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;
import coalre.network.Network;
import coalre.network.NetworkNode;


@Description("Prior over set of taxa, useful for defining monophyletic constraints and "
        + "distributions over MRCA times or (sets of) tips of trees")
public class TipPrior extends Distribution {
    public final Input<Network> networkInput = new Input<>("network", "the network containing the taxon set", Validate.REQUIRED);
    public final Input<TaxonSet> taxonsetInput = new Input<>("taxonset",
            "set of taxa for which prior information is available");
    public final Input<Boolean> isMonophyleticInput = new Input<>("monophyletic",
            "whether the taxon set is monophyletic (forms a clade without other taxa) or nor. Default is false.", false);
    public final Input<ParametricDistribution> distInput = new Input<>("distr",
            "distribution used to calculate prior over MRCA time, "
                    + "e.g. normal, beta, gamma. If not specified, monophyletic must be true");
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree associated with network.",
            new ArrayList<>());
    public Input<RealParameter> dateOffsetInput = new Input<>("dateOffset", 
    		"keeps track of how much the dates have change", Validate.REQUIRED);


//    public final Input<Boolean> onlyUseTipsInput = new Input<>("tipsonly",
//            "flag to indicate tip dates are to be used instead of the MRCA node. " +
//                    "If set to true, the prior is applied to the height of all tips in the taxonset " +
//                    "and the monophyletic flag is ignored. Default is false.", false);
//    public final Input<Boolean> useOriginateInput = new Input<>("useOriginate", "Use parent of clade instead of clade. Cannot be used with tipsonly, or on the root.", false);

    /**
     * shadow members *
     */
    ParametricDistribution dist;
    Network network;
    
    // number of taxa in taxon set
    int nrOfTaxa = -1;
    // array of flags to indicate which taxa are in the set
    Set<String> isInTaxaSet = new LinkedHashSet<>();

    // array of indices of taxa
    int[] taxonIndex;
    // stores time to be calculated
    double MRCATime = -1;
    double storedMRCATime = -1;
    
    boolean initialised = false;
    
    NetworkNode operatingNode;

    @Override
    public void initAndValidate() {
        dist = distInput.get();
        network = networkInput.get();
        
        if (taxonsetInput.get().asStringList().size()!= 1) {
        	throw new IllegalArgumentException("TipPrior expects the number of tips to be 1");
        }

        for (final NetworkNode taxon : network.getLeafNodes()) {
        	if (taxon.getTaxonLabel()==taxonsetInput.get().getTaxonId(0)){
        		operatingNode = taxon;
        		break;
        	}
        }    
        
        initialised = false;
    }

    boolean [] nodesTraversed;
    int nseen;


    // A lightweight version for finding the most recent common ancestor of a group of taxa.
    // return the node-ref of the MRCA.

    @Override
    public double calculateLogP() {
    	if (!initialised) {
    		initialise();
    	}
        logP = 0;
        
        // tip date
    	if (dist == null) {
    		return logP;
    	}
    	
        for (final NetworkNode taxon : network.getLeafNodes()) {
        	if (taxon.getTaxonLabel()==taxonsetInput.get().getTaxonId(0)){
        		operatingNode = taxon;
                MRCATime = segmentTreesInput.get().get(0).getDate(operatingNode.getHeight() + dateOffsetInput.get().getArrayValue());
                logP += dist.logDensity(MRCATime);
        		break;
        	}
        }    
        
//        System.out.println(logP);

        return logP;
    }
    
    public void initialise() {
        dist = distInput.get();
        network = networkInput.get();
        
        
        final List<String> taxaNames = new ArrayList<>();
        for (final NetworkNode taxon : network.getLeafNodes()) {
            taxaNames.add(taxon.getTaxonLabel());
        }
        
        // determine nr of taxa in taxon set
        List<String> set = null;
        if (taxonsetInput.get() != null) {
            set = taxonsetInput.get().asStringList();
            nrOfTaxa = set.size();
        } else {
            // assume all taxa
            nrOfTaxa = taxaNames.size();
        }

        if (nrOfTaxa == 1) {
        }else{
        	throw new IllegalArgumentException("TipPrior expects the number of tips to be 1");
        }
        
        initialised = false;
    }



    @Override
    public void store() {
        storedMRCATime = MRCATime;
        // don't need to store m_bIsMonophyletic since it is never reported
        // explicitly, only logP and MRCA time are (re)stored
        super.store();
    }

    @Override
    public void restore() {
        MRCATime = storedMRCATime;
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation();
    }


    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
    	if (!initialised) {
    		initialise();
    	}
        if (dist != null) {
            out.print("logP(mrca(" + getID() + "))\t");
        }
        out.print("height(" + operatingNode.getTaxonLabel() + ")\t");
        
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        if (dist != null) {
            out.print(getCurrentLogP() + "\t");
        }
        out.print(MRCATime + "\t");
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }

    /**
     * Valuable interface implementation follows, first dimension is log likelihood, second the time *
     */
    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public double getArrayValue() {
    	if (Double.isNaN(logP)) {
    		try {
    			calculateLogP();
    		}catch (Exception e) {
    			logP  = Double.NaN;
    		}
    	}
        return logP;
    }

    @Override
    public double getArrayValue(final int dim) {
    	if (Double.isNaN(logP)) {
    		try {
    			calculateLogP();
    		}catch (Exception e) {
    			logP  = Double.NaN;
    		}
    	}
        switch (dim) {
            case 0:
                return logP;
            case 1:
                return MRCATime;
            default:
                return 0;
        }
    }

    @Override
    public void sample(final State state, final Random random) {
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }
}