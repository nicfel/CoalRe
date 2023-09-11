package coalre.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.List;


/**
 * Implements the subnet slide move. General workflow:
 * 1. Choose an edge to move and a child it will carry
 * 2. Make a copy of subnet starting with this child edge
 * 3. Attach a new coppy to the new parent with randomly drawn height
 * 4. Rearrange segments
 * 5. Delete subnet starting at the child in the old position
 */
@Description("Moves the height of an internal node along the branch. " +
        "If it moves up, it can exceed the root and become a new root. " +
        "If it moves down, it may need to make a choice which branch to " +
        "slide down into.")
public class TipReheight extends NetworkOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 0.1", 0.1);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    
    public final Input<TaxonSet> taxonsetInput = new Input<>("taxonset",
            "set of taxa for which prior information is available");
    public final Input<Boolean> isMonophyleticInput = new Input<>("monophyletic",
            "whether the taxon set is monophyletic (forms a clade without other taxa) or nor. Default is false.", false);
    public final Input<ParametricDistribution> distInput = new Input<>("distr",
            "distribution used to calculate prior over MRCA time, "
                    + "e.g. normal, beta, gamma. If not specified, monophyletic must be true");
    public Input<RealParameter> dateOffset = new Input<>("dateOffset", "keeps track of how much the dates have change", Validate.REQUIRED);

    // shadows size
    double size;
    NetworkNode operatingNode;
    List<Node> segmentTreeNodes;


	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        
        if (taxonsetInput.get().asStringList().size()!= 1) {
        	throw new IllegalArgumentException("TipPrior expects the number of tips to be 1");
        }
        
        for (final NetworkNode taxon : network.getLeafNodes()) {
        	if (taxon.getTaxonLabel().equals(taxonsetInput.get().getTaxonId(0))){
        		operatingNode = taxon;
        		break;
        	}
        }    
	}

	@Override
	public double networkProposal() {

		double logHR = 0.0;
		network.startEditing(this);
		
        for (final NetworkNode taxon : network.getLeafNodes()) {
        	if (taxon.getTaxonLabel().equals(taxonsetInput.get().getTaxonId(0))){
        		operatingNode = taxon;

                // 2. choose a delta to move
                final double delta = getDelta();
                final double oldHeight = operatingNode.getHeight();
                final double newHeight = oldHeight + delta;
                                
                // check if ne height is valid
                if (operatingNode.getParentEdges().get(0).parentNode.getHeight()<=newHeight)
                	return Double.NEGATIVE_INFINITY;
                
                // set the height of the leaf node, even if the height drops below 0
                operatingNode.setHeight(newHeight);               
                
                // if the height drops below 0, reheight the whole network (probably inefficient, due to making nodes dirty)
                if (newHeight < 0){
                	double diff = newHeight;
                	for (NetworkNode node : network.getNodes())
                		node.setHeight(node.getHeight()-diff);  
                	
                	dateOffset.get().startEditing(this);
                	dateOffset.get().setValue(dateOffset.get().getValue()+diff);
                	
                }else if (oldHeight==0){
                	// get the second lowest height    
                	double minHeight = Double.POSITIVE_INFINITY;
                	for (NetworkNode node : network.getLeafNodes()){
                		if (!node.equals(operatingNode)){
                			if (node.getHeight()<minHeight)
                				minHeight = node.getHeight();
                		}  
                	}
                	
                	// rescale all internal nodes relative to the newest most recently sampled individual
                	if (newHeight>minHeight){
                    	for (NetworkNode node : network.getNodes())
                    		node.setHeight(node.getHeight()-minHeight);                	
                	}
                	
                	dateOffset.get().startEditing(this);
                	dateOffset.get().setValue(dateOffset.get().getValue()+minHeight);
                }

        		break;
        	}
        }    

		return logHR;
	}
	
		

    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }
    
    
    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
            size = f;
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }

}
