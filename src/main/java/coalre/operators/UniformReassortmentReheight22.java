package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class UniformReassortmentReheight22 extends NetworkOperator {
	
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "magnitude factor used for scaling", Validate.REQUIRED);

        
    final public Input<Boolean> optimiseInput = new Input<>(
    		"optimise", 
    		"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", 
    		true);

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 10.0);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 0.0);
    

    double scaleFactor;
    boolean scaleRootOnly;
    List<RealParameter> upParameters, downParameters;
    double upper, lower;
   


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        scaleFactor = scaleFactorInput.get();
    	
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();


    }
    
    int c1=0,c2=0;
    
    
    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    }

    @Override
    public double networkProposal() {

        int count = 0;
   
        final double f = getScaler();


        network.startEditing(this);
        
        // pick a random coalescent node
        List<NetworkNode> node = network.getInternalNodes().
        		stream().collect(Collectors.toList());
        
        // pick a random node
        NetworkNode randomNode = node.get(Randomizer.nextInt(node.size()));
        
        double logHR = Math.log(1.0/node.size());

        count += scaleNode(randomNode, f);	

        for (NetworkEdge e : network.getEdges()) {
        	if (e.isRootEdge())
        		continue;
        	if (e.getLength() < 0.0) {
        		return Double.NEGATIVE_INFINITY;
            }
        }
        
        
        
        logHR += Math.log(f)*(count);
        return logHR;
    }
    
    private int scaleNode(NetworkNode n, double f) {
    	if (n.isLeaf()) {
    		return 0;
    	}
    	int count = 1;
    	n.setHeight(n.getHeight()*f);
    	for (NetworkEdge edge : n.getChildEdges()) {
			count += scaleNode(edge.childNode, f);
    	}
    	
    	return count;  			
    	
    }
    
    

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


}
