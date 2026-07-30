package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

public class NetworkScaleSkyline extends NetworkOperator {

    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "magnitude factor used for scaling", Validate.REQUIRED);

    public Input<List<RealParameter>> upParametersInput = new Input<>(
            "upParameter",
            "Parameters to scale in the SAME direction as the network.",
            new ArrayList<>());

    public Input<List<RealParameter>> downParametersInput = new Input<>(
            "downParameter",
            "Parameters to scale in the OPPOSITE direction as the network.",
            new ArrayList<>());
    
    public Input<RealParameter> NeInput = new Input<>(
            "Ne",
            "Ne input to scale in the same direction as the network, but with different scale factor.");
    
    final public Input<Boolean> optimiseInput = new Input<>(
    		"optimise", 
    		"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", 
    		true);

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 10.0);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 0.0);
    
    public final Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());

    protected KernelDistribution kernelDistribution;


    double scaleFactor;
    List<RealParameter> upParameters, downParameters;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        scaleFactor = scaleFactorInput.get();
        upParameters = upParametersInput.get();
        downParameters = downParametersInput.get();
    	kernelDistribution = kernelDistributionInput.get();

    }
    
    int c1=0,c2=0;
    
    
	protected double getScaler(int i) {
		return kernelDistribution.getScaler(i, Double.NaN, getCoercableParameterValue());
	}

    @Override
    public double networkProposal() {

        int count = 0;
   
        final double f = getScaler(0);


        network.startEditing(this);

        for (NetworkNode node : network.getInternalNodes()) {
            node.setHeight(node.getHeight() * f);
            count += 1;
        }

        if (f < 1.0) {
            for (NetworkNode leaf : network.getLeafNodes()) {
                if (leaf.getParentEdges().get(0).parentNode.getHeight() < leaf.getHeight()) {
                    c2++;

                    return Double.NEGATIVE_INFINITY;
                }
            }
        }


        // Scale parameters
        try {

            for (RealParameter param : upParameters) {
                param.startEditing(this);
                count += param.scale(f);
            }

            for (RealParameter param : downParameters) {
                param.startEditing(this);
                count -= param.scale(1.0 / f);
            }

        } catch (IllegalArgumentException ex) {
            return Double.NEGATIVE_INFINITY;
        }
        return Math.log(f)*(count-2);
    }
    
    /** 
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value,scaleUpperLimit.get()),scaleLowerLimit.get());
    }

    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }

    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }
    

}
