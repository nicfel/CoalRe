package coalre.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

public class NetworkScaleOperator extends NetworkOperator {

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Scale factor tuning parameter.",
            0.8);

    public Input<Boolean> scaleRootOnlyInput = new Input<>(
            "scaleRootOnly",
            "Scale only the root node.",
            false);

    public Input<List<RealParameter>> upParametersInput = new Input<>(
            "upParameter",
            "Parameters to scale in the SAME direction as the network.",
            new ArrayList<>());

    public Input<List<RealParameter>> downParametersInput = new Input<>(
            "downParameter",
            "Parameters to scale in the OPPOSITE direction as the network.",
            new ArrayList<>());
    
    final public Input<Boolean> optimiseInput = new Input<>(
    		"optimise", 
    		"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", 
    		true);


    final public Input<Double> scaleUpperLimit = new Input<>(
    		"upper", 
    		"Upper Limit of scale factor", 
    		1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>(
    		"lower", 
    		"Lower limit of scale factor", 
    		1e-8);


    double scaleFactor;
    boolean scaleRootOnly;
    List<RealParameter> upParameters, downParameters;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        scaleFactor = scaleFactorInput.get();
        scaleRootOnly = scaleRootOnlyInput.get();
        upParameters = upParametersInput.get();
        downParameters = downParametersInput.get();
    }

    @Override
    public double networkProposal() {

        int count = 0;

        double f = scaleFactor + Randomizer.nextDouble()*(1.0/scaleFactor - scaleFactor);

        network.startEditing(this);

        if (scaleRootOnly) {

            // Scale root

            NetworkNode rootNode = network.getRootEdge().childNode;

            rootNode.setHeight(rootNode.getHeight() * f);
            count += 1;

            if (f<1.0) {

                for (NetworkEdge childEdge : rootNode.getChildEdges())
                    if (rootNode.getHeight() < childEdge.childNode.getHeight())
                        return Double.NEGATIVE_INFINITY;

            }

        } else {

            // Scale network nodes

            for (NetworkNode node : network.getInternalNodes()) {
                node.setHeight(node.getHeight() * f);
                count += 1;
            }

            if (f < 1.0) {
                for (NetworkNode leaf : network.getLeafNodes()) {
                    if (leaf.getParentEdges().get(0).parentNode.getHeight() < leaf.getHeight())
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
    	scaleFactor = Math.max(Math.min(value, scaleUpperLimit.get()), scaleLowerLimit.get());
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
