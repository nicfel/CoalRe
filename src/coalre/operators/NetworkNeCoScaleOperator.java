package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class NetworkNeCoScaleOperator extends NetworkOperator {

	final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "magnitude factor used for scaling",
			Validate.REQUIRED);

	public Input<Boolean> scaleRootOnlyInput = new Input<>("scaleRootOnly", "Scale only the root node.", false);

	public Input<List<RealParameter>> upParametersInput = new Input<>("upParameter",
			"Parameters to scale in the SAME direction as the network.", new ArrayList<>());

	public Input<List<RealParameter>> downParametersInput = new Input<>("downParameter",
			"Parameters to scale in the OPPOSITE direction as the network.", new ArrayList<>());

	final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);

	final public Input<Boolean> scaleIntervalsIndependently = new Input<>("scaleIntervalsIndependently",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			false);

//    final public Input<Boolean> useprime = new Input<>(
//    		"useprime",
//    		"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
//    		true);

	final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 10.0);
	final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 0.0);

	double scaleFactor;
	boolean scaleRootOnly;
	List<RealParameter> upParameters, downParameters;
	double upper, lower;

	Map<NetworkNode, Double> oldLengths;
	Map<NetworkNode, Boolean> updated;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		scaleFactor = scaleFactorInput.get();
		scaleRootOnly = scaleRootOnlyInput.get();
		upParameters = upParametersInput.get();
		downParameters = downParametersInput.get();

		upper = scaleUpperLimit.get();
		lower = scaleLowerLimit.get();

	}

	int c1 = 0, c2 = 0;

	protected double getScaler() {
		return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
	}

	@Override
    public double networkProposal() {

        int count = 0;
   
        final double f = getScaler();


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
        	c1++;

//        	double fprime;
//        	if (useprime.get()) {
//        		double networkLength = 0.0;
//            	double sumLeafHeights = 0.0;
//        		if (segmentTreesInput.get()!=null) {
//        			for (Tree t : segmentTreesInput.get()) {
//        				networkLength += TreeUtils.getTreeLength(t, t.getRoot());
//                        // Scale network nodes
//                    	for (Node n : t.getInternalNodes()) {
//                    		if (n.isLeaf())
//                    			sumLeafHeights += n.getHeight();
//                    	}
//        			}
//        		}else {
//                    // Scale network nodes
//                	for (NetworkNode n : network.getLeafNodes()) {
//                		sumLeafHeights += n.getHeight();
//                	}
//        			networkLength = coalre.statistics.NetworkStatsLogger.getTotalEdgeLength(network);
//        		}
//
//        		fprime = (f*networkLength + sumLeafHeights)/(networkLength + sumLeafHeights);
//        	}else {
//        		fprime = f;
//        	}
        	
        	
        	// calculate the network length
//        	List<Double> l_before = new ArrayList<>();
//        	// loop over all segments
//			for (int s = 0; s < network.getSegmentCount(); s++) {
//				l_before.add(0.0);
//			}
//			for (NetworkEdge edge : network.getEdges()) {
//				if (edge.isRootEdge() )
//					continue;
//				for (int s = 0; s < network.getSegmentCount(); s++) {
//					if (edge.hasSegments.get(s))
//						l_before.set(s, l_before.get(s) + edge.getLength());
//				}				
//				
//			}
			
        	
        	
        	
        	oldLengths = new HashMap<>(); 
        	updated = new HashMap<>();

            for (NetworkNode node : network.getInternalNodes()) {
            	if (node.isLeaf())
            		continue;
            
            	if (node.isReassortment()) 
            		oldLengths.put(node, node.getChildEdges().get(0).getLength());
            	else
            		oldLengths.put(node, node.getHeight() - Math.max(node.getChildEdges().get(0).childNode.getHeight(), node.getChildEdges().get(1).childNode.getHeight()));            		
            	
            	updated.put(node, false);
            }
//            System.out.println(network);
            if (scaleIntervalsIndependently.get()) {
            	double lhr = resampleNodeHeight(network.getRootEdge().childNode);
            	return lhr;
            }else {
            	count += scaleNodes(network.getRootEdge().childNode, f);
            }
            
//            List<Double> l_after = new ArrayList<>();
//            // loop over all segments
//            for (int s = 0; s < network.getSegmentCount(); s++) {
//            	l_after.add(0.0);
//            }
//            for (NetworkEdge edge : network.getEdges()) {
//                if (edge.isRootEdge() )
//                    continue;
//                for (int s = 0; s < network.getSegmentCount(); s++) {
//                    if (edge.hasSegments.get(s))
//                        l_after.set(s, l_after.get(s) + edge.getLength());
//                }
//                
//            }
//
//        	System.out.println(network);
//
//			System.out.print(" scale factor: " + f);
//			for (int s = 0; s < network.getSegmentCount(); s++)
//                System.out.print(" segment " + s + ": " + l_after.get(s)/l_before.get(s));
            
			
			
            
//            System.out.println(network);
//            System.out.println(f);
//            System.out.println(count);
//            System.exit(0);
            

//            if (f < 1.0) {
//                for (NetworkNode leaf : network.getLeafNodes()) {
//                    if (leaf.getParentEdges().get(0).parentNode.getHeight() < leaf.getHeight()) {
//                    	c2++;
//                    	
//                        return Double.NEGATIVE_INFINITY;
//                    }
//                }
//            }            
        }

	// Scale parameters
	try

	{

		for (RealParameter param : upParameters) {
			param.startEditing(this);
			count += param.scale(f);
		}

		for (RealParameter param : downParameters) {
			param.startEditing(this);
			count -= param.scale(1.0 / f);
		}

	}catch(
	IllegalArgumentException ex)
	{
		return Double.NEGATIVE_INFINITY;
	}return Math.log(f)*(count-2);
	}

	private int scaleNodes(NetworkNode n, double scaler) {
    	int count = 0;
		if (n.isLeaf())
			return 0;
		else {
			if (n.isReassortment()) {
				// check if the edge below the reassortment edge hast already been scaled
				if (!updated.get(n)) {
					count += scaleNodes(n.getChildEdges().get(0).childNode, scaler);
					double newLength = oldLengths.get(n)*scaler;
					n.setHeight(newLength + n.getChildEdges().get(0).childNode.getHeight());
					updated.put(n, true);				
					count++;

				}
				// update the 
			}
			else {
				for (NetworkEdge e : n.getChildEdges()) {
					count += scaleNodes(e.childNode, scaler);
				}
				double newLength = oldLengths.get(n)*scaler;
				n.setHeight(newLength + Math.max(n.getChildEdges().get(0).childNode.getHeight(), n.getChildEdges().get(1).childNode.getHeight()));
				updated.put(n, true);	
				count++;
			}
			// scale the node height such that
		}
    	return count;
    }

	private double resampleNodeHeight(NetworkNode n) {
		double logHR = 0.0;
		if (n.isLeaf()) {
			return logHR;
		} else {
			if (n.isReassortment()) {
				// check if the edge below the reassortment edge hast already been scaled
				if (!updated.get(n)) {
					double scaler = getScalerExp();
					logHR += Math.log(scaler);

					logHR += resampleNodeHeight(n.getChildEdges().get(0).childNode);				
					double newLength = oldLengths.get(n)*scaler;
					n.setHeight(newLength + n.getChildEdges().get(0).childNode.getHeight());
					updated.put(n, true);								

				}
				// update the 
			}
			else {
				double scaler = getScalerExp();
				logHR += Math.log(scaler);
				for (NetworkEdge e : n.getChildEdges()) {
					logHR += resampleNodeHeight(e.childNode);
				}
				double newLength = oldLengths.get(n)*scaler;
				n.setHeight(newLength + Math.max(n.getChildEdges().get(0).childNode.getHeight(), n.getChildEdges().get(1).childNode.getHeight()));
				updated.put(n, true);	
			}
			// scale the node height such that
		}
    	return logHR;
	}

	protected double getScalerExp() {
		return Math.exp(Randomizer.nextGaussian()*(1-scaleFactor));
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
