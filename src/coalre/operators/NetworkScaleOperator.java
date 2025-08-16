package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class NetworkScaleOperator extends NetworkOperator {

	final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "magnitude factor used for scaling",
			Validate.REQUIRED);

	public Input<Boolean> scaleRootOnlyInput = new Input<>("scaleRootOnly", "Scale only the root node.", false);

	public Input<List<RealParameter>> upParametersInput = new Input<>("upParameter",
			"Parameters to scale in the SAME direction as the network.", new ArrayList<>());

	public Input<List<RealParameter>> upLogScaledParametersInput = new Input<>("upLogScaledParameter",
			"Parameters to scale in the SAME direction as the network, but with a log-scaled proposal distribution.",
			new ArrayList<>());

	public Input<List<RealParameter>> downParametersInput = new Input<>("downParameter",
			"Parameters to scale in the OPPOSITE direction as the network.", new ArrayList<>());

	public Input<List<RealParameter>> downLogScaledParametersInput = new Input<>("downLogScaledParameter",
			"Parameters to scale in the OPPOSITE direction as the network, but with a log-scaled proposal distribution.",
			new ArrayList<>());

	final public Input<Boolean> doubleDiscountInput = new Input<>("doubleDiscount",
			"flag to indicate that the scale factor of the downLogScaledParametersInput is discounted twice (default false)",
			false);

	final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);

	final public Input<Boolean> scaleIntervalsIndependently = new Input<>("scaleIntervalsIndependently",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			false);

	final public Input<Double> lastIntervalInput = new Input<>("lastInterval",
			"the last interval in which the network is scaled");

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

		double logHR = 0.0;

		final double f = getScaler();

		network.startEditing(this);

		if (scaleRootOnly) {

			oldLengths = new HashMap<>();
			updated = new HashMap<>();
			
			double oriLociMRCA = NetworkStatsLogger.getLociMRCA(network);


			for (NetworkNode node : network.getInternalNodes()) {
				if (node.isLeaf() || node.getHeight() < oriLociMRCA)
					continue;

				if (node.isReassortment())
					oldLengths.put(node, node.getChildEdges().get(0).getLength());
				else
					oldLengths.put(node, node.getHeight() - Math.max(node.getChildEdges().get(0).childNode.getHeight(),
							node.getChildEdges().get(1).childNode.getHeight()));

				updated.put(node, false);
			}
			
			return resampleNodeHeight(network.getRootEdge().childNode, oriLociMRCA);
			


		} else {
			c1++;

			oldLengths = new HashMap<>();
			updated = new HashMap<>();

			for (NetworkNode node : network.getInternalNodes()) {
				if (node.isLeaf())
					continue;

				if (node.isReassortment())
					oldLengths.put(node, node.getChildEdges().get(0).getLength());
				else
					oldLengths.put(node, node.getHeight() - Math.max(node.getChildEdges().get(0).childNode.getHeight(),
							node.getChildEdges().get(1).childNode.getHeight()));

				updated.put(node, false);
			}

			if (scaleIntervalsIndependently.get()) {
				double lhr = resampleNodeHeight(network.getRootEdge().childNode);
				return lhr;
			} else if (lastIntervalInput.get() != null) {
//				int count = resampleNodeHeightLastInterval(network.getRootEdge().childNode, lastIntervalInput.get(), f);
//				logHR += Math.log(f) * (count - 2);

				for (RealParameter param : upLogScaledParametersInput.get()) {
					param.startEditing(this);
					double logScaler = Math.log(f);
					// scale only the last interval, written to generalize in the gfuture
					for (int i = param.getDimension() - 1; i < param.getDimension(); i++) {
						param.setValue(i, param.getValue(i) + logScaler);
						logHR += logScaler;
					}
				}
				for (RealParameter param : downLogScaledParametersInput.get()) {
					param.startEditing(this);
					// convert scaler to log space
					double logScaler = doubleDiscountInput.get() ? 2 * Math.log(f) : Math.log(f);
					for (int i = param.getDimension() - 1; i < param.getDimension(); i++) {
						param.setValue(i, param.getValue(i) - logScaler);
						logHR -= doubleDiscountInput.get() ? logScaler/2 : logScaler;
					}
				}
				return logHR;

			} else {
				int count = scaleNodes(network.getRootEdge().childNode, f);
				logHR += Math.log(f) * (count - 2);
		        segmentsChanged.set(0, network.getSegmentCount(), false);

			}

		}

		// Scale parameters
		try

		{

			for (RealParameter param : upParameters) {
				param.startEditing(this);
				int count = param.scale(f);
				logHR += Math.log(f) * count;
			}

			for (RealParameter param : downParameters) {
				param.startEditing(this);
				int count = param.scale(1.0 / f);
				logHR -= Math.log(f) * count;
			}

			for (RealParameter param : upLogScaledParametersInput.get()) {
				param.startEditing(this);
				// convert scaler to log space
				double logScaler = Math.log(f);

				for (int i = 0; i < param.getDimension(); i++) {
					param.setValue(i, param.getValue(i) + logScaler);
					logHR += logScaler;
				}

			}
			for (RealParameter param : downLogScaledParametersInput.get()) {
				param.startEditing(this);
				// convert scaler to log space
				double logScaler = doubleDiscountInput.get() ? 2 * Math.log(f) : Math.log(f);
				for (int i = 0; i < param.getDimension(); i++) {
					param.setValue(i, param.getValue(i) - logScaler);
					logHR -= doubleDiscountInput.get() ? logScaler/2 : logScaler;
				}

			}

		} catch (IllegalArgumentException ex) {
			return Double.NEGATIVE_INFINITY;
		}

		return logHR;
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
					double newLength = oldLengths.get(n) * scaler;
					n.setHeight(newLength + n.getChildEdges().get(0).childNode.getHeight());
					updated.put(n, true);
					count++;

				}
				// update the
			} else {
				for (NetworkEdge e : n.getChildEdges()) {
					count += scaleNodes(e.childNode, scaler);
				}
				double newLength = oldLengths.get(n) * scaler;
				n.setHeight(newLength + Math.max(n.getChildEdges().get(0).childNode.getHeight(),
						n.getChildEdges().get(1).childNode.getHeight()));
				
                if (n.segmentIndices!=null && segmentTreesInput.get().size()>0) {
                	for (int i =0; i < n.segmentIndices.length; i++) {
                		if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i))
                			segmentTreesInput.get().get(i).getNode(n.segmentIndices[i]).setHeight(n.getHeight());
                    }
                }
                
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
					double newLength = oldLengths.get(n) * scaler;
					n.setHeight(newLength + n.getChildEdges().get(0).childNode.getHeight());
					
			        if (n.isCoalescence()) {
				        if (n.segmentIndices!=null && segmentTreesInput.get().size()>0) {
				        	for (int i =0; i < n.segmentIndices.length; i++) {
				        		if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i))
				        			segmentTreesInput.get().get(i).getNode(n.segmentIndices[i]).setHeight(n.getHeight());
							}
				        }
			        }

					updated.put(n, true);

				}
				// update the
			} else {
				double scaler = getScalerExp();
				logHR += Math.log(scaler);
				for (NetworkEdge e : n.getChildEdges()) {
					logHR += resampleNodeHeight(e.childNode);
				}
				double newLength = oldLengths.get(n) * scaler;
				n.setHeight(newLength + Math.max(n.getChildEdges().get(0).childNode.getHeight(),
						n.getChildEdges().get(1).childNode.getHeight()));
				
		        if (n.isCoalescence()) {
			        if (n.segmentIndices!=null && segmentTreesInput.get().size()>0) {
			        	for (int i =0; i < n.segmentIndices.length; i++) {
			        		if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i))
			        			segmentTreesInput.get().get(i).getNode(n.segmentIndices[i]).setHeight(n.getHeight());
						}
			        }
		        }

				updated.put(n, true);
			}
			// scale the node height such that
		}
		return logHR;
	}

	private double resampleNodeHeight(NetworkNode n, double minHeight) {
		double logHR = 0.0;
		
		if (n.getHeight() < minHeight) {
			return logHR;
		} else {
			if (n.isReassortment()) {
				// check if the edge below the reassortment edge hast already been scaled
				if (!updated.get(n)) {
					double scaler = getScalerExp();
					logHR += Math.log(scaler);

					logHR += resampleNodeHeight(n.getChildEdges().get(0).childNode, minHeight);
					double newLength = oldLengths.get(n) * scaler;
					n.setHeight(newLength + n.getChildEdges().get(0).childNode.getHeight());
					
			        if (n.isCoalescence()) {
				        if (n.segmentIndices!=null && segmentTreesInput.get().size()>0) {
				        	for (int i =0; i < n.segmentIndices.length; i++) {
				        		if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i))
				        			segmentTreesInput.get().get(i).getNode(n.segmentIndices[i]).setHeight(n.getHeight());
							}
				        }
			        }

					updated.put(n, true);

				}
				// update the
			} else {
				double scaler = getScalerExp();
				logHR += Math.log(scaler);
				for (NetworkEdge e : n.getChildEdges()) {
					logHR += resampleNodeHeight(e.childNode, minHeight);
				}
				double newLength = oldLengths.get(n) * scaler;
				n.setHeight(newLength + Math.max(n.getChildEdges().get(0).childNode.getHeight(),
						n.getChildEdges().get(1).childNode.getHeight()));
				
		        if (n.isCoalescence()) {
			        if (n.segmentIndices!=null && segmentTreesInput.get().size()>0) {
			        	for (int i =0; i < n.segmentIndices.length; i++) {
			        		if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i))
			        			segmentTreesInput.get().get(i).getNode(n.segmentIndices[i]).setHeight(n.getHeight());
						}
			        }
		        }

				updated.put(n, true);
			}
			// scale the node height such that
		}
		return logHR;
	}


	private int resampleNodeHeightLastInterval(NetworkNode n, double lastInterval, double scaler) {
		int count = 0;

		
		if (n.getHeight() < lastInterval) {
			// if the node is below the last interval, just return 0
			return 0;
		}

		if (n.isLeaf()) {
			return 0;
		} else {
			if (n.isReassortment()) {
				// check if the edge below the reassortment edge hast already been scaled
				if (!updated.get(n)) {

					count += resampleNodeHeightLastInterval(n.getChildEdges().get(0).childNode, lastInterval, scaler);
                    double oldLength = oldLengths.get(n);
					// adjust the scaler such that only the part above lastInterval is scaled
					if (n.getChildEdges().get(0).childNode.getHeight() < lastInterval) {
						// discount the part that is below the interval
						oldLength -= (lastInterval - n.getChildEdges().get(0).childNode.getHeight());
					}
					count++;
					
					double newLength = oldLength * scaler;
					if (n.getChildEdges().get(0).childNode.getHeight() < lastInterval) {
						n.setHeight(newLength + lastInterval);
					}else {
						n.setHeight(newLength + n.getChildEdges().get(0).childNode.getHeight());
					}
					updated.put(n, true);

				}
				// update the
			} else {

				for (NetworkEdge e : n.getChildEdges()) {
					count += resampleNodeHeightLastInterval(e.childNode, lastInterval, scaler);
				}
				// adjust the scaler such that only the part above lastInterval is scaled
				// adjust the scaler such that only the part above lastInterval is scaled
                double oldLength = oldLengths.get(n);

				if (Math.max(n.getChildEdges().get(0).childNode.getHeight(),
						n.getChildEdges().get(1).childNode.getHeight()) < lastInterval) {
					// discount the part that is below the interval
					oldLength -= (lastInterval - Math.max(n.getChildEdges().get(0).childNode.getHeight(),
							n.getChildEdges().get(1).childNode.getHeight()));
				}

				count++;
				double newLength = oldLength * scaler;
				
				if (Math.max(n.getChildEdges().get(0).childNode.getHeight(),
						n.getChildEdges().get(1).childNode.getHeight()) < lastInterval) {
					n.setHeight(newLength + lastInterval);
				}else {
					n.setHeight(newLength + Math.max(n.getChildEdges().get(0).childNode.getHeight(),
							n.getChildEdges().get(1).childNode.getHeight()));
				}

				
		        if (n.segmentIndices!=null && segmentTreesInput.get().size()>0) {
		        	for (int i =0; i < n.segmentIndices.length; i++) {
		        		if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i))
		        			segmentTreesInput.get().get(i).getNode(n.segmentIndices[i]).setHeight(n.getHeight());
					}
		        }

				updated.put(n, true);
			}
			// scale the node height such that
		}
		return count;
	}

	protected double getScalerExp() {
		return Math.exp(Randomizer.nextGaussian() * (1 - scaleFactor));
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
		if (ratio > 2.0)
			ratio = 2.0;
		if (ratio < 0.5)
			ratio = 0.5;

		// new scale factor
		final double sf = Math.pow(scaleFactor, ratio);

		final DecimalFormat formatter = new DecimalFormat("#.###");
		if (prob < 0.10) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else if (prob > 0.40) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else
			return "";
	}

}
