package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class UniformReassortmentReheight extends NetworkOperator {
	

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);
    // shadows size
    double size;
    private double limit;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
	}
    
    private double getDelta() {
        return (Randomizer.nextDouble() * size) - (size / 2.0);
    }
    

    @Override
    public double networkProposal() {

        final double delta = getDelta();


        network.startEditing(this);
        
        // pick a random coalescent node
        List<NetworkNode> node = network.getInternalNodes().
        		stream().collect(Collectors.toList());
        
        // pick a random node
        NetworkNode randomNode = node.get(Randomizer.nextInt(node.size()));
//        double logHR = Math.log(1.0/node.size());
                
        double toHeight = randomNode.getHeight() + Math.abs(2*delta);
        double fromHeight = randomNode.getHeight() - Math.abs(2*delta);
                
                        
//        System.out.println(randomNode.getHeight());
//        System.out.println(network + ";");
        scaleNodesAbove(randomNode, toHeight, fromHeight, delta);
        scaleNodesBelow(randomNode, toHeight, fromHeight, delta);
        
        for (NetworkEdge e : network.getEdges()) {
        	if (e.isRootEdge())
        		continue;
        	if (e.getLength() < 0.0) {
        		return Double.NEGATIVE_INFINITY;
            }
        }
        
//        System.out.println(network + ";");

        

        return 0.0;
    }
    
	private void scaleNodesAbove(NetworkNode n, double toHeight, double fromHeight, double delta) {	
		if (n == null) {
			return;
		}
		double currentPosition = (n.getHeight() - fromHeight) / (toHeight - fromHeight);
		if (currentPosition < 0.0 || currentPosition > 1.0) {
			return;
		}
		if (delta>0) {
//			System.out.println("currentHeight " + n.getHeight() + " new Height " + (fromHeight + delta + currentPosition * delta));
			n.setHeight(fromHeight + delta + currentPosition * delta);
		}else
			n.setHeight(fromHeight + currentPosition * Math.abs(delta));
		for (NetworkEdge edge : n.getParentEdges()) {
			scaleNodesAbove(edge.parentNode, toHeight, fromHeight, delta);
		}
    }
	
	private void scaleNodesBelow(NetworkNode n, double toHeight, double fromHeight, double delta) {
		if (n.isLeaf()) {
			return;
		}
		double currentPosition = (n.getHeight() - fromHeight) / (toHeight - fromHeight);
		if (currentPosition < 0.0 || currentPosition > 1.0) {
			return;
		}
		if (delta>0)
			n.setHeight(fromHeight + delta + currentPosition * delta);
		else
			n.setHeight(fromHeight + currentPosition  * Math.abs(delta));
		
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
            if( limit > 0 ) {
                final Network network = networkInput.get();
                final double h = network.getRootEdge().childNode.getHeight();
                final double k = Math.log(network.getLeafNodes().size()) / Math.log(2);
                final double lim = (h / k) * limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
               size = f;
            }
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
