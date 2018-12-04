package coalre.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.sun.xml.internal.bind.v2.runtime.output.StAXExStreamWriterOutput;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import jdk.nashorn.internal.runtime.regexp.joni.Warnings;


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
public class SubNetworkSlide extends DivertSegmentOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
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

	@Override
	public double networkProposal() {

		double logHR = 0.0;
		network.startEditing(this);
		// 1. Choose a random edge, avoiding root
		List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
		
		final List<NetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> !e.parentNode.isReassortment())
				.collect(Collectors.toList());
		
		final int possibleNodes = possibleEdges.size();
        if (possibleNodes < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        logHR -= Math.log(1.0/possibleNodes);        

		final NetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleNodes));
		final NetworkNode i = iEdge.childNode;
		final NetworkNode ip = iEdge.parentNode;
		final NetworkEdge ipEdge = ip.getParentEdges().get(0); //ip not a reassortment, has only one parent
		final NetworkEdge jEdge = getSisterEdge(iEdge);
		final NetworkNode j = jEdge.childNode;
		final NetworkNode pip = ipEdge.parentNode; 
		NetworkEdge newChildEdge = null;
		
		
        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = ip.getHeight();
        final double newHeight = oldHeight + delta;
        

        
        // 3. if the move is up
        if (delta > 0) {

        	// 3.1 if the topology will change
            if (pip != null &&  pip.getHeight() < newHeight) {

                NetworkNode newParent = pip;
                newChildEdge = ipEdge;
                while (newParent.getHeight() < newHeight) {
                	
                    if (newParent.isReassortment()) {
                    	if (Randomizer.nextBoolean()) {
                    		newChildEdge = newParent.getParentEdges().get(0);
                    		newParent = newChildEdge.parentNode;
                    	}
                    	else {
                    		newChildEdge = newParent.getParentEdges().get(1);
                    		newParent = newChildEdge.parentNode;
                    	}
                    	logHR -= Math.log(0.5);
                    } else {
                        // if not reassortment, has only one parent
                        newChildEdge = newParent.getParentEdges().get(0);
                        newParent = newChildEdge.parentNode;
                    }

                    if (newParent == null) break;
                }
                // the moved node 'p' would become a child of 'newParent'

                // 3.1.1 if creating a new root
                if (newChildEdge.isRootEdge()) {

                	NetworkNode destNode = new NetworkNode();
                	destNode.setHeight(newHeight);
                	
                	NetworkEdge newIpEdge = new NetworkEdge();
                	newIpEdge.hasSegments = (BitSet)newChildEdge.hasSegments.clone();
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);
                	
                	NetworkEdge iEdgeCopy = new NetworkEdge();
                	iEdgeCopy.hasSegments = (BitSet)iEdge.hasSegments.clone();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
                	
                	BitSet segsToAdd = (BitSet)iEdgeCopy.hasSegments.clone();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).hasSegments);
                	logHR -= addSegmentsToAncestors(newIpEdge, segsToAdd);
                	
                	BitSet segsToRemove = (BitSet)iEdge.hasSegments.clone();
                	segsToRemove.andNot(getSisterEdge(iEdge).hasSegments);
                	logHR += removeSegmentsFromAncestors(ipEdge, segsToRemove);
                	
                	
                	pip.removeChildEdge(ipEdge);
                	ip.removeChildEdge(jEdge);
                	pip.addChildEdge(jEdge);
                	
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);
                	ip.removeParentEdge(ipEdge);
                	
                	
                	newIpEdge.parentNode = null;
                	network.setRootEdge(newIpEdge);
                	

                }
                // 3.1.2 no new root
                else {
                	NetworkNode destNode = new NetworkNode();
                	destNode.setHeight(newHeight);
                	
                	NetworkEdge newIpEdge = new NetworkEdge();
                	newIpEdge.hasSegments = (BitSet)newChildEdge.hasSegments.clone();
                	newParent.removeChildEdge(newChildEdge);
                	newParent.addChildEdge(newIpEdge);
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);
                	
                	NetworkEdge iEdgeCopy = new NetworkEdge();
                	iEdgeCopy.hasSegments = (BitSet)iEdge.hasSegments.clone();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
          	
                	BitSet segsToAdd = (BitSet)iEdgeCopy.hasSegments.clone();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).hasSegments);
                	logHR -= addSegmentsToAncestors(newIpEdge, segsToAdd);

                	
                	BitSet segsToRemove = (BitSet)iEdge.hasSegments.clone();
                	segsToRemove.andNot(getSisterEdge(iEdge).hasSegments);
                	logHR += removeSegmentsFromAncestors(ipEdge, segsToRemove);
                	
                	pip.removeChildEdge(ipEdge);
                	ip.removeChildEdge(jEdge);
                	pip.addChildEdge(jEdge);
                	
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);
                	ip.removeParentEdge(ipEdge);
                }
                
                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChildEdge, oldHeight, null);

                logHR += Math.log(1.0/possibleSources);

            } else {
                ip.setHeight(newHeight);
            }
        }
        
        // 4 if we are sliding the subnet down.
        else {
        	if (i.getHeight() > newHeight) {
        		return Double.NEGATIVE_INFINITY;
        	}
  	
        	
            // 4.1 will the move change the topology
            if (j.getHeight() > newHeight) {	

            	final List<NetworkEdge> possibleChildEdges = new ArrayList<>();
            	final int possibleDestinations = intersectingEdges(jEdge, newHeight, possibleChildEdges);
            	
                // if no valid destinations then return a failure
                if (possibleChildEdges.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }
                
                logHR -= Math.log(1.0/possibleDestinations); 
                
                // pick a random parent/child destination edge uniformly from options
                final int childEdgeIndex = Randomizer.nextInt(possibleChildEdges.size());
                newChildEdge = possibleChildEdges.get(childEdgeIndex);
                NetworkNode newParent = newChildEdge.parentNode;
            	
                
                // 4.1.1 if p was root
                if (ipEdge.isRootEdge()) {
                	NetworkNode destNode = new NetworkNode();
                	destNode.setHeight(newHeight);
                	
                   	NetworkEdge newIpEdge = new NetworkEdge();
                	newIpEdge.hasSegments = (BitSet)newChildEdge.hasSegments.clone();
                	newParent.removeChildEdge(newChildEdge);
                	newParent.addChildEdge(newIpEdge);
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);
                	
                	
                	NetworkEdge iEdgeCopy = new NetworkEdge();
                	iEdgeCopy.hasSegments = (BitSet)iEdge.hasSegments.clone();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
                	BitSet segsToAdd = (BitSet)iEdgeCopy.hasSegments.clone();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).hasSegments);
                	logHR -= addSegmentsToAncestors(newIpEdge, segsToAdd);
                	
                	BitSet segsToRemove = (BitSet)iEdge.hasSegments.clone();
                	segsToRemove.andNot(getSisterEdge(iEdge).hasSegments);
                	logHR += removeSegmentsFromAncestors(ipEdge, segsToRemove);
                	
                	ip.removeChildEdge(jEdge);
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);
                	ip.removeParentEdge(ipEdge);
                	
                	jEdge.parentNode = null;
                	network.setRootEdge(jEdge);

                } else {
                	NetworkNode destNode = new NetworkNode();
                	destNode.setHeight(newHeight);
                	
                	NetworkEdge newIpEdge = new NetworkEdge();
                	newIpEdge.hasSegments = (BitSet)newChildEdge.hasSegments.clone();
                	newParent.removeChildEdge(newChildEdge);
                	newParent.addChildEdge(newIpEdge);
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);

                	NetworkEdge iEdgeCopy = new NetworkEdge();
                	iEdgeCopy.hasSegments = (BitSet)iEdge.hasSegments.clone();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
                	BitSet segsToAdd = (BitSet)iEdgeCopy.hasSegments.clone();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).hasSegments);
                	logHR -= addSegmentsToAncestors(newIpEdge, segsToAdd);
                	
                	BitSet segsToRemove = (BitSet)iEdge.hasSegments.clone();
                	segsToRemove.andNot(getSisterEdge(iEdge).hasSegments);
                	logHR += removeSegmentsFromAncestors(ipEdge, segsToRemove);
                	
                	pip.removeChildEdge(ipEdge);
                	ip.removeChildEdge(jEdge);
                	pip.addChildEdge(jEdge);
                	
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);

                }
                
                ip.setHeight(newHeight);
            	
            } else {
            	
                ip.setHeight(newHeight);
            }
        }
        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;
        

        
		List<NetworkEdge> networkEdgesAfter = new ArrayList<>(network.getEdges());
		
		final List<NetworkEdge> possibleEdgesAfter = networkEdgesAfter.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> !e.parentNode.isReassortment())
				.collect(Collectors.toList());
		
		final int possibleNodesAfter = possibleEdgesAfter.size();
        if (possibleNodesAfter < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        
        logHR += Math.log(1.0/possibleNodesAfter);  

		return logHR;
	}
	
	 
    public double checkAncestralSegments(NetworkEdge edge) {
    	double logHR = 0.0;
    	
    	NetworkNode parentNode = edge.parentNode;
    	if (parentNode == null) {
    		return logHR;
    	}
    	
    	if (parentNode.getParentCount() > 1) {
    		NetworkEdge parentEdge = parentNode.getParentEdges().get(0);
    		NetworkEdge spouseEdge = parentNode.getParentEdges().get(1);
    		
    		BitSet segsToAdd = (BitSet)edge.hasSegments.clone();
    		segsToAdd.andNot(parentEdge.hasSegments);
    		segsToAdd.andNot(spouseEdge.hasSegments);
    		
    		if (!segsToAdd.isEmpty()) {
    			if (Randomizer.nextBoolean()) {
    				logHR -= addSegmentsToAncestors(parentEdge, segsToAdd);
    			} else {
    				logHR -= addSegmentsToAncestors(spouseEdge, segsToAdd);
    			}
    		}
    		
    		BitSet segsToRemoveParent = (BitSet)parentEdge.hasSegments.clone();
    		segsToRemoveParent.andNot(edge.hasSegments);
    		if (!segsToRemoveParent.isEmpty()) {
    			logHR += removeSegmentsFromAncestors(parentEdge, segsToRemoveParent);
    		}
    		
    		BitSet segsToRemoveSpouse= (BitSet)spouseEdge.hasSegments.clone();
    		segsToRemoveSpouse.andNot(edge.hasSegments);
    		if (!segsToRemoveSpouse.isEmpty()) {
    			logHR += removeSegmentsFromAncestors(spouseEdge, segsToRemoveSpouse);
    		}
	
    	} else {
    		NetworkEdge parentEdge = parentNode.getParentEdges().get(0);
    		BitSet segsToAdd = (BitSet)edge.hasSegments.clone();
    		segsToAdd.andNot(parentEdge.hasSegments);
    		
    		if (!segsToAdd.isEmpty()) {
    			logHR -= addSegmentsToAncestors(parentEdge, segsToAdd);
    		}
    		
    		BitSet segsToRemoveParent = (BitSet)parentEdge.hasSegments.clone();
    		segsToRemoveParent.andNot(edge.hasSegments);
    		if (!segsToRemoveParent.isEmpty()) {
    			logHR += removeSegmentsFromAncestors(parentEdge, segsToRemoveParent);
    		}
    	}
    	return logHR;
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


    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }
    
    private int intersectingEdges(NetworkEdge edge, double height, List<NetworkEdge> directChildEdges) {
        final NetworkNode node = edge.childNode;
    	final NetworkNode parent = edge.parentNode;
    	
        if (parent.getHeight() < height) return 0;


        if (node.getHeight() < height) {
            if (directChildEdges != null) directChildEdges.add(edge);
            return 1;
        }

        if (node.getChildCount() > 1 && node.getChildEdges().get(0).childNode == node.getChildEdges().get(1).childNode) {
        	int count = 0;
        	if (node.getChildEdges().get(0).childNode.getHeight() < height) {
        		if (directChildEdges != null) {
        			directChildEdges.add(node.getChildEdges().get(0));
        			directChildEdges.add(node.getChildEdges().get(1));
        		}
        		return 2;
        	}

        	count += intersectingEdges(node.getChildEdges().get(0).childNode.getChildEdges().get(0), height, directChildEdges);        	
        	return count;
        }
        if (node.isLeaf()) {
            // TODO: verify that this makes sense
            return 0;
        } else {
            int count = intersectingEdges(node.getChildEdges().get(0), height, directChildEdges);
            if (node.getChildCount() > 1) count += 
            		intersectingEdges(node.getChildEdges().get(1), height, directChildEdges);          
            return count;
        }
    }  
    
    /**
     * Simple (but probably too expensive) check for a kind of invalid network
     * which can result from an edge deletion operation: one in which the
     * network posesses nontrivial structure above the MRCA. (I.e. the MRCA
     * is not the root.)
     *
     * @return true if the network terminates at the true MRCA. (I.e. is valid.)
     */
    protected boolean networkTerminatesAtMRCA() {
        List<NetworkNode> sortedNodes = new ArrayList<>(network.getNodes());
        sortedNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));
        List<NetworkNode> sampleNodes = sortedNodes.stream().filter(NetworkNode::isLeaf).collect(Collectors.toList());
        double maxSampleHeight = sampleNodes.get(sampleNodes.size()-1).getHeight();

        int lineages = 0;
        for (NetworkNode node : sortedNodes) {
            switch(node.getChildEdges().size()) {
                case 2:
                    // Coalescence

                    lineages -= 1;
                    break;

                case 1:
                    // Reassortment

                    if (lineages < 2 && node.getHeight() > maxSampleHeight)
                        return false;

                    lineages += 1;
                    break;

                case 0:
                    // Sample

                    lineages += 1;
                    break;
            }
        }

        return true;
    }
}
