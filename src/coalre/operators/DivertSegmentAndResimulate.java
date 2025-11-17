package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class DivertSegmentAndResimulate extends NetworkOperator {

	public Input<Boolean> divertOneSegmentInput = new Input<>("divertOneSegment",
			"If true, only one segment is diverted", true);

	public Input<CoalescentWithReassortment> coalescentDistrInput = new Input<>("coalescentWithReassortment",
			"Coalescent distribution for sampling reassortment events.",
			Input.Validate.OPTIONAL);

	protected CoalescentWithReassortment coalescentDistr;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		// set the value of lambda to be 0 in the super class
		
		
		if (coalescentDistrInput.get() != null) {
			coalescentDistr = coalescentDistrInput.get();
		}
		// if the super class has a lambda parameter, set its value to zero
	}
	
	boolean exitOnError = false;
	String networkBefore ="";
	int eventsBefore = 0;
	int eventsAfter = 0;
	
	List<NetworkEdge> networkEdges;
	
	@Override
    public double proposal() {
		double logHR = 0.0;
		
		

		if (network.segmentTreesNeedInit) {			
			System.out.println("restore segment trees in first move");
			// update all segment trees			
            for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
          		network.updateSegmentTree(segmentTrees.get(segmentTreeMap[segIdx]), segmentTreeMap[segIdx]); 
            
            network.segmentTreesNeedInit = false;
			return Double.POSITIVE_INFINITY;
		}
		
		network.startEditing(this);
		networkEdges = new ArrayList<>(network.getEdges());
		
//		int edgesBefore = networkEdges.size();
//		String networkBefore = network.getExtendedNewick();
		logHR = networkProposal();
		
//		int edgesAfter = network.getEdges().size();
//		System.out.println((edgesAfter-edgesBefore) + " " + logHR + " " + getID());
		if (Double.isNaN(logHR)) {
			throw new IllegalArgumentException("NaN in proposal logHR");
		}
		if (logHR ==Double.POSITIVE_INFINITY)
			throw new IllegalArgumentException("Infinity in proposal logHR");
		// check if any edge is empty
//		for (NetworkEdge edge : network.getEdges()) {
//			if (edge.hasSegments.isEmpty()) {
//				throw new IllegalArgumentException("Empty edge created!");
//			}
//		}
//		if (Math.abs(logHR)>40)
//			System.out.println("DivertSegmentAndResimulate: " + logHR);
		return logHR;
    }
	
	@Override
	public double networkProposal() {
		double logHR = 0.0;
		
		List<Integer> sourceEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			if (networkEdges.get(i).childNode.isReassortment() && networkEdges.get(i).hasSegments.cardinality()>1) {
				sourceEdges.add(i);
			}
		}

		if (sourceEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
		

		logHR -= Math.log(1.0 / sourceEdges.size());

		NetworkEdge sourceEdge = networkEdges.get(sourceEdges.get(Randomizer.nextInt(sourceEdges.size())));
		NetworkEdge destEdge = getSpouseEdge(sourceEdge);

		BitSet segsToDivert;

		if (divertOneSegmentInput.get()) {
			segsToDivert = getRandomSegment(sourceEdge.hasSegments);
			logHR -= Math.log(1. / sourceEdge.hasSegments.cardinality());
		} else {
			segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
			logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, 0.5);
		}

		if (segsToDivert.cardinality() == 0)
			return Double.NEGATIVE_INFINITY;

		double logDivert = divertSegments(destEdge, sourceEdge, segsToDivert);
		
		logHR += logDivert;


		if (divertOneSegmentInput.get()) {
			logHR += Math.log(1.0 / destEdge.hasSegments.cardinality());
		} else {
//			System.err.println("this part of the operator has an error");
			logHR += getLogConditionedSubsetProb(destEdge.hasSegments, segsToDivert, 0.5);
		}
		
		
		int reverseSourceEdgeCount = 0;
		for (NetworkEdge edge : networkEdges){
			if (edge.childNode.isReassortment() && edge.hasSegments.cardinality()>1){
                reverseSourceEdgeCount++;
            }
		}
		logHR += Math.log(1.0 / reverseSourceEdgeCount);

		return logHR;
	}
	

	protected double divertSegments(NetworkEdge destEdge, NetworkEdge sourceEdge, BitSet segsToDivert) {
		double logHR = 0.0;
		
//		System.out.println("------- ");
//		System.out.println("------- ");
//		System.out.println("------- ");
//		System.out.println("------- ");
//		
//		System.out.println(network);

		// Get network events to determine lineages at different times, they will contain any inactive Edges, but not newly added active ones
		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];

		getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);
		exitOnError=false;
		
		BitSet segsToDivertCopy = (BitSet) segsToDivert.clone();

		// make a defensive implementation by selection one segment at a time to divert, choose the order randomly
		for (int i = 0; i < segsToDivert.cardinality(); i++) {
			// select a random segment from segsToDivertCopy
			int segIdx = -1;
			int randIdx = Randomizer.nextInt(segsToDivertCopy.cardinality());
			for (int j = 0; j <= randIdx; j++) {
				segIdx = segsToDivertCopy.nextSetBit(segIdx + 1);
			}
			segsToDivertCopy.set(segIdx, false);
			
			BitSet currentSegsToDivert = new BitSet();
			currentSegsToDivert.set(segIdx);
			
			BitSet segsToAdd = (BitSet) currentSegsToDivert.clone();
			// startwith keeping a list of active edges and the segs to Add
			List<NetworkEdge> activeEdges = new ArrayList<>();
			List<BitSet> segsToAddList = new ArrayList<>();
			
			List<NetworkEdge> inactiveEdges = new ArrayList<>();
			List<BitSet> segsToRemoveList = new ArrayList<>();
			
			activeEdges.add(destEdge);
			segsToAddList.add(segsToAdd);
			
			inactiveEdges.add(sourceEdge);
			segsToRemoveList.add((BitSet) currentSegsToDivert.clone());
			
			coalescentDistr.intervals.update();
			List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();

			
//			System.out.println("Diverting segment " + segIdx + " from edge " + sourceEdge + " to edge " + destEdge);
//			System.out.println(network);
	
			// Track the current time (start at the bottom of the lowest edge)
			double currentTime = destEdge.childNode.getHeight();
			logHR += divertSegmentsToAncestors(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, currentTime, networkEventList, false, false);
			
//			System.out.println("After diverting segment " + network);
			if (reconnectSegmentTrees(treeChildNodeList, destEdge, currentSegsToDivert))
				return Double.NEGATIVE_INFINITY;
			
//			System.out.println(network.getExtendedNewick(0));
			
			cleanEmptyEdgesTopDown();

		}
//		System.out.println(network);
//		System.out.println("Diverted segments " + segsToDivert);

//		if (segsToDivert.cardinality() > 1)
//			System.exit(0);
		
		
		
		return logHR;
//		return (logPbefore - logPafter);
	}

	protected boolean reconnectSegmentTrees(Integer[] treeChildNodeList, NetworkEdge destEdge, BitSet segsToDivert) {

		List<NetworkNode> targetNodeIndices = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++)
			targetNodeIndices.add(null);
		Integer[] newTargetNode = new Integer[network.getSegmentCount()];
		Double[] newNodeHeights = new Double[network.getSegmentCount()];

		getTreeNodes(destEdge, segsToDivert, targetNodeIndices, newTargetNode, newNodeHeights, treeChildNodeList);

		// update the corresponding tree
		for (int i = 0; i < network.getSegmentCount(); i++) {

			if (treeChildNodeList[i] != null && newTargetNode[i] != null) {
				Node startChild = segmentTrees.get(i).getNode(treeChildNodeList[i]);
				if (startChild.isRoot())
					continue;
				
				Node oldParent = startChild.getParent();
				Node newSibling = segmentTrees.get(i).getNode(newTargetNode[i]);

				if (oldParent.isRoot()) {

					if (oldParent==newSibling) {
						oldParent.setHeight(newNodeHeights[i]);
					}else if (newSibling.getParent()==oldParent) {
						oldParent.setHeight(newNodeHeights[i]);
					}else {

						// the node that will become the new root
						Node oldSibling = oldParent.getChild(0).getNr() == treeChildNodeList[i] ? oldParent.getChild(1) : oldParent.getChild(0);
						
						int prevNumberNewRoot = oldSibling.getNr();
						int rootNumber = segmentTrees.get(i).getRoot().getNr();
						double newRootHeight = oldSibling.getHeight();
						
						// update the segmentIndices of the network node
						if (targetNodeIndices.get(i).segmentIndices == null)
							targetNodeIndices.get(i).segmentIndices = new int[network.getSegmentCount()];
						targetNodeIndices.get(i).segmentIndices[i] = prevNumberNewRoot;

						updatePreviousRootNumber(i, network.getRootEdge().childNode, rootNumber, newRootHeight);
						
						Node newGrandParent = newSibling.getParent();
						newGrandParent.removeChild(newSibling);
						newGrandParent.addChild(oldParent);
						oldParent.addChild(newSibling);
						// do the corresponding tree move
						oldParent.removeChild(oldSibling);
						oldParent.setHeight(newNodeHeights[i]);
						oldSibling.setParent(null);
						
						segmentTrees.get(i).setRoot(oldSibling);
					}
					
				}else {
					Node oldGrandParent = oldParent.getParent();
					Node oldSibling = oldParent.getChild(0).getNr() == treeChildNodeList[i] ? oldParent.getChild(1) : oldParent.getChild(0);
					Node newGrandParent = newSibling.getParent();
					
					if (newGrandParent==oldParent) {
						oldParent.setHeight(newNodeHeights[i]);
					}else if (oldGrandParent==newGrandParent) {
						if (oldParent==newSibling) {
							throw new IllegalArgumentException("oldParent==newSibling " + oldParent.getHeight() + " " + newSibling.getHeight());
						}
						
						newGrandParent.removeChild(newSibling);
						oldParent.removeChild(oldSibling);
						newGrandParent.addChild(oldSibling);
						
						oldParent.addChild(newSibling);
						oldParent.setHeight(newNodeHeights[i]);
						
					}else {
						
						oldGrandParent.removeChild(oldParent);
						oldGrandParent.addChild(oldSibling);
						if (newGrandParent==null) {
							// changing the root node also changes the node number
							// which is an issue here 
							oldParent.setParent(null);
							
							int prevRootNumber = oldParent.getNr();
							double oldRootHeight = segmentTrees.get(i).getRoot().getHeight();
							
							// find the network node corresponding to the old root height
							// and change the corresponding segment index
							updatePreviousRootNumber(i,targetNodeIndices.get(i), prevRootNumber, oldRootHeight);
													oldParent.removeChild(oldSibling);
							oldParent.addChild(newSibling);
							oldParent.setHeight(newNodeHeights[i]);

							segmentTrees.get(i).setRoot(oldParent);
							
							// update the segmentIndices of the network node
							if (targetNodeIndices.get(i).segmentIndices == null)
								targetNodeIndices.get(i).segmentIndices = new int[network.getSegmentCount()];
							targetNodeIndices.get(i).segmentIndices[i] = oldParent.getNr();
	
						}else {
							newGrandParent.removeChild(newSibling);
							oldParent.removeChild(oldSibling);
							newGrandParent.addChild(oldParent);
							oldParent.addChild(newSibling);
							oldParent.setHeight(newNodeHeights[i]);
						}
					}
				}
			}

		}	
		return false;
	}

	private void updatePreviousRootNumber(int i, NetworkNode n, int prevRootNumber, double oldRootHeight) {
		if (n.getHeight()<oldRootHeight) {
			return;
		}
		
		if (n.isCoalescence()) {
			if (n.getHeight() == oldRootHeight) {
				// update the segment index
				if (n.segmentIndices == null)
					n.segmentIndices = new int[network.getSegmentCount()];
				n.segmentIndices[i] = prevRootNumber;
				return;
			}
			
			if (n.getChildEdges().get(0).hasSegments.get(i))
				updatePreviousRootNumber(i, n.getChildEdges().get(0).childNode, prevRootNumber, oldRootHeight);
			if (n.getChildEdges().get(1).hasSegments.get(i))
				updatePreviousRootNumber(i, n.getChildEdges().get(1).childNode, prevRootNumber, oldRootHeight);
		} else if (n.isReassortment()) {
			updatePreviousRootNumber(i, n.getChildEdges().get(0).childNode, prevRootNumber, oldRootHeight);
		}
		return;
	}

	void getTreeNodesDown(NetworkEdge edge, BitSet segsToRemove, Integer[] treeNodeList) {
		if (segmentTrees.size() == 0) {
			return;
		}
			
		segsToRemove = (BitSet) segsToRemove.clone();
		segsToRemove.and(edge.hasSegments);

		if (segsToRemove.isEmpty())
			return;
		
		if (edge.childNode.isReassortment()) {
			getTreeNodesDown(edge.childNode.getChildEdges().get(0), segsToRemove, treeNodeList);
		} else if (edge.childNode.isCoalescence()) { // check if this node corresponds to the tree node
			if (segmentTrees.size() > 0 && edge.childNode.segmentIndices != null) {
				for (int segIdx = 0; segIdx < segsToRemove.length(); segIdx++) {
					// get the corresponding tree node
					if (segsToRemove.get(segIdx) && 
							edge.childNode.getChildEdges().get(0).hasSegments.get(segIdx) && 
							edge.childNode.getChildEdges().get(1).hasSegments.get(segIdx)) {
						
						treeNodeList[segIdx] = edge.childNode.segmentIndices[segIdx];

						// set bitset to false
						segsToRemove.set(segIdx, false);
					}
				}
			}

			getTreeNodesDown(edge.childNode.getChildEdges().get(0), segsToRemove, treeNodeList);
			getTreeNodesDown(edge.childNode.getChildEdges().get(1), segsToRemove, treeNodeList);

		} else {
			// loop over the segsToRemove Left and update the Tree Nodes
			for (int segIdx = 0; segIdx < segsToRemove.length(); segIdx++) {
				if (segsToRemove.get(segIdx)) {
					treeNodeList[segIdx] = edge.childNode.segmentIndices[segIdx];
				}
			}
		}

		return;
	}

	void getTreeNodes(NetworkEdge edge, BitSet segsToRemove, List<NetworkNode> targetNodeIndices,
			Integer[] treeNodeList, Double[] nodeHeightList, Integer[] treeChildNodeList) {
		segsToRemove = (BitSet) segsToRemove.clone();
		segsToRemove.and(edge.hasSegments);

		if (segsToRemove.isEmpty())
			return;

		if (edge.isRootEdge())
			return;

		if (edge.parentNode.isReassortment()) {

			getTreeNodes(edge.parentNode.getParentEdges().get(0), segsToRemove, targetNodeIndices, treeNodeList,
					nodeHeightList, treeChildNodeList);
			getTreeNodes(edge.parentNode.getParentEdges().get(1), segsToRemove, targetNodeIndices, treeNodeList,
					nodeHeightList, treeChildNodeList);

		} else {
			if (segmentTrees.size() > 0) {
				BitSet sibSegs = getSisterEdge(edge).hasSegments;

				if (segmentTrees.size() > 0) {
					for (int segIdx = 0; segIdx < segsToRemove.length(); segIdx++) {
						
						if (segsToRemove.get(segIdx) && sibSegs.get(segIdx)) {
							nodeHeightList[segIdx] = edge.parentNode.getHeight();
							treeNodeList[segIdx] = getTreeNodeIndex(getSisterEdge(edge), segIdx);
							targetNodeIndices.set(segIdx, edge.parentNode);
							
							if (edge.parentNode.segmentIndices == null)
								edge.parentNode.segmentIndices = new int[network.getSegmentCount()];
							try {
							edge.parentNode.segmentIndices[segIdx] = segmentTrees.get(segIdx)
									.getNode(treeChildNodeList[segIdx]).getParent().getNr();
							}catch (Exception e) {
								System.out.println(segmentTrees.get(segIdx)
									.getNode(treeChildNodeList[segIdx]).getHeight());
								System.out.println(treeChildNodeList[segIdx]);
								System.out.println("Error in segment tree " + segIdx + " " + edge.parentNode.getHeight());
								System.out.println(network.getExtendedNewick(segIdx));
								System.out.println(segmentTrees.get(segIdx) + ";");
								System.exit(0);
							}
							segsToRemove.set(segIdx, false);
						}
					}
				}

			}

			getTreeNodes(edge.parentNode.getParentEdges().get(0), segsToRemove, targetNodeIndices, treeNodeList,
					nodeHeightList, treeChildNodeList);

		}

		return;
	}

	private Integer getTreeNodeIndex(NetworkEdge edge, int segIdx) {
		if (edge.childNode.isLeaf()) {
			return edge.childNode.segmentIndices[segIdx];
		} else if (edge.childNode.isReassortment()) {
			return getTreeNodeIndex(edge.childNode.getChildEdges().get(0), segIdx);
		} else {
			if (edge.childNode.segmentIndices != null && 
					edge.childNode.getChildEdges().get(0).hasSegments.get(segIdx) && 
					edge.childNode.getChildEdges().get(1).hasSegments.get(segIdx)) {
				return edge.childNode.segmentIndices[segIdx];
			}
			return edge.childNode.getChildEdges().get(0).hasSegments.get(segIdx)
					? getTreeNodeIndex(edge.childNode.getChildEdges().get(0), segIdx)
					: getTreeNodeIndex(edge.childNode.getChildEdges().get(1), segIdx);
		}
	}

	int count=0;
	double val = 0.0;
	double val2 = 0.0;
	/**
	 * Add segments to this edge and ancestors, potentially sampling new reassortment events.
	 * @param c 
	 * @param b 
	 *
	 * @param edge      edge at which to start addition
	 * @param segsToAdd segments to add to the edge and ancestors
	 * @param eventTimes 
	 * @return log probability of operation
	 */
	protected double divertSegmentsToAncestors(List<NetworkEdge> activeEdges, 
			List<NetworkEdge> inactiveEdges, 
			List<BitSet> segsToAddList, 
			List<BitSet> segsToRemoveList, 
			double currentTime, List<NetworkEvent> networkEventList, 
			boolean zeroActiveReassortment, boolean zeroInactiveReassortment) {
		double logHR = 0.0;
		
		NetworkEdge rootEdge = network.getRootEdge();
		List<NetworkEdge> edgesAdded = new ArrayList<>();
//		System.out.println(network.getExtendedNewick(0));
				
		// now, sample the time to the next reassortment or coalescent event on the active edges,
		// i.e. essentially simulate the history of these edges and segments until there is nothing left. The
		// goal is to avoid the creation of empty
		while (!activeEdges.isEmpty() || !inactiveEdges.isEmpty()) {
//			System.out.println(".... new");
			// check if the only active edge is the root edge
			if (activeEdges.size() == 1 && activeEdges.get(0) == rootEdge && inactiveEdges.isEmpty()) {
				// simply add the segments to the root edge and finish
				activeEdges.get(0).hasSegments.or(segsToAddList.get(0));
				break;
			}
			
			int coalLines = 0;
			int negLines = 0;			
			
			List<NetworkEdge> excludedEdges = new ArrayList<>();
			List<NetworkEdge> excludedEdgesReverse = new ArrayList<>(); // excluded partners for calculation with inactive edges
			double[] reassortmentObsProb = new double[activeEdges.size()];
			double totalReassortmentProb = 0.0;
			double timeToNextEdgeEvent = Double.POSITIVE_INFINITY;
			int nextEventIndex = -1;
			int nextInactiveEventIndex = -1;
			for (int i = 0; i < activeEdges.size(); i++) {
				if (!activeEdges.get(i).hasSegments.isEmpty() || zeroActiveReassortment){
					int m = activeEdges.get(i).hasSegments.cardinality();
					int k = segsToAddList.get(i).cardinality();

//					// Check if edge is also being removed from (both active and inactive)
//					if (inactiveEdges.contains(activeEdges.get(i))) {
//						m -= segsToRemoveList.get(inactiveEdges.indexOf(activeEdges.get(i))).cardinality();
////						if (m==0) {
////							reassortmentObsProb[i] = 1.0 - 2.0 * Math.pow(0.5, k);
////						}else {
//							reassortmentObsProb[i] = 2.0 * (Math.pow(0.5, m)*(1-Math.pow(0.5, k)));
////						}
//					} else {
						reassortmentObsProb[i] = 2.0 * (Math.pow(0.5, m)*(1-Math.pow(0.5, k)));
//					}

					totalReassortmentProb += reassortmentObsProb[i];

					if (!activeEdges.get(i).isRootEdge() &&  timeToNextEdgeEvent > activeEdges.get(i).parentNode.getHeight()) {
						timeToNextEdgeEvent = activeEdges.get(i).parentNode.getHeight();
						nextEventIndex = i;
					}
				}else{
					// For newly created edges (not existing), they can coalesce
					coalLines++;
					// exlude newly created edges from reverse coalescence calculation
                    excludedEdgesReverse.add(activeEdges.get(i));
					// Calculate reassortment obs prob for the segments being added
					int nSegs = segsToAddList.get(i).cardinality();
					double obsProb = 1.0 - 2.0 * Math.pow(0.5, nSegs);
					reassortmentObsProb[i] = obsProb;
					totalReassortmentProb += reassortmentObsProb[i];
				}
			}

			double totalReverseReassortmentProb = 0.0;
			double[] reassortmentObsProbInactive = new double[inactiveEdges.size()];
			for (int i = 0; i < inactiveEdges.size(); i++) {
				if (!inactiveEdges.get(i).isRootEdge()
						&& timeToNextEdgeEvent > inactiveEdges.get(i).parentNode.getHeight()) {
					timeToNextEdgeEvent = inactiveEdges.get(i).parentNode.getHeight();
					nextInactiveEventIndex = i;
				}
					
				// check if the lineage would be empty after removing the segments
				if (inactiveEdges.get(i).hasSegments.cardinality() == segsToRemoveList.get(i).cardinality()) {
					int k = segsToRemoveList.get(i).cardinality();
					if (activeEdges.contains(inactiveEdges.get(i))) { 
						// the edge is simoultaneously being filled with new segments that can only undergo their seperate reassortment events.
						int m = segsToAddList.get(activeEdges.indexOf(inactiveEdges.get(i))).cardinality();
	                    reassortmentObsProbInactive[i] = 2.0 * (Math.pow(0.5, m)*(1-Math.pow(0.5, k)));
					}else { // truly inactive
	                    negLines++;
	                    excludedEdges.add(inactiveEdges.get(i));
						reassortmentObsProbInactive[i] = 1.0 - 2.0 * Math.pow(0.5, k);
					}
					totalReverseReassortmentProb += reassortmentObsProbInactive[i];
				}else {
					int m = inactiveEdges.get(i).hasSegments.cardinality();

					// Check if edge is also being added to (both active and inactive)
					if (activeEdges.contains(inactiveEdges.get(i))) {
						// the edge is simoultaneously being filled with other segments
						m += segsToAddList.get(activeEdges.indexOf(inactiveEdges.get(i))).cardinality();
					}

					int k = segsToRemoveList.get(i).cardinality();
					// reverse reassortment: Probability of observable reassortment for m segments
					// that will be on the edge after removing k segments, but none
					// of the m-k segments reassorted
					double obsProb = 2.0 * (Math.pow(0.5, m-k)*(1-Math.pow(0.5, k)));
					reassortmentObsProbInactive[i] = obsProb;
					totalReverseReassortmentProb += obsProb;
				}
			}
			// Sample time to next reassortment event
			double timeToNextReassortment = Double.POSITIVE_INFINITY;
			
			// clean up if there is supposed to be zero active or inactive reassortment, meaning that
			// until the first node is reached, we shouldn't allow for new branches
			if (totalReassortmentProb > 0 && !zeroActiveReassortment) {
	            if (coalescentDistr.timeVaryingReassortmentRates!= null) {
	                double currentTransformedReaTime = coalescentDistr.timeVaryingReassortmentRates.getIntensity(currentTime);
	                double transformedTimeToNextRea = Randomizer.nextExponential(totalReassortmentProb);
	            	timeToNextReassortment = coalescentDistr.timeVaryingReassortmentRates.getInverseIntensity(
	            			transformedTimeToNextRea + currentTransformedReaTime) - currentTime;
	            }else {
	            	timeToNextReassortment = Randomizer.nextExponential(totalReassortmentProb*coalescentDistr.reassortmentRateInput.get().getArrayValue());
	            }
			}	
			
			// Sample time to next coalescence using network intervals approach
			// Like AddRemoveReassortmentCoalescent, use the network event list to track lineages
			double timeToNextCoalescence = Double.POSITIVE_INFINITY;
			double coalRate = 0.0;
			double reverseCoalRate = 0.0;
			if (coalLines >= 1 || negLines >= 1) {
				
				double checkTime = currentTime;
				NetworkEvent prevEvent = null;
				
				// Calculate the earliest stopping time (reassortment or edge event)
				double stopTime = Math.min(currentTime + timeToNextReassortment, timeToNextEdgeEvent);
				
//				System.out.println("stop Time " + stopTime + " " + timeToNextReassortment + " " + timeToNextEdgeEvent + " " + stopTime + " " + coalLines + " " + negLines);
//				System.out.println(checkTime);
				// Find the appropriate interval and sample coalescence time
				for (NetworkEvent event : networkEventList) {
//					if (prevEvent != null)
//						System.out.println("exiss " + event.time + " " + prevEvent.time);
					if (event.time > checkTime) {
						// Check if we've reached the stop time (reassortment or edge event)
						if (checkTime >= stopTime) {
							break;
						}
//						System.out.println("ev " + event.time);
						
						// prevEvent.lineages is the number of existing lineages in the network
						// Add our new lineages (coalLines) to get total
//						int totalLineages = prevEvent.lineages + coalLines - negLines; 
						coalRate =  coalLines * (coalLines - 1) + coalLines * prevEvent.lineages;
						
//						int totalLineagesReverse = prevEvent.lineages;
						reverseCoalRate =  negLines * (negLines - 1) + negLines * (prevEvent.lineages - negLines);
																		
						// Sample coalescence time in this interval
						if (zeroActiveReassortment) {
							coalRate = 0.0;
						}
						if (zeroInactiveReassortment) {
							reverseCoalRate = 0.0;
						}
						double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(checkTime);
						double transformedTimeToNextCoal = Randomizer.nextExponential(coalRate);
						double coalescenceTime = coalescentDistr.populationFunction.getInverseIntensity(
								transformedTimeToNextCoal + currentTransformedTime);
												
						// calculate the hastings ratio contribution from here until min(event.time, stopTime)
						double minNextTime = Math.min(coalescenceTime, Math.min(event.time, stopTime));
						
						double integral = coalescentDistr.populationFunction.getIntegral(checkTime, minNextTime);
						double coalSurv = -coalRate * integral;
						double reverseCoalSurv = -reverseCoalRate * integral;
						logHR -= coalSurv;
						logHR += reverseCoalSurv;
						
//						System.out.println("coalTime " + coalescenceTime + " " + prevEvent.lineages + " " + coalLines + " " + checkTime);
						
						
						// Check if coalescence happens before next event AND before stop time
						if (coalescenceTime == minNextTime) {
							// add the event contribution
							double coalDens = Math.log(coalRate/coalescentDistr.populationFunction.getPopSize(coalescenceTime));
							logHR -= coalDens;
							timeToNextCoalescence = coalescenceTime - currentTime;
//							System.out.println("leaving " + currentTime);
							break;
						}
						
						// Move to next interval, but don't go past stop time
						checkTime = Math.min(event.time, stopTime);
					}
					prevEvent = event;
				}
//				System.out.println("coalTime " + timeToNextCoalescence);

//				System.out.println(stopTime + " " + checkTime + " " + coalLines + " " + negLines + " " + network.getRootEdge().childNode.getHeight());
				// If we've gone past all events and haven't reached stop time, coalescence can still happen above the root
				if (stopTime > checkTime && checkTime >= network.getRootEdge().childNode.getHeight() && coalLines > 0) {
//					int totalLineages = coalLines + 1;  // Only new lineages left (plus the existing root)
//					coalRate = 0.5 * coalLines * (totalLineages - 1);
					
//					System.out.println("root calculation");
					coalRate =  coalLines * (coalLines - 1) + coalLines * 1;
					if (zeroActiveReassortment)
						coalRate = 0.0;

					double startTime = Math.max(checkTime, network.getRootEdge().childNode.getHeight());
					
					double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(startTime);
					double transformedTimeToNextCoal = Randomizer.nextExponential(coalRate);
					double coalescenceTime = coalescentDistr.populationFunction.getInverseIntensity(
							transformedTimeToNextCoal + currentTransformedTime);
					
					
					double integral = coalescentDistr.populationFunction.getIntegral(startTime,coalescenceTime);
					
					double coalSurv2 = -coalRate * integral;
					logHR -= coalSurv2;

					// Reverse survival: negLines lineages among themselves
					if (negLines > 0) {
						throw new IllegalArgumentException("Negative lineages above root not allowed");
					}

					
					// Only use this if it happens before stop time
					if (coalescenceTime < stopTime) {
						timeToNextCoalescence = coalescenceTime - currentTime;
						double coalDens2 = Math.log(coalRate/coalescentDistr.populationFunction.getPopSize(coalescenceTime));
						logHR -= coalDens2;
					}
				}
			}
			
			
//			System.out.println(currentTime );
//			System.out.println(Arrays.toString(reassortmentObsProbInactive) + " " + Arrays.toString(reassortmentObsProb));
////			System.out.println(nextEventIndex + " " + nextInactiveEventIndex);
//			System.out.println(logHR + " " + currentTime + " " + timeToNextCoalescence + " " + timeToNextReassortment + " " + timeToNextEdgeEvent);
//			System.out.println(coalLines + " " + negLines + " " + coalRate + " " + reverseCoalRate + " " + totalReassortmentProb + " " + totalReverseReassortmentProb);
//			System.out.println(activeEdges.size() + " " + inactiveEdges.size());
//			
//			if (timeToNextCoalescence==Double.POSITIVE_INFINITY) {
//                throw new IllegalArgumentException("Coalescence time infinite " + coalLines + " " + currentTime);
//			}
			
			// Determine which event happens first
			double treeEventDiff = timeToNextEdgeEvent - currentTime;
			double timeUntilNextEvent = Math.min(Math.min(timeToNextCoalescence, timeToNextReassortment),
					treeEventDiff);
			
			boolean nextIsTreeEvent = timeUntilNextEvent == treeEventDiff;
			if (zeroActiveReassortment) {
				totalReassortmentProb = 0.0;
				if (nextIsTreeEvent && nextInactiveEventIndex == -1) {
					zeroActiveReassortment = false;
				}
			}
			if (zeroInactiveReassortment) {
				totalReverseReassortmentProb = 0.0;
				if (nextIsTreeEvent) {
					if (nextInactiveEventIndex != -1) { // loops are immediately ended, so do not need to be considered
						zeroInactiveReassortment = false;
					}
				}
			}
//			System.out.println(nextIsTreeEvent + " " + nextEventIndex + " " + nextInactiveEventIndex + " " + timeUntilNextEvent + " " + timeToNextCoalescence + " " + timeToNextReassortment + " " + (timeToNextEdgeEvent - currentTime));

			
			
//			System.out.println(zeroActiveReassortment + " " + zeroInactiveReassortment + " " + coalLines + " " + negLines);
//			System.out.println(totalReassortmentProb + " " + totalReverseReassortmentProb + " " + logHR + " " + currentTime + " " + (currentTime + timeUntilNextEvent));
			// Store the previous and new event times for logHR calculation
			double prevTime = currentTime;
			double eventTime = currentTime + timeUntilNextEvent;

			// Update current time
			currentTime = eventTime;
			

			// 2. Reassortment contribution (probability of no reassortment until eventTime, or reassortment at eventTime if that's the event)
			if (totalReassortmentProb > 0 || totalReverseReassortmentProb > 0) {
				double reassortSurvFwd = 0.0;
				double reassortSurvRev = 0.0;
				if (coalescentDistr.timeVaryingReassortmentRates != null) {
					double integral = coalescentDistr.timeVaryingReassortmentRates.getIntegral(prevTime, eventTime);
					reassortSurvFwd = -totalReassortmentProb * integral;
					reassortSurvRev = -totalReverseReassortmentProb * integral;
					logHR -= reassortSurvFwd;
					logHR += reassortSurvRev;
				} else {
					reassortSurvFwd = -totalReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue() * timeUntilNextEvent;
					reassortSurvRev = -totalReverseReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue() * timeUntilNextEvent;
					logHR -= reassortSurvFwd;
					logHR += reassortSurvRev;
				}
			}
//			System.out.println("logHR "+logHR + " " + currentTime);


			// Determine what type of event occurred and handle it
			if (timeUntilNextEvent == timeToNextCoalescence) {
				// add the event contribution
				double coalEventHR = handleCoalescenceEvent(coalLines, activeEdges, excludedEdges, segsToAddList, currentTime, rootEdge, edgesAdded);
				logHR += coalEventHR;
//				System.out.println(currentTime + " Coalescence event");
				if (currentTime==Double.POSITIVE_INFINITY)
					throw new IllegalArgumentException("Current time infinite");
			} else if (timeUntilNextEvent == timeToNextReassortment) {
				// add the event contribution
				double reassortDens = 0.0;
				if (coalescentDistr.timeVaryingReassortmentRates != null) {
					reassortDens = Math.log(totalReassortmentProb*coalescentDistr.timeVaryingReassortmentRates.getPopSize(currentTime));
				}else {
					reassortDens = Math.log(totalReassortmentProb*coalescentDistr.reassortmentRateInput.get().getArrayValue());
				}
				logHR -= reassortDens;
					
				double reassortEventHR = handleReassortmentEvent(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, reassortmentObsProb, totalReassortmentProb, currentTime, edgesAdded);
				logHR += reassortEventHR;
//				System.out.println(currentTime + " Reassortment event");
			} else {
				double edgeTraversalHR = handleEdgeTraversalEvent(nextInactiveEventIndex, nextEventIndex, inactiveEdges, segsToRemoveList, activeEdges, 
						segsToAddList, reverseCoalRate, totalReverseReassortmentProb, reassortmentObsProbInactive, negLines, excludedEdgesReverse);
				logHR += edgeTraversalHR;
//				System.out.println(currentTime + " Edge event ");
			}
		}
//		System.exit(0);
		//		networkEdges.addAll(edgesAdded);
		return logHR;
	}
	
	/**
	 * Handle a coalescence event among the active edges that are coal lines.
	 *
	 * @param coalLines     number of coal lines
	 * @param activeEdges   list of active edges
	 * @param excludedEdges list of excluded edges
	 * @param segsToAddList list of segments to add for each active edge
	 * @param currentTime   time of the coalescence event
	 * @param rootEdge      root edge of the network
	 * @param edgesAdded    list to track newly added edges
	 * @return log probability of the coalescence event
	 */
	private double handleCoalescenceEvent(int coalLines, List<NetworkEdge> activeEdges, List<NetworkEdge> excludedEdges, List<BitSet> segsToAddList,
			double currentTime, NetworkEdge rootEdge, List<NetworkEdge> edgesAdded) {
		double logHR = 0.0;
//		System.out.println("Handling coalescence event at time " + currentTime + " with " + coalLines + " coal lines.");
		
		// Step 1: sample which of the coalLines will coalesce
		int coalLineIdx1 = Randomizer.nextInt(coalLines);
		// Forward proposal probability: choosing which coal line
		logHR -= Math.log(1.0 / coalLines);
		
		// find the corresponding index within activeEdges
		int idx = -1;
		int count = 0;
		for (int i = 0; i < activeEdges.size(); i++) {
			if (activeEdges.get(i).hasSegments.isEmpty()) {
				if (count == coalLineIdx1) {
					idx = i;
					break;
				}
				count++;
			}
		}

		NetworkEdge edge1 = activeEdges.get(idx);


		// Step 2: sample which lineage it will coalesce with, can be other coalLines or any of the existing lineages,
		// but not any of the excludedEdges
		List<NetworkEdge> coexistingLineages = new ArrayList<>();
		for (NetworkEdge e : networkEdges) {
			if (e.isRootEdge()) {
				if (e.childNode.getHeight() < currentTime && !excludedEdges.contains(e)) {
					coexistingLineages.add(e);
				}
			}else {
			    if (e.parentNode.getHeight() > currentTime && 
			    		e.childNode.getHeight() < currentTime &&
			    		!excludedEdges.contains(e)) {
			        coexistingLineages.add(e);
			    }
			}
		}
		
		// remove the coalLineIdx1 from the coexistingLineages list
		coexistingLineages.remove(edge1);

//		// add the coalLines that are not coalLineIdx1 to the coexistingLineages list
//		count = 0;
//		for (int i = 0; i < activeEdges.size(); i++) {
//			if (activeEdges.get(i).hasSegments.isEmpty()) {
//				if (count != coalLineIdx1) {
//					coexistingLineages.add(activeEdges.get(i));
//				}
//				count++;
//			}
//		}
		if (coexistingLineages.contains(edge1)) {
			throw new IllegalArgumentException("Error: coexisting lineages contains edge1: this should not happen");
		}
		
//		for (NetworkEdge e : coexistingLineages) {
//			if (e.parentNode!=null)
//				System.out.println(e.childNode.getHeight() + " to " + e.parentNode.getHeight());
//			else
//				System.out.println(e.childNode.getHeight() + " to null");
//		}
//
//		for (NetworkEdge e : excludedEdges) {
//			if (e.parentNode!=null)
//				System.out.println(e.childNode.getHeight() + " ex " + e.parentNode.getHeight());
//			else
//				System.out.println(e.childNode.getHeight() + " ex null");
//			
//		}
//		System.out.println(network);
		

		int coalLineIdx2 = Randomizer.nextInt(coexistingLineages.size());
		// Forward proposal probability: choosing which partner lineage
		logHR -= Math.log(1.0 / coexistingLineages.size());

		// Step 3: create a new coalescent node and add it to the network
		NetworkNode coalescentNode = new NetworkNode();
		coalescentNode.setHeight(currentTime);

		NetworkEdge edge2 = coexistingLineages.get(coalLineIdx2);
		if (edge2.hasSegments.isEmpty()) { // there isn't an established edge and we need to add a new node
			// add the corresponding segs to add to both edges
			edge1.hasSegments.or(segsToAddList.get(idx));

			// find which index this edge corresponds to in the activeEdges list
			int idx2 = -1;
			for (int i = 0; i < activeEdges.size(); i++) {
				if (activeEdges.get(i) == edge2) {
					idx2 = i;
					break;
				}
			}
			edge2.hasSegments.or(segsToAddList.get(idx2));

			BitSet parentSegs = (BitSet) edge1.hasSegments.clone();
			parentSegs.or(edge2.hasSegments);

			coalescentNode.addChildEdge(edge2);
			coalescentNode.addChildEdge(edge1);
            NetworkEdge parentEdge = new NetworkEdge();
            parentEdge.hasSegments = new BitSet();
            coalescentNode.addParentEdge(parentEdge);

			activeEdges.remove(Math.max(idx, idx2));
			activeEdges.remove(Math.min(idx, idx2));
			segsToAddList.remove(Math.max(idx, idx2));
			segsToAddList.remove(Math.min(idx, idx2));

			activeEdges.add(parentEdge);
			segsToAddList.add(parentSegs);

			networkEdges.add(parentEdge);


		}else if (activeEdges.contains(edge2)) {
			int idx2 = -1;
			for (int i = 0; i < activeEdges.size(); i++) {
				if (activeEdges.get(i) == edge2) {
					idx2 = i;
					break;
				}
			}

			NetworkNode parent = new NetworkNode();
			parent.setHeight(currentTime);
			NetworkNode grandParent = edge2.parentNode;
			NetworkEdge parentEdge = new NetworkEdge();

			BitSet parentSegsToAdd = (BitSet) segsToAddList.get(idx).clone();
			parentSegsToAdd.andNot(edge2.hasSegments);
			parentSegsToAdd.or(segsToAddList.get(idx2));

			// Set parentEdge.hasSegments
			parentEdge.hasSegments = (BitSet) edge2.hasSegments.clone();

			edge1.hasSegments.or(segsToAddList.get(idx));
			edge2.hasSegments.or(segsToAddList.get(idx2));


			if (grandParent != null) {
				grandParent.removeChildEdge(edge2);
				grandParent.addChildEdge(parentEdge);
			}else if (!edge2.hasSegments.isEmpty()) {
				// we are making a new root edge
				// was root edge
				network.setRootEdge(parentEdge);
				rootEdge = parentEdge;
			}

			parent.addChildEdge(edge2);
			parent.addChildEdge(edge1);
			parent.addParentEdge(parentEdge);

			activeEdges.remove(Math.max(idx, idx2));
			activeEdges.remove(Math.min(idx, idx2));
			segsToAddList.remove(Math.max(idx, idx2));
			segsToAddList.remove(Math.min(idx, idx2));

			if (!parentSegsToAdd.isEmpty()){
				activeEdges.add(parentEdge);
				segsToAddList.add(parentSegsToAdd);
			}
			networkEdges.add(parentEdge);
		}else {
			NetworkNode grandParent = edge2.parentNode;
			if (grandParent != null) {
				grandParent.removeChildEdge(edge2);
			}
			coalescentNode.addChildEdge(edge2);
			coalescentNode.addChildEdge(edge1);

			edge1.hasSegments.or(segsToAddList.get(idx));


			NetworkEdge parentEdge = new NetworkEdge();
			parentEdge.hasSegments = (BitSet) edge2.hasSegments.clone();
			coalescentNode.addParentEdge(parentEdge);
			if (grandParent != null) {
				grandParent.addChildEdge(parentEdge);
			}else {
				// was root edge
				network.setRootEdge(parentEdge);
				rootEdge = parentEdge;
			}

			segsToAddList.get(idx).andNot(edge2.hasSegments);

			if (segsToAddList.get(idx).isEmpty()) {
				activeEdges.remove(idx);
				segsToAddList.remove(idx);
			}else {
				activeEdges.set(idx, parentEdge);
			}

			networkEdges.add(parentEdge);
		}

		return logHR;
	}

	private double handleReassortmentEvent(List<NetworkEdge> activeEdges, List<NetworkEdge> inactiveEdges,
			List<BitSet> segsToAddList, List<BitSet> segsToRemoveList,
			double[] reassortmentObsProb, double totalReassortmentProb, double currentTime, List<NetworkEdge> edgesAdded) {
		double logHR = 0.0;
//		System.out.println("Handling reassortment event at time " + currentTime + " with " + activeEdges.size() + " active edges.");
		// Step 1: sample which edge the reassortment occurs on (weighted by reassortmentObsProb)
		double randomChoice = Randomizer.nextDouble() * totalReassortmentProb;
		double cumulative = 0.0;
		int reassortEdgeIdx = -1;
		for (int i = 0; i < activeEdges.size(); i++) {
			cumulative += reassortmentObsProb[i];
			if (randomChoice <= cumulative) {
				reassortEdgeIdx = i;
				break;
			}
		}
		
		// Forward proposal probability: choosing which edge
		logHR -= Math.log(reassortmentObsProb[reassortEdgeIdx] / totalReassortmentProb);

		// Sample which segments will reassort (each with 0.5 probability, conditional on observable reassortment)
		BitSet segsToReassort = segsToAddList.get(reassortEdgeIdx);
		BitSet segsToStay = new BitSet();
		BitSet segsToGo = new BitSet();
		
		NetworkEdge parentToStay = new NetworkEdge();
		parentToStay.hasSegments = (BitSet) activeEdges.get(reassortEdgeIdx).hasSegments.clone();
		
		
		int segsBeingRemoved = 0;
		if (inactiveEdges.contains(activeEdges.get(reassortEdgeIdx))){
			int idx = inactiveEdges.indexOf(activeEdges.get(reassortEdgeIdx));
			activeEdges.get(reassortEdgeIdx).hasSegments.andNot(segsToRemoveList.get(idx));
			inactiveEdges.set(idx, parentToStay);
			segsBeingRemoved = segsToRemoveList.get(idx).cardinality();
		}

		
		// Sample segments with 0.5 probability each, conditional on at least one on each side
		if (activeEdges.get(reassortEdgeIdx).hasSegments.isEmpty() && segsBeingRemoved == 0) {
			do {
				segsToStay.clear();
				segsToGo.clear();
				for (int segIdx = segsToReassort.nextSetBit(0); segIdx != -1; segIdx = segsToReassort.nextSetBit(segIdx + 1)) {
					if (Randomizer.nextBoolean()) {
						segsToStay.set(segIdx);
					} else {
						segsToGo.set(segIdx);
					}
				}
			} while (segsToGo.cardinality() == 0 || segsToStay.cardinality() == 0);
			int k = segsToReassort.cardinality();

			// Forward proposal: k segments randomly assigned
			// Conditioning on at least one on each side is already in reassortmentObsProb (line 467)
			logHR -= k*Math.log(0.5);
//			logHR -= k*Math.log(Math.pow(2,  k-1)/(Math.pow(2,k)-2));


		} else {
			do {
				segsToStay.clear();
				segsToGo.clear();
				for (int segIdx = segsToReassort.nextSetBit(0); segIdx != -1; segIdx = segsToReassort.nextSetBit(segIdx + 1)) {
					if (Randomizer.nextBoolean()) {
						segsToStay.set(segIdx);
					} else {
						segsToGo.set(segIdx);
					}
				}
			} while (segsToGo.cardinality() == 0);
			int k = segsToReassort.cardinality();

			// Forward proposal: k segments randomly assigned, conditioned on at least one going
			// Base: -k*log(0.5), Conditioning: +log(2^k - 1)
			logHR -= k*Math.log(Math.pow(2,  k-1)/(Math.pow(2,k)-1));
			
		}

		// Step 2: create a new reassortment node and add it to the network
		NetworkNode newGrandParent = activeEdges.get(reassortEdgeIdx).parentNode;
		NetworkNode reassortmentNode = new NetworkNode();
		reassortmentNode.setHeight(currentTime);
		
		if (newGrandParent!=null) {
			newGrandParent.removeChildEdge(activeEdges.get(reassortEdgeIdx));
			newGrandParent.addChildEdge(parentToStay);
		}

		NetworkEdge parentToGo = new NetworkEdge();
		parentToGo.hasSegments = new BitSet();

		activeEdges.get(reassortEdgeIdx).hasSegments.or(segsToReassort);


		reassortmentNode.addParentEdge(parentToStay);
		reassortmentNode.addParentEdge(parentToGo);
		reassortmentNode.addChildEdge(activeEdges.get(reassortEdgeIdx));

		// replace the active edge with the parentToStay edge
		if (segsToStay.isEmpty()) {
			activeEdges.remove(reassortEdgeIdx);
			segsToAddList.remove(reassortEdgeIdx);
        } else {
			activeEdges.set(reassortEdgeIdx, parentToStay);
			segsToAddList.set(reassortEdgeIdx, segsToStay);
        }

		networkEdges.add(parentToGo);
		networkEdges.add(parentToStay);

		activeEdges.add(parentToGo);
		segsToAddList.add(segsToGo);

		return logHR;
	}

	private double handleEdgeTraversalEvent(int nextInactiveEventIndex, int nextEventIndex,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList,
			List<NetworkEdge> activeEdges, List<BitSet> segsToAddList, double reverseCoalRate, double totalReverseReassortmentProb, double[] reassortmentObsProbInactive, int reverseCoalLines, List<NetworkEdge> excludedEdgesReverse) {
		double logHR = 0.0;
		if (nextInactiveEventIndex != -1) {
			logHR += handleInactiveEdgeTraversal(nextInactiveEventIndex, inactiveEdges, segsToRemoveList, activeEdges, reverseCoalRate, totalReverseReassortmentProb, reassortmentObsProbInactive, reverseCoalLines, excludedEdgesReverse);
		} else { // the edge is both an active and an inactive edge
			logHR += handleActiveEdgeTraversal(nextEventIndex, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, reverseCoalRate, totalReverseReassortmentProb, reassortmentObsProbInactive, reverseCoalLines, excludedEdgesReverse);
		}
		return logHR;
	}

	private double handleInactiveEdgeTraversal(int idx, List<NetworkEdge> inactiveEdges, 
			List<BitSet> segsToRemoveList, List<NetworkEdge> activeEdges, double reverseCoalRate, double totalReverseReassortmentProb, double[] reassortmentObsProbInactive, int reverseCoalLines, List<NetworkEdge> excludedEdgesReverse) {
		double logHR = 0.0;
		if (inactiveEdges.get(idx).parentNode.isReassortment()) {
			logHR += handleReassortment(activeEdges, segsToRemoveList, inactiveEdges, segsToRemoveList, inactiveEdges.get(idx), totalReverseReassortmentProb, reassortmentObsProbInactive);
		} else {
			NetworkEdge sisterEdge = getSisterEdge(inactiveEdges.get(idx));
			logHR += handleCoalescence(activeEdges, segsToRemoveList, inactiveEdges, segsToRemoveList, inactiveEdges.get(idx), sisterEdge, reverseCoalRate, reverseCoalLines, excludedEdgesReverse);
		}
		return logHR;
	}

	private double handleActiveEdgeTraversal(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, double reverseCoalRate, double totalReverseReassortmentProb, double[] reassortmentObsProbInactive, int reverseCoalLines, List<NetworkEdge> excludedEdgesReverse) {
		double logHR = 0.0;
		if (activeEdges.get(idx).parentNode.isReassortment()) {
			logHR += handleReassortment(activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, activeEdges.get(idx), totalReverseReassortmentProb, reassortmentObsProbInactive);
		} else {
			NetworkEdge sisterEdge = getSisterEdge(activeEdges.get(idx));
			logHR += handleCoalescence(activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, activeEdges.get(idx), sisterEdge, reverseCoalRate, reverseCoalLines, excludedEdgesReverse);
		}
		return logHR;
	}

	private double handleReassortment(List<NetworkEdge> activeEdges, List<BitSet> segsToAddList,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, NetworkEdge edge, double totalReverseReassortmentProb, double[] reassortmentObsProbInactive) {

		double logHR = 0.0;
		
		int idx = -1;
		int iidx = -1;

		for (int i = 0; i < activeEdges.size(); i++) {
			if (activeEdges.get(i) == edge) {
				idx = i;
				break;
			}
		}

		for (int i = 0; i < inactiveEdges.size(); i++) {
			if (inactiveEdges.get(i) == edge) {
				iidx = i;
				break;
			}
		}
		
		BitSet segsToGoLeft = new BitSet();
		BitSet segsToGoRight = new BitSet();
		// randomly assign each segment to go left or right
		if (idx != -1) {
			for (int segIdx = segsToAddList.get(idx).nextSetBit(0); segIdx != -1; segIdx = segsToAddList.get(idx).nextSetBit(segIdx + 1)) {
				if (Randomizer.nextBoolean()) {
					segsToGoLeft.set(segIdx);
				} else {
					segsToGoRight.set(segIdx);
				}
			}
			// Forward proposal probability: For each segment, choose a parent at random
			logHR -= Math.log(0.5) * segsToAddList.get(idx).cardinality();
			
			edge.hasSegments.or(segsToAddList.get(idx));
			activeEdges.remove(idx);
			segsToAddList.remove(idx);

		}
		
		// now, check what to remove
		BitSet segsToRemoveLeft = new BitSet();
		BitSet segsToRemoveRight = new BitSet();
		

		// remove from inactive edges
		if (iidx != -1) {
			int segsCard = segsToRemoveList.get(iidx).cardinality();
			segsToRemoveLeft.or(segsToRemoveList.get(iidx));
			segsToRemoveRight.or(segsToRemoveList.get(iidx));
			
			segsToRemoveLeft.and(edge.parentNode.getParentEdges().get(0).hasSegments);
			segsToRemoveRight.and(edge.parentNode.getParentEdges().get(1).hasSegments);
			// Reverse proposal probability: For each segment, choose a parent at random

			boolean createdEmptyEdge = false;
			if (inactiveEdges.get(iidx).hasSegments.cardinality() == segsCard){// move creates two empty parents
				createdEmptyEdge = true;
				// Reverse proposal: k segments randomly assigned, conditioned on at least one on each side
				// Base: +k*log(0.5), Conditioning: -log(2^k - 2)
				logHR += segsCard * Math.log(0.5);

			}else if (segsToRemoveLeft.cardinality()==edge.parentNode.getParentEdges().get(0).hasSegments.cardinality() && segsToGoLeft.isEmpty()){ // move creates one empty parent
				createdEmptyEdge = true;
				// contional draws, conditioning on at least least one segment having originated from the other parent
				logHR += segsCard*Math.log(Math.pow(2,  segsCard-1)/(Math.pow(2,segsCard)-1));

			}else if (segsToRemoveRight.cardinality()==edge.parentNode.getParentEdges().get(1).hasSegments.cardinality() && segsToGoRight.isEmpty()){  // move creates one empty parent
				createdEmptyEdge = true;

				// Reverse proposal: k segments randomly assigned, conditioned on at least one going
				// Base: +k*log(0.5), Conditioning: -log(2^k - 1)
				logHR += segsCard*Math.log(Math.pow(2,  segsCard-1)/(Math.pow(2,segsCard)-1));

			}
//			System.out.println("Reassortment traversal before: " + logHR + " " + createdEmptyEdge);
//			System.out.println(totalReverseReassortmentProb);
			// the totalReverseReassortmentProb condition is to ignore the first event
			if (totalReverseReassortmentProb>0.0) {
				if (createdEmptyEdge) {
					// Reverse proposal probability: choosing which edge
					logHR += Math.log(reassortmentObsProbInactive[iidx] / totalReverseReassortmentProb);
					
					if (coalescentDistr.timeVaryingReassortmentRates != null) {
						double popSize = coalescentDistr.timeVaryingReassortmentRates.getPopSize(edge.parentNode.getHeight());
						logHR += Math.log(totalReverseReassortmentProb * popSize);
					} else {
						logHR += Math.log(totalReverseReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue());
					}			
				}else {
					logHR += Math.log(0.5) * segsCard;
				}
			}
//			System.out.println("rafter: " + logHR);

			
			edge.hasSegments.andNot(segsToRemoveList.get(iidx));
			inactiveEdges.remove(iidx);
			segsToRemoveList.remove(iidx);



		}

		// now, check if we add
		if (!segsToGoLeft.isEmpty()) {
			activeEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToAddList.add(segsToGoLeft);
		}
		if (!segsToGoRight.isEmpty()) {
			activeEdges.add(edge.parentNode.getParentEdges().get(1));
			segsToAddList.add(segsToGoRight);
		}
		
		if (!segsToRemoveLeft.isEmpty()) {
			inactiveEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToRemoveList.add(segsToRemoveLeft);
		}
		if (!segsToRemoveRight.isEmpty()) {
			inactiveEdges.add(edge.parentNode.getParentEdges().get(1));
			segsToRemoveList.add(segsToRemoveRight);
		}
		
//		System.out.println("After reassortment traversal: " + logHR);

		return logHR;
	}
	
	private double handleCoalescence(List<NetworkEdge> activeEdges, List<BitSet> segsToAddList,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, NetworkEdge edge,  NetworkEdge sisterEdge, double reverseCoalRate, int reverseCoalLines, List<NetworkEdge> excludedEdgesReverse) {
		
		double logHR = 0.0;
		
		int idx1 = -1;
		int idx2 = -1;
		
		int iidx1 = -1;
		int iidx2 = -1;
		
		BitSet segsToAdd1 = new BitSet();
		BitSet segsToAdd2 = new BitSet();
		
		BitSet segsToRemove1 = new BitSet();
		BitSet segsToRemove2 = new BitSet();
		
		BitSet segsBefore1 = (BitSet) edge.hasSegments.clone();
		BitSet segsBefore2 = (BitSet) sisterEdge.hasSegments.clone();
		
		for (int i = 0; i < activeEdges.size(); i++) {
			if (activeEdges.get(i) == edge) {
				idx1 = i;
				segsToAdd1.or(segsToAddList.get(i));
			} else if (activeEdges.get(i) == sisterEdge) {
				idx2 = i;
				segsToAdd2.or(segsToAddList.get(i));
			}
		}
		
		for (int i = 0; i < inactiveEdges.size(); i++) {
			if (inactiveEdges.get(i) == edge) {
				iidx1 = i;
				segsToRemove1.or(segsToRemoveList.get(i));
			} else if (inactiveEdges.get(i) == sisterEdge) {
				iidx2 = i;
				segsToRemove2.or(segsToRemoveList.get(i));
			}
		}
		
//		System.out.println(reverseCoalRate);

		if (reverseCoalRate>0.0) {
			boolean calculateReverseProb = false;
			if (iidx1 != -1 && segsToRemoveList.get(iidx1).cardinality()==inactiveEdges.get(iidx1).hasSegments.cardinality() ) {
				calculateReverseProb = true;
			}
			if (iidx2 != -1 && segsToRemoveList.get(iidx2).cardinality() == inactiveEdges.get(iidx2).hasSegments.cardinality()) {
				calculateReverseProb = true;
			}
			if (calculateReverseProb) {

				logHR += Math.log(reverseCoalRate/(coalescentDistr.populationFunction.getPopSize(edge.parentNode.getHeight())));
//				System.out.println("a " + logHR);
				
				// Count coexisting lineages at the time of coalescence (edge.parentNode.getHeight())
				List<NetworkEdge> reverseCoexistingLineages = new ArrayList<>();
				for (NetworkEdge e : networkEdges) {
					if (e.isRootEdge()) {
						if (e.childNode.getHeight() < edge.parentNode.getHeight() && !excludedEdgesReverse.contains(e)) {
							reverseCoexistingLineages.add(e);
						}
					} else {
						if (e.parentNode.getHeight() > edge.parentNode.getHeight() && 
						    e.childNode.getHeight() < edge.parentNode.getHeight() && !excludedEdgesReverse.contains(e)) {
							reverseCoexistingLineages.add(e);
						}
					}
				}
				
				// Check if sister edge is also becoming empty
				boolean sisterAlsoEmpty = false;
				if (iidx1 != -1 && iidx2 != -1) {
				    // Both edges becoming empty - sister is already in negLines
				    sisterAlsoEmpty = true;
				}

				if (reverseCoalLines > 0) {
//					System.out.println(reverseCoalLines);
				    logHR += Math.log(1.0 / reverseCoalLines);
//				    int partnerCount = sisterAlsoEmpty ? reverseCoexistingLineages.size() : reverseCoexistingLineages.size() + 1;

				    int partnerCount = reverseCoexistingLineages.size() + 1;
//				    partnerCount = Math.max(1, partnerCount);
                    
				    logHR += Math.log(1.0 / partnerCount);
//				    System.out.println("c " + logHR);
				}else {
					throw new RuntimeException("Error calculating reverse coalescence probability, reverseCoalLines: " + reverseCoalLines + ", reverseCoexistingLineages: " + reverseCoexistingLineages.size());
				}
			}
		}
		
		
		edge.hasSegments.or(segsToAdd1);
		sisterEdge.hasSegments.or(segsToAdd2);
		
		edge.hasSegments.andNot(segsToRemove1);
		sisterEdge.hasSegments.andNot(segsToRemove2);
		
		// now, remove
		segsToAdd1.andNot(segsBefore2);
		segsToAdd2.andNot(segsBefore1);
		
		segsToRemove1.andNot(sisterEdge.hasSegments);
		segsToRemove2.andNot(edge.hasSegments);
		
		
		// remove from inactive edges
		if (iidx1 != -1 &&  iidx2 != -1) {
			inactiveEdges.remove(Math.max(iidx1, iidx2));
			segsToRemoveList.remove(Math.max(iidx1, iidx2));
			inactiveEdges.remove(Math.min(iidx1, iidx2));
			segsToRemoveList.remove(Math.min(iidx1, iidx2));

		} else if (iidx1 != -1 || iidx2 != -1) {
			int removeIdx = iidx1 != -1 ? iidx1 : iidx2;
			inactiveEdges.remove(removeIdx);
			segsToRemoveList.remove(removeIdx);
			// check the number of potential partners
		}
		

		
		// remove from active edges
		if (idx1 != -1 && idx2 != -1) {
			activeEdges.remove(Math.max(idx1, idx2));
			segsToAddList.remove(Math.max(idx1, idx2));
			activeEdges.remove(Math.min(idx1, idx2));
			segsToAddList.remove(Math.min(idx1, idx2));
		} else if (idx1 != -1 || idx2 != -1) {
			int removeIdx = idx1 != -1 ? idx1 : idx2;
			activeEdges.remove(removeIdx);
			segsToAddList.remove(removeIdx);
		}
		
		// now, check if we add
		segsToAdd1.or(segsToAdd2);
		segsToRemove1.or(segsToRemove2);

		if (!segsToAdd1.isEmpty()) {
			activeEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToAddList.add(segsToAdd1);
		}
		if (!segsToRemove1.isEmpty()) {
			inactiveEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToRemoveList.add(segsToRemove1);
		}
//		System.out.println("After coale traversal: " + logHR);

		return logHR;
	}


	/**
	 * Clean up empty reassortment edges.
	 * Uses the same pattern as EmptyEdgesNetworkOperator.RemoveAllEmptyNetworkSegments()
	 */
	protected void cleanEmptyEdgesTopDown() {

		// Create a local list of network edges to track during removal
//		System.out.println(network);
		List<Integer> removableEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge() && edge.childNode.isReassortment()
					&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() == 0) {
				removableEdges.add(i);
			}
		}

		// Remove all empty edges iteratively
		while (removableEdges.size() > 0) {
			int edgeInd = 0; // deterministically remove edges as their random contribution is already accounted for

			removeEmptyReassortmentEdgeAdapted(networkEdges, networkEdges.get(removableEdges.get(edgeInd)));

			// Rebuild the list of removable edges
			removableEdges = new ArrayList<>();
			for (int i = 0; i < networkEdges.size(); i++) {
				NetworkEdge edge = networkEdges.get(i);
				if (!edge.isRootEdge() && edge.childNode.isReassortment()
						&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() == 0) {
					removableEdges.add(i);
				}
			}
		}

	}

	/**
	 * Remove an empty reassortment edge, following the exact pattern from EmptyEdgesNetworkOperator.
	 * Updates the localNetworkEdges list to track removed edges.
	 */
	private void removeEmptyReassortmentEdgeAdapted(List<NetworkEdge> localNetworkEdges, NetworkEdge edgeToRemove) {

		NetworkNode nodeToRemove = edgeToRemove.childNode;
		NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
		NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

		// Update local edges list (same as parent class lines 355-356)
		localNetworkEdges.remove(edgeToRemove);
		localNetworkEdges.remove(edgeToRemoveSpouse);

		// Remove edge and associated nodes
		NetworkEdge edgeToExtend = nodeToRemove.getChildEdges().get(0);
		nodeToRemove.removeChildEdge(edgeToExtend);
		nodeToRemove.removeParentEdge(edgeToRemove);
		nodeToRemove.removeParentEdge(edgeToRemoveSpouse);
		edgeToRemoveSpouseParent.removeChildEdge(edgeToRemoveSpouse);
		edgeToRemoveSpouseParent.addChildEdge(edgeToExtend);

		NetworkNode secondNodeToRemove = edgeToRemove.parentNode;
		NetworkEdge secondEdgeToExtend = getSisterEdge(edgeToRemove);

		secondNodeToRemove.removeChildEdge(secondEdgeToExtend);
		secondNodeToRemove.removeChildEdge(edgeToRemove);

		// Update local edges list (same as parent class line 372)
		localNetworkEdges.remove(secondNodeToRemove.getParentEdges().get(0));
		
		networkEdges.remove(secondNodeToRemove.getParentEdges().get(0));
		networkEdges.remove(edgeToRemove);
		networkEdges.remove(edgeToRemoveSpouse);

		if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
			network.setRootEdge(secondEdgeToExtend);
		} else {
			NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
			NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
			secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
			secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);
			secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
		}
	}

	protected BitSet getRandomUnconditionedSubset(BitSet sourceSegments) {
		BitSet destSegments = new BitSet();

		destSegments.clear();

		for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1; segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

			if (Randomizer.nextBoolean())
				destSegments.set(segIdx);
		}

		return destSegments;
	}

	protected double getLogUnconditionedSubsetProb(BitSet sourceSegments) {
		return sourceSegments.cardinality() * Math.log(0.5);
	}
	
    protected double getLogConditionedSubsetProb(BitSet sourceSegments, BitSet chosenSubset, double binomProb) {

        if (sourceSegments.cardinality() < 2) {
            return Double.NEGATIVE_INFINITY;
        }

        int n = sourceSegments.cardinality();
        int k = chosenSubset.cardinality();

        double logProb = k * Math.log(binomProb)
                       + (n - k) * Math.log(1.0 - binomProb)
                       - Math.log(1.0 - Math.pow(binomProb, n) - Math.pow(1.0 - binomProb, n));

        return logProb;
    }
    
    protected BitSet getRandomConditionedSubset(BitSet sourceSegments, double binomProb) {

        if (sourceSegments.cardinality() < 2) {
            return null;
        }

        BitSet destSegments = new BitSet();

        do {
            destSegments.clear();

            for (int segIdx = sourceSegments.nextSetBit(0);
                 segIdx != -1;
                 segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

                if (Randomizer.nextDouble() < binomProb) {
                    destSegments.set(segIdx);
                }
            }

        } while (destSegments.cardinality() == 0
                || destSegments.cardinality() == sourceSegments.cardinality());
        return destSegments;
    }

	protected BitSet getRandomSegment(BitSet hasSegments) {
		// pick a random segment to divert.
		BitSet segsToDivert = new BitSet();
		segsToDivert.clear();

		int segToDivert = Randomizer.nextInt(hasSegments.cardinality());
		int segIdx = hasSegments.nextSetBit(0);
		for (int i = 0; i < segToDivert; i++) {
			segIdx = hasSegments.nextSetBit(segIdx + 1);
		}
		segsToDivert.set(segIdx);
		return segsToDivert;
	}

}

