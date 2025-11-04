package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class DivertSegmentOperatorWithReassortment extends EmptyEdgesNetworkOperator {

	public Input<Boolean> divertOneSegmentInput = new Input<>("divertOneSegment",
			"If true, only one segment is diverted", false);

	public Input<CoalescentWithReassortment> coalescentDistrInput = new Input<>("coalescentWithReassortment",
			"Coalescent distribution for sampling reassortment events.",
			Input.Validate.OPTIONAL);

	private CoalescentWithReassortment coalescentDistr;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		if (coalescentDistrInput.get() != null) {
			coalescentDistr = coalescentDistrInput.get();
		}
	}
	
	boolean exitOnError = false;

	@Override
	public double networkProposal() {
		double logHR = 0.0;
//		System.out.println(network.getExtendedNewickVerbose(0));
		
		List<Integer> sourceEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			if (networkEdges.get(i).childNode.isReassortment() 
					&& networkEdges.get(i).hasSegments.cardinality() > 0) {
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
			// throw not implemetned error
			throw new IllegalArgumentException("Divert one segment not implemented yet");
		} else {
			segsToDivert = getRandomUnconditionedSubset(sourceEdge.hasSegments);
			logHR -= getLogUnconditionedSubsetProb(sourceEdge.hasSegments);
		}

		if (segsToDivert.cardinality() == 0)
			return Double.NEGATIVE_INFINITY;

		logHR += divertSegments(destEdge, sourceEdge, segsToDivert);


		if (divertOneSegmentInput.get()) {
			logHR += Math.log(1.0 / destEdge.hasSegments.cardinality());
		} else {
			logHR += getLogUnconditionedSubsetProb(destEdge.hasSegments);
		}
		
		
		int reverseSourceEdgeCount = 0;
		for (NetworkEdge edge : network.getEdges()){
			if (edge.childNode.isReassortment()){
                reverseSourceEdgeCount++;
            }
		}
//		for (int i = 0; i < .size(); i++) {
//			if (networkEdges.get(i).childNode.isReassortment() && networkEdges.get(i).hasSegments.cardinality() > 0) {
//				reverseSourceEdgeCount++;
//			}
//		}

		logHR += Math.log(1.0 / reverseSourceEdgeCount);
//		System.out.println(network.getExtendedNewickVerbose(0) + "\n");

		return logHR;
	}

	protected double divertSegments(NetworkEdge destEdge, NetworkEdge sourceEdge, BitSet segsToDivert) {
		double logHR = 0.0;

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];

		getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);

//		System.out.println("source Edge c p heights " + sourceEdge.childNode.getHeight() + " " + sourceEdge.parentNode.getHeight());
//		System.out.println("dest   Edge c p heights " + destEdge.childNode.getHeight() + " " + destEdge.parentNode.getHeight());
//		System.out.println(segsToDivert);
		exitOnError=false;

		List<NetworkEdge> edgesAdded = new ArrayList<>();
//		System.out.println(segsToDivert + " " + destEdge.childNode.getHeight() + " " + destEdge.parentNode.getHeight());
//		System.out.println(network.getExtendedNewickVerbose(0));
		
		List<Double> eventTimes = new ArrayList<>();
		
		logHR -= divertSegmentsToAncestors(sourceEdge, destEdge, segsToDivert, edgesAdded, eventTimes);
//		System.out.println(network.getExtendedNewickVerbose(0));
//		logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert);
//		System.out.println(network.getExtendedNewickVerbose(0));

		if (reconnectSegmentTrees(treeChildNodeList, destEdge, segsToDivert))
			return Double.NEGATIVE_INFINITY;

		logHR += cleanEmptyEdgesTopDown();
		
		
		
//		System.out.println(network.getExtendedNewickVerbose(0));
//
//		if (exitOnError) {
//			System.exit(0);
//		}
		
		// Validation: Check that all true coalescent nodes have segmentIndices for their segments
		for (NetworkEdge edge : network.getEdges()) {
			if (!edge.isRootEdge() && edge.getLength()<0) {
				System.out.println(eventTimes);
				System.out.println(network.getExtendedNewickVerbose(0));
				System.exit(0);
			}
			
			if (edge.parentNode != null && edge.parentNode.isCoalescence()) {
				NetworkEdge sisterEdge = getSisterEdge(edge);
				// This is a true coalescent node if both child edges have segments
				for (int segIdx = 0; segIdx < network.getSegmentCount(); segIdx++) {
					if (edge.hasSegments.get(segIdx) && sisterEdge.hasSegments.get(segIdx)) {
						// Both children have this segment, so parent must be a coalescent node for it
						if (edge.parentNode.segmentIndices == null) {
							System.err.println("ERROR: Coalescent node at height " + edge.parentNode.getHeight() +
									" has no segmentIndices for segment " + segIdx);
							System.err.println("Edge: " + edge.childNode.getHeight() + " -> " + edge.parentNode.getHeight());
							System.err.println("Sister: " + sisterEdge.childNode.getHeight() + " -> " + sisterEdge.parentNode.getHeight());
							System.exit(0);
						}else {
							// check that the segment index is set correctly
							
						}
					}
				}
			}
		}

//		if (edgesAdded.size()== 3) {
//			System.out.println(network.getExtendedNewickVerbose(0));
//			System.exit(0);
//		}
//			
		
//		System.out.println(destEdge.hasSegments + " " + sourceEdge.hasSegments);
		
		

		return logHR;
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

	/**
	 * Add segments to this edge and ancestors, potentially sampling new reassortment events.
	 *
	 * @param edge      edge at which to start addition
	 * @param segsToAdd segments to add to the edge and ancestors
	 * @param eventTimes 
	 * @return log probability of operation
	 */
	protected double divertSegmentsToAncestors(NetworkEdge sourceEdge, NetworkEdge destEdge, BitSet segsToAdd, List<NetworkEdge> edgesAdded, List<Double> eventTimes) {
		double logHR = 0.0;
		
		segsToAdd = (BitSet) segsToAdd.clone();
		// startwith keeping a list of active edges and the segs to Add
		List<NetworkEdge> activeEdges = new ArrayList<>();
		List<BitSet> segsToAddList = new ArrayList<>();
		
		List<NetworkEdge> inactiveEdges = new ArrayList<>();
		List<BitSet> segsToRemoveList = new ArrayList<>();
		
		activeEdges.add(destEdge);
		segsToAddList.add(segsToAdd);
		
		inactiveEdges.add(sourceEdge);
		segsToRemoveList.add((BitSet) segsToAdd.clone());

		// Track the current time (start at the bottom of the lowest edge)
		double currentTime = destEdge.childNode.getHeight();
		NetworkEdge rootEdge = network.getRootEdge();
		
		// Get network events to determine lineages at different times, they will contain any inactive Edges, but not newly added active ones
		List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();



		// now, sample the time to the next reassortment or coalescent event on the active edges,
		// i.e. essentially simulate the history of these edges and segments until there is nothing left. The
		// goal is to avoid the creation of empty
		while (!activeEdges.isEmpty() || !inactiveEdges.isEmpty()) {
			// check if the only active edge is the root edge
			if (activeEdges.size() == 1 && activeEdges.get(0) == rootEdge && inactiveEdges.isEmpty()) {
				// simply add the segments to the root edge and finish
				activeEdges.get(0).hasSegments.or(segsToAddList.get(0));
				break;
			}
			
			int coalLines = 0;
			int negLines = 0;
			List<NetworkEdge> excludedEdges = new ArrayList<>();
			double[] reassortmentObsProb = new double[activeEdges.size()];
			double totalReassortmentProb = 0.0;
			double timeToNextEdgeEvent = Double.POSITIVE_INFINITY;
			int nextEventIndex = -1;
			int nextInactiveEventIndex = -1;
			for (int i = 0; i < activeEdges.size(); i++) {
				if (!activeEdges.get(i).hasSegments.isEmpty()){
					int m = activeEdges.get(i).hasSegments.cardinality();
					int k = segsToAddList.get(i).cardinality();

					double obsProb = Math.pow(2, 1 - m) - Math.pow(2, 1 - (m + k));

					reassortmentObsProb[i] = obsProb;
					totalReassortmentProb += reassortmentObsProb[i];
					
					if (!activeEdges.get(i).isRootEdge() &&  timeToNextEdgeEvent > activeEdges.get(i).parentNode.getHeight()) {
						timeToNextEdgeEvent = activeEdges.get(i).parentNode.getHeight();
						nextEventIndex = i;
					}
				}else{
					// For newly created edges (not existing), they can coalesce
					coalLines++;
					// Calculate reassortment obs prob for the segments being added
					int nSegs = segsToAddList.get(i).cardinality();
					double obsProb = 1.0 - 2.0 * Math.pow(0.5, nSegs);
					reassortmentObsProb[i] = nSegs * obsProb;
					totalReassortmentProb += reassortmentObsProb[i];
				}
			}

			for (int i = 0; i < inactiveEdges.size(); i++) {
				if (!inactiveEdges.get(i).isRootEdge()
						&& timeToNextEdgeEvent > inactiveEdges.get(i).parentNode.getHeight()) {
					timeToNextEdgeEvent = inactiveEdges.get(i).parentNode.getHeight();
					nextInactiveEventIndex = i;
					
					// check if the lineage would be empty after removing the segments
					if (inactiveEdges.get(i).hasSegments.cardinality() == segsToRemoveList.get(i).cardinality() && !activeEdges.contains(inactiveEdges.get(i))) {
                        negLines++;
                        excludedEdges.add(inactiveEdges.get(i));
					}
				}
			}
			
			// Sample time to next reassortment
			double timeToNextReassortment = Double.POSITIVE_INFINITY;
			if (totalReassortmentProb > 0) {

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
			if (coalLines >= 1) {
				
				double checkTime = currentTime;
				NetworkEvent prevEvent = null;
				double sampledCoalTime = 0.0;
				
				// Calculate the earliest stopping time (reassortment or edge event)
				double stopTime = Math.min(currentTime + timeToNextReassortment, timeToNextEdgeEvent);
				
				// Find the appropriate interval and sample coalescence time
				for (NetworkEvent event : networkEventList) {
					if (event.time > checkTime) {
						// Check if we've reached the stop time (reassortment or edge event)
						if (checkTime >= stopTime) {
							break;
						}
						
						// prevEvent.lineages is the number of existing lineages in the network
						// Add our new lineages (coalLines) to get total
						int totalLineages = prevEvent.lineages + coalLines - negLines; 
						coalRate = 0.5 * coalLines * (totalLineages - 1); // TODO actually needs to represent how many pairs of combinations with the coalLines exist
						
						// Sample coalescence time in this interval
						double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(checkTime);
						double transformedTimeToNextCoal = Randomizer.nextExponential(coalRate);
						double coalescenceTime = coalescentDistr.populationFunction.getInverseIntensity(
								transformedTimeToNextCoal + currentTransformedTime);
						
						sampledCoalTime = coalescenceTime;
						
						// Check if coalescence happens before next event AND before stop time
						if (coalescenceTime < Math.min(event.time, stopTime)) {
							timeToNextCoalescence = coalescenceTime - currentTime;
							break;
						}
						
						// Move to next interval, but don't go past stop time
						checkTime = Math.min(event.time, stopTime);
					}
					prevEvent = event;
				}
				
				// If we've gone past all events and haven't reached stop time, coalescence can still happen above the root
				if ((sampledCoalTime > network.getRootEdge().childNode.getHeight() || timeToNextCoalescence == Double.POSITIVE_INFINITY) 
						&& checkTime < stopTime) {
					int totalLineages = coalLines+1;  // Only new lineages left
					coalRate = 0.5 * coalLines * (totalLineages - 1);
					
					double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(
							Math.max(checkTime, network.getRootEdge().childNode.getHeight()));
					double transformedTimeToNextCoal = Randomizer.nextExponential(coalRate);
					double coalescenceTime = coalescentDistr.populationFunction.getInverseIntensity(
							transformedTimeToNextCoal + currentTransformedTime);
					
					// Only use this if it happens before stop time
					if (coalescenceTime < stopTime) {
						timeToNextCoalescence = coalescenceTime - currentTime;
					}
				}
			}
			
			// Determine which event happens first
			double timeUntilNextEvent = Math.min(Math.min(timeToNextCoalescence, timeToNextReassortment),
					timeToNextEdgeEvent - currentTime);

			// Store the previous and new event times for logHR calculation
			double prevTime = currentTime;
			double eventTime = currentTime + timeUntilNextEvent;

			// Update current time
			currentTime = eventTime;
			eventTimes.add(currentTime);
						
			// Add logHR contributions for all competing processes up to the event time

			// 1. Coalescence contribution (probability of no coalescence until eventTime, or coalescence at eventTime if that's the event)
			if (coalRate > 0) {
				logHR -= -coalRate * coalescentDistr.populationFunction.getIntegral(prevTime, eventTime);
				if (timeUntilNextEvent == timeToNextCoalescence) {
					// If coalescence occurred, add the density at this time
					logHR -= Math.log(coalRate / coalescentDistr.populationFunction.getPopSize(eventTime));
				}
			}

			// 2. Reassortment contribution (probability of no reassortment until eventTime, or reassortment at eventTime if that's the event)
			if (totalReassortmentProb > 0) {
				if (coalescentDistr.timeVaryingReassortmentRates != null) {
					double prevTransformedTime = coalescentDistr.timeVaryingReassortmentRates.getIntensity(prevTime);
					double currTransformedTime = coalescentDistr.timeVaryingReassortmentRates.getIntensity(eventTime);
					logHR -= -totalReassortmentProb * (currTransformedTime - prevTransformedTime);
					if (timeUntilNextEvent == timeToNextReassortment) {
						// If reassortment occurred, add the density at this time
						logHR -= Math.log(totalReassortmentProb * coalescentDistr.timeVaryingReassortmentRates.getPopSize(eventTime));
					}
				} else {
					logHR -= -totalReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue() * timeUntilNextEvent;
					if (timeUntilNextEvent == timeToNextReassortment) {
						// If reassortment occurred, add the density at this time
						logHR -= Math.log(totalReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue());
					}
				}
			}

			// Determine what type of event occurred and handle it
			if (timeUntilNextEvent == timeToNextCoalescence) {
//				System.out.println("Coalescence event at time " + currentTime + " with " + coalLines + " coal lines");
				logHR += handleCoalescenceEvent(coalLines, activeEdges, excludedEdges, segsToAddList, currentTime, rootEdge, edgesAdded);
//				return Double.NEGATIVE_INFINITY;
			} else if (timeUntilNextEvent == timeToNextReassortment) {
//				System.out.println("Reassortment event at time " + currentTime);
				logHR += handleReassortmentEvent(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, reassortmentObsProb, totalReassortmentProb, currentTime, edgesAdded);
//				return Double.NEGATIVE_INFINITY;
			} else {
//				System.out.println("Edge traversal event at time " + currentTime + " " + nextInactiveEventIndex + " " + nextEventIndex);
				handleEdgeTraversalEvent(nextInactiveEventIndex, nextEventIndex, inactiveEdges, segsToRemoveList, activeEdges, segsToAddList);
			}
//			System.out.println("active Edges size : " + activeEdges.size());
		}

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
		
		// Step 1: sample which of the coalLines will coalesce
		int coalLineIdx1 = Randomizer.nextInt(coalLines);
		
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
		for (NetworkEdge e : network.getEdges()) {
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

//		// if higher than root time, add root edge
//		if (rootEdge.childNode.getHeight() < currentTime) {
//			coexistingLineages.add(rootEdge);
//		}

		// add the coalLines that are not coalLineIdx1 to the coexistingLineages list
		count = 0;
		for (int i = 0; i < activeEdges.size(); i++) {
			if (activeEdges.get(i).hasSegments.isEmpty()) {
				if (count != coalLineIdx1) {
					coexistingLineages.add(activeEdges.get(i));
				}
				count++;
			}else{
				// also add active edges with segments
//						coexistingLineages.add(activeEdges.get(i));
			}
		}
//				System.out.println(activeEdges.size() + " " + currentTime);

		int coalLineIdx2 = Randomizer.nextInt(coexistingLineages.size());

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

//					System.out.println(segsToAddList.get(idx) + " " + segsToAddList.get(idx2));


			BitSet parentSegs = (BitSet) edge1.hasSegments.clone();
			parentSegs.or(edge2.hasSegments);

			coalescentNode.addChildEdge(edge2);
			coalescentNode.addChildEdge(edge1);
            NetworkEdge parentEdge = new NetworkEdge();
            parentEdge.hasSegments = new BitSet();
            coalescentNode.addParentEdge(parentEdge);

//                    System.out.println(activeEdges.size() + " " + idx + " " + idx2);

			activeEdges.remove(Math.max(idx, idx2));
			activeEdges.remove(Math.min(idx, idx2));
			segsToAddList.remove(Math.max(idx, idx2));
			segsToAddList.remove(Math.min(idx, idx2));

			activeEdges.add(parentEdge);
			segsToAddList.add(parentSegs);

			edgesAdded.add(parentEdge);


		}else if (activeEdges.contains(edge2)) {
//					System.out.println("Coalescence between coal and active " + currentTime + " " + edge2.isRootEdge());
			int idx2 = -1;
			for (int i = 0; i < activeEdges.size(); i++) {
				if (activeEdges.get(i) == edge2) {
					idx2 = i;
					break;
				}
			}
//					System.out.println(segsToAddList.get(idx) + " " + segsToAddList.get(idx2));
//					System.out.println(activeEdges.get(idx2).parentNode.getHeight());
//					System.out.println(activeEdges.get(idx2).childNode.getHeight() + " " + activeEdges.get(idx).childNode.getHeight());
//					System.out.println(activeEdges.get(idx2).hasSegments + " " + activeEdges.get(idx).hasSegments);
//					System.out.println(segsToAddList.get(idx2) + " " + segsToAddList.get(idx));
			exitOnError=true;
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


//					System.out.println(grandParent == null);
//
//					System.out.println(parentEdge.childNode.getChildEdges().get(0).parentNode.getHeight());
//					System.out.println(parentEdge.childNode.getChildEdges().get(1).parentNode.getHeight());


//					System.out.println(activeEdges.size());


			activeEdges.remove(Math.max(idx, idx2));
			activeEdges.remove(Math.min(idx, idx2));
			segsToAddList.remove(Math.max(idx, idx2));
			segsToAddList.remove(Math.min(idx, idx2));

			if (!parentSegsToAdd.isEmpty()){
				activeEdges.add(parentEdge);
				segsToAddList.add(parentSegsToAdd);
			}
//					System.out.println(parentSegsToAdd);
//					System.out.println(parentEdge.hasSegments);
//					System.out.println(parent.getHeight());
//					System.out.println(parentEdge.childNode.getChildEdges().get(0).parentNode.getHeight());
//							System.out.println(		" " + parentEdge.childNode.getChildEdges().get(1).parentNode.getHeight());

			edgesAdded.add(parentEdge);

//					System.out.println(activeEdges.size());

		}else {
//					System.out.println("Coalescence between coal line and existing lineage at time " + currentTime+ " " + edge2.isRootEdge());
//					System.out.println(edge2.hasSegments);
			NetworkNode grandParent = edge2.parentNode;
			if (grandParent != null) {
//						System.out.println(edge2.parentNode.getHeight() + " " + edge2.childNode.getHeight());
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
//						System.out.println("replace " + parentEdge.parentNode==null);
			}

			edgesAdded.add(parentEdge);
//					System.out.println(activeEdges.size());

		}

		// Reverse operation probabilities
		// 1. Probability of choosing which coal line to remove: 1/coalLines
		logHR += Math.log(1.0 / coalLines);

		// 2. Probability of choosing which partner lineage: 1/coexistingLineages.size()
		logHR += Math.log(1.0 / coexistingLineages.size());

		return logHR;
	}

	private double handleReassortmentEvent(List<NetworkEdge> activeEdges, List<NetworkEdge> inactiveEdges,
			List<BitSet> segsToAddList, List<BitSet> segsToRemoveList,
			double[] reassortmentObsProb, double totalReassortmentProb, double currentTime, List<NetworkEdge> edgesAdded) {
		double logHR = 0.0;

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

				// Sample which segments will reassort (each with 0.5 probability, conditional on observable reassortment)
				BitSet segsToReassort = segsToAddList.get(reassortEdgeIdx);
				BitSet segsToStay = new BitSet();
				BitSet segsToGo = new BitSet();

				// Sample segments with 0.5 probability each, conditional on at least one on each side
				if (activeEdges.get(reassortEdgeIdx).hasSegments.isEmpty()) {
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

				}

//				System.out.println("Segs to stay: " + segsToStay + ", Segs to go: " + segsToGo + " segs to reassort " + segsToReassort);
//				System.out.println(activeEdges.get(reassortEdgeIdx).hasSegments);

				// Step 2: create a new reassortment node and add it to the network
				NetworkNode newGrandParent = activeEdges.get(reassortEdgeIdx).parentNode;
				NetworkNode reassortmentNode = new NetworkNode();
				reassortmentNode.setHeight(currentTime);

				NetworkEdge parentToStay = new NetworkEdge();
				parentToStay.hasSegments = (BitSet) activeEdges.get(reassortEdgeIdx).hasSegments.clone();

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
				
				
				if (inactiveEdges.contains(activeEdges.get(reassortEdgeIdx))){
					int idx = inactiveEdges.indexOf(activeEdges.get(reassortEdgeIdx));
					activeEdges.get(reassortEdgeIdx).hasSegments.andNot(segsToRemoveList.get(idx));
					inactiveEdges.set(idx, parentToStay);
				}



				// replace the active edge with the parentToStay edge
				if (segsToStay.isEmpty()) {
					activeEdges.remove(reassortEdgeIdx);
					segsToAddList.remove(reassortEdgeIdx);
                } else {
    				activeEdges.set(reassortEdgeIdx, parentToStay);
    				segsToAddList.set(reassortEdgeIdx, segsToStay);
                }

				edgesAdded.add(parentToGo);
				edgesAdded.add(parentToStay);



				activeEdges.add(parentToGo);
				segsToAddList.add(segsToGo);


//				System.out.println(activeEdges.size());

		// Reverse operation probabilities
		// 1. Probability of choosing which edge had the reassortment: weighted by reassortmentObsProb
		logHR += Math.log(reassortmentObsProb[reassortEdgeIdx] / totalReassortmentProb);

		// 2. Probability of the segment split (conditional on observable reassortment)
		// Calculate total segments that were split
		BitSet allReassortSegs = (BitSet) segsToStay.clone();
		allReassortSegs.or(segsToGo);
		int k = allReassortSegs.cardinality();

		// For reverse: probability of splitting segsToStay and segsToGo
		if (activeEdges.get(reassortEdgeIdx).hasSegments.isEmpty()) {
			// Both sides must be non-empty (was a coal line)
			double prob = Math.pow(0.5, k) / (1.0 - 2.0 * Math.pow(0.5, k));
			logHR += Math.log(prob);
		} else {
			// Only segsToGo must be non-empty (already had segments)
			double prob = Math.pow(0.5, k) / (1.0 - Math.pow(0.5, k));
			logHR += Math.log(prob);
		}

		return logHR;
	}

	private void handleEdgeTraversalEvent(int nextInactiveEventIndex, int nextEventIndex,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList,
			List<NetworkEdge> activeEdges, List<BitSet> segsToAddList) {
		if (nextInactiveEventIndex != -1) {
			handleInactiveEdgeTraversal(nextInactiveEventIndex, inactiveEdges, segsToRemoveList, activeEdges);
		} else { // the edge is both an active and an inactive edge
			handleActiveEdgeTraversal(nextEventIndex, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList);
		}
	}

	private void handleInactiveEdgeTraversal(int idx, List<NetworkEdge> inactiveEdges, 
			List<BitSet> segsToRemoveList, List<NetworkEdge> activeEdges) {
//		System.out.println("handleInactiveEdgeTraversal");
		if (inactiveEdges.get(idx).parentNode.isReassortment()) {
			handleReassortment(activeEdges, segsToRemoveList, inactiveEdges, segsToRemoveList, inactiveEdges.get(idx));
//			handleInactiveReassortmentTraversal(idx, inactiveEdges, segsToRemoveList);
		} else {
			NetworkEdge sisterEdge = getSisterEdge(inactiveEdges.get(idx));
			handleCoalescence(activeEdges, segsToRemoveList, inactiveEdges, segsToRemoveList, inactiveEdges.get(idx), sisterEdge);
			
//			if (!activeEdges.contains(sisterEdge) && !inactiveEdges.contains(sisterEdge)) {
//				handleInactiveCoalescenceSisterNotInLists(idx, inactiveEdges, segsToRemoveList, sisterEdge);
//			} else if (inactiveEdges.contains(sisterEdge) && !activeEdges.contains(sisterEdge)) {
//				handleInactiveCoalescenceSisterInInactive(idx, inactiveEdges, segsToRemoveList, sisterEdge);
//			} else {
//				// missing implementation
//				throw new UnsupportedOperationException("Not yet implemented");
//			}
		}
	}

	private void handleActiveEdgeTraversal(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList) {
//		System.out.println("handleActiveEdgeTraversal");
//		System.out.println(segsToAddList.get(idx));
		if (activeEdges.get(idx).parentNode.isReassortment()) {
			handleReassortment(activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, activeEdges.get(idx));
//			handleActiveReassortmentTraversal(idx, activeEdges, segsToAddList);
		} else {
			NetworkEdge sisterEdge = getSisterEdge(activeEdges.get(idx));
			handleCoalescence(activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, activeEdges.get(idx), sisterEdge);
		}
	}

	private void handleReassortment(List<NetworkEdge> activeEdges, List<BitSet> segsToAddList,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, NetworkEdge edge) {

//		System.out.println("handleReassortment");

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
		
//		System.out.println("Reassortment - idx in active: " + idx + ", idx in inactive: " + iidx);
//		System.out.println("Edge height: " + edge.childNode.getHeight() + " -> " + edge.parentNode.getHeight());
//		System.out.println("Edge hasSegments: " + edge.hasSegments);
//		if (idx != -1) {
//			System.out.println("segsToAdd: " + segsToAddList.get(idx));
//		}
//		if (iidx != -1) {
//			System.out.println("segsToRemove: " + segsToRemoveList.get(iidx));
//		}
		
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
		}
		
		// now, check what to remove
		BitSet segsToRemoveLeft = new BitSet();
		BitSet segsToRemoveRight = new BitSet();


		// remove from inactive edges
		if (iidx != -1) {
			segsToRemoveLeft.or(segsToRemoveList.get(iidx));
			segsToRemoveRight.or(segsToRemoveList.get(iidx));
			
			edge.hasSegments.andNot(segsToRemoveList.get(iidx));
			inactiveEdges.remove(iidx);
			segsToRemoveList.remove(iidx);
			
			segsToRemoveLeft.and(edge.parentNode.getParentEdges().get(0).hasSegments);
			segsToRemoveRight.and(edge.parentNode.getParentEdges().get(1).hasSegments);
		}

		// remove from active edges
		if (idx != -1) {
			edge.hasSegments.or(segsToAddList.get(idx));
			activeEdges.remove(idx);
			segsToAddList.remove(idx);
		}

//		System.out.println("segsToGoLeft: " + segsToGoLeft + ", segsToGoRight: " + segsToGoRight);
//		System.out.println("segsToRemoveLeft: " + segsToRemoveLeft + ", segsToRemoveRight: " + segsToRemoveRight);
//		
		// now, check if we add
		if (!segsToGoLeft.isEmpty()) {
			activeEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToAddList.add(segsToGoLeft);
//			System.out.println("Added left parent to activeEdges");
		}
		if (!segsToGoRight.isEmpty()) {
			activeEdges.add(edge.parentNode.getParentEdges().get(1));
			segsToAddList.add(segsToGoRight);
//			System.out.println("Added right parent to activeEdges");
		}
		
		if (!segsToRemoveLeft.isEmpty()) {
			inactiveEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToRemoveList.add(segsToRemoveLeft);
//			System.out.println("Added left parent to inactiveEdges");
		}
		if (!segsToRemoveRight.isEmpty()) {
			inactiveEdges.add(edge.parentNode.getParentEdges().get(1));
			segsToRemoveList.add(segsToRemoveRight);
//			System.out.println("Added right parent to inactiveEdges");
		}
		
//		System.out.println("After handleReassortment - activeEdges size: " + activeEdges.size() + ", inactiveEdges size: " + inactiveEdges.size());
		
	}
	
	private void handleCoalescence(List<NetworkEdge> activeEdges, List<BitSet> segsToAddList,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, NetworkEdge edge,  NetworkEdge sisterEdge) {
		
//		System.out.println("handleCoalescence");
		
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
		
//		System.out.println("Coalescence - edge idx in active: " + idx1 + ", in inactive: " + iidx1);
//		System.out.println("Coalescence - sister idx in active: " + idx2 + ", in inactive: " + iidx2);
//		System.out.println("Edge height: " + edge.childNode.getHeight() + " -> " + edge.parentNode.getHeight());
//		System.out.println("Sister height: " + sisterEdge.childNode.getHeight() + " -> " + sisterEdge.parentNode.getHeight());
//		System.out.println("Edge hasSegments before: " + segsBefore1 + ", sister hasSegments before: " + segsBefore2);
//		System.out.println("segsToAdd1: " + segsToAdd1 + ", segsToAdd2: " + segsToAdd2);
//		System.out.println("segsToRemove1: " + segsToRemove1 + ", segsToRemove2: " + segsToRemove2);
		
		edge.hasSegments.or(segsToAdd1);
		sisterEdge.hasSegments.or(segsToAdd2);
		
		edge.hasSegments.andNot(segsToRemove1);
		sisterEdge.hasSegments.andNot(segsToRemove2);
		
		// now, remove
		segsToAdd1.andNot(segsBefore2);
		segsToAdd2.andNot(segsBefore1);
		
		segsToRemove1.andNot(sisterEdge.hasSegments);
		segsToRemove2.andNot(edge.hasSegments);
		
//		System.out.println(segsBefore1 + " " + segsBefore2);
		
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
		
//		System.out.println("After processing - edge hasSegments: " + edge.hasSegments + ", sister hasSegments: " + sisterEdge.hasSegments);
//		System.out.println("Combined segsToAdd: " + segsToAdd1 + ", combined segsToRemove: " + segsToRemove1);
		
		if (!segsToAdd1.isEmpty()) {
			activeEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToAddList.add(segsToAdd1);
//			System.out.println("Added parent to activeEdges");
		}
		if (!segsToRemove1.isEmpty()) {
			inactiveEdges.add(edge.parentNode.getParentEdges().get(0));
			segsToRemoveList.add(segsToRemove1);
//			System.out.println("Added parent to inactiveEdges");
		}
		
//		System.out.println("After handleCoalescence - activeEdges size: " + activeEdges.size() + ", inactiveEdges size: " + inactiveEdges.size());
	}

	private void handleActiveReassortmentTraversal(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList) {
		System.out.println("handleActiveReassortmentTraversal");
		NetworkEdge e = activeEdges.get(idx);
		// Step 2: sample which of the segsToAdd will have originated from which parent
		BitSet segsToAddLeft = new BitSet();
		BitSet segsToAddRight = new BitSet();
		for (int segIdx = segsToAddList.get(idx).nextSetBit(0); segIdx != -1; segIdx = segsToAddList.get(idx).nextSetBit(segIdx + 1)) {
			if (Randomizer.nextBoolean()) {
				segsToAddLeft.set(segIdx);
			} else {
				segsToAddRight.set(segIdx);
			}
		}

		e.hasSegments.or(segsToAddList.get(idx));

		// remove the current active edge and the segs to add
		activeEdges.remove(idx);
		segsToAddList.remove(idx);

//					System.out.println(segsToAddLeft+" " + segsToAddRight);

		if (!segsToAddLeft.isEmpty()) {
			activeEdges.add(e.parentNode.getParentEdges().get(0));
			segsToAddList.add(segsToAddLeft);
		}
		if (!segsToAddRight.isEmpty()) {
			activeEdges.add(e.parentNode.getParentEdges().get(1));
			segsToAddList.add(segsToAddRight);
		}
	}

	private void handleActiveCoalescenceNotInInactive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, 
			NetworkEdge sisterEdge) {
		if (!activeEdges.contains(sisterEdge)) {
			handleActiveCoalescenceSisterNotInActive(idx, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, sisterEdge);
		} else if (activeEdges.contains(sisterEdge)) {
			handleActiveCoalescenceSisterInActive(idx, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, sisterEdge);
		}
	}

	private void handleActiveCoalescenceSisterNotInActive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, 
			NetworkEdge sisterEdge) {
		System.out.println("handleActiveCoalescenceSisterNotInActive");
		activeEdges.get(idx).hasSegments.or(segsToAddList.get(idx));

		segsToAddList.get(idx).andNot(sisterEdge.hasSegments);

		
		if (inactiveEdges.contains(sisterEdge)) {

			int inactiveIdx = -1;
			for (int i = 0; i < inactiveEdges.size(); i++) {
				if (inactiveEdges.get(i) == sisterEdge) {
					inactiveIdx = i;
					break;
				}
			}

			BitSet segsToRemoveParent = (BitSet) segsToRemoveList.get(inactiveIdx).clone();
			segsToRemoveParent.andNot(activeEdges.get(idx).hasSegments);

			inactiveEdges.get(inactiveIdx).hasSegments.andNot(segsToRemoveList.get(inactiveIdx));

			if (segsToRemoveParent.isEmpty()) {
				inactiveEdges.remove(inactiveIdx);
				segsToRemoveList.remove(inactiveIdx);
			} else {
				inactiveEdges.set(inactiveIdx,
						inactiveEdges.get(inactiveIdx).parentNode.getParentEdges().get(0));
				segsToRemoveList.set(inactiveIdx, segsToRemoveParent);
			}									
		}
		if (segsToAddList.get(idx).isEmpty()) {
			activeEdges.remove(idx);
			segsToAddList.remove(idx);
		} else {
			activeEdges.set(idx, activeEdges.get(idx).parentNode.getParentEdges().get(0));
		}
	}

	private void handleActiveCoalescenceSisterInActive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, 
			NetworkEdge sisterEdge) {
		System.out.println("handleActiveCoalescenceSisterInActive");
		// find idx2 for sisterEdge
		int idx2 = -1;
		for (int i = 0; i < activeEdges.size(); i++) {
			if (activeEdges.get(i) == sisterEdge) {
				idx2 = i;
				break;
			}
		}

		BitSet segsToAddParent = (BitSet) segsToAddList.get(idx).clone();
		segsToAddParent.or(segsToAddList.get(idx2));
		segsToAddParent.andNot(activeEdges.get(idx).hasSegments);
		segsToAddParent.andNot(activeEdges.get(idx2).hasSegments);

		activeEdges.get(idx).hasSegments.or(segsToAddList.get(idx));
		activeEdges.get(idx2).hasSegments.or(segsToAddList.get(idx2));


		
		if (inactiveEdges.contains(sisterEdge)) {
			int inactiveIdx = -1;
			for (int i = 0; i < inactiveEdges.size(); i++) {
				if (inactiveEdges.get(i) == sisterEdge) {
					inactiveIdx = i;
					break;
				}
			}
			
			BitSet segsToRemoveParent = (BitSet) segsToRemoveList.get(inactiveIdx).clone();
			segsToRemoveParent.andNot(activeEdges.get(idx).hasSegments);
			
			activeEdges.get(idx).hasSegments.andNot(segsToRemoveList.get(inactiveIdx));
			
			if (segsToRemoveList.get(inactiveIdx).isEmpty()) {
				inactiveEdges.remove(inactiveIdx);
				segsToRemoveList.remove(inactiveIdx);
			} else {
				inactiveEdges.set(inactiveIdx, inactiveEdges.get(inactiveIdx).parentNode.getParentEdges().get(0));
				segsToRemoveList.set(inactiveIdx, segsToRemoveParent);
			}

			
		}
		
		// add parent edge
		activeEdges.add(activeEdges.get(idx).parentNode.getParentEdges().get(0));
		segsToAddList.add(segsToAddParent);

		// remove both edges
		activeEdges.remove(Math.max(idx, idx2));
		activeEdges.remove(Math.min(idx, idx2));
		segsToAddList.remove(Math.max(idx, idx2));
		segsToAddList.remove(Math.min(idx, idx2));
	}

	private void handleActiveCoalescenceInInactive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, 
			NetworkEdge sisterEdge) {
		if (!activeEdges.contains(sisterEdge)) {
			handleActiveInInactiveSisterNotInActive(idx, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, sisterEdge);
		} else {
			// both edges are in active and inactive edges
			activeEdges.get(idx).hasSegments.or(segsToAddList.get(idx));
			segsToAddList.get(idx).andNot(sisterEdge.hasSegments);

			
			// missing implementation
			throw new UnsupportedOperationException("Not yet implemented2");
		}
	}

	private void handleActiveInInactiveSisterNotInActive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, 
			NetworkEdge sisterEdge) {
		System.out.println("handleActiveInInactiveSisterNotInActive");
		activeEdges.get(idx).hasSegments.or(segsToAddList.get(idx));
		segsToAddList.get(idx).andNot(sisterEdge.hasSegments);

		if (inactiveEdges.contains(sisterEdge)) {
			handleActiveInInactiveBothSistersInInactive(idx, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList, sisterEdge);
		} else {
			handleActiveInInactiveOnlyCurrentInInactive(idx, activeEdges, segsToAddList, inactiveEdges, segsToRemoveList);
		}
		
		
		if (segsToAddList.get(idx).isEmpty()) {
			activeEdges.remove(idx);
			segsToAddList.remove(idx);
		} else {
			activeEdges.set(idx, activeEdges.get(idx).parentNode.getParentEdges().get(0));
		}
	}

	private void handleActiveInInactiveBothSistersInInactive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList, 
			NetworkEdge sisterEdge) {
		System.out.println("handleActiveInInactiveBothSistersInInactive");
		// inadctive Edges contains both edges
		int inactiveIdx1 = -1;
		for (int i = 0; i < inactiveEdges.size(); i++) {
			if (inactiveEdges.get(i) == sisterEdge) {
				inactiveIdx1 = i;
				break;
			}
		}
		int inactiveIdx2 = -1;
		for (int i = 0; i < inactiveEdges.size(); i++) {
			if (inactiveEdges.get(i) == activeEdges.get(idx)) {
				inactiveIdx2 = i;
				break;
			}
		}
		
		BitSet segsToRemoveParent = (BitSet) segsToRemoveList.get(inactiveIdx1).clone();
		segsToRemoveParent.or(segsToRemoveList.get(inactiveIdx2));
		segsToRemoveParent.andNot(inactiveEdges.get(inactiveIdx1).hasSegments);
		segsToRemoveParent.andNot(inactiveEdges.get(inactiveIdx2).hasSegments);

		inactiveEdges.get(inactiveIdx1).hasSegments.andNot(segsToRemoveList.get(inactiveIdx1));
		inactiveEdges.get(inactiveIdx2).hasSegments.andNot(segsToRemoveList.get(inactiveIdx2));


		// add parent edge
		if (!segsToRemoveParent.isEmpty()) {
			inactiveEdges.add(inactiveEdges.get(idx).parentNode.getParentEdges().get(0));
			segsToRemoveList.add(segsToRemoveParent);
		}

		// remove both edges
		inactiveEdges.remove(Math.max(inactiveIdx1, inactiveIdx2));
		inactiveEdges.remove(Math.min(inactiveIdx1, inactiveIdx2));
		segsToRemoveList.remove(Math.max(inactiveIdx1, inactiveIdx2));
		segsToRemoveList.remove(Math.min(inactiveIdx1, inactiveIdx2));
	}

	private void handleActiveInInactiveOnlyCurrentInInactive(int idx, List<NetworkEdge> activeEdges, 
			List<BitSet> segsToAddList, List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList) {
		System.out.println("handleActiveInInactiveOnlyCurrentInInactive");
		int inactiveIdx = -1;
		for (int i = 0; i < inactiveEdges.size(); i++) {
			if (inactiveEdges.get(i) == activeEdges.get(idx)) {
				inactiveIdx = i;
				break;
			}
		}
		
		BitSet segsToRemoveParent = (BitSet) segsToRemoveList.get(inactiveIdx).clone();
		segsToRemoveParent.andNot(getSisterEdge(activeEdges.get(idx)).hasSegments);
		
		activeEdges.get(idx).hasSegments.andNot(segsToRemoveList.get(inactiveIdx));
		
		if (segsToRemoveList.get(inactiveIdx).isEmpty()) {
			inactiveEdges.remove(inactiveIdx);
			segsToRemoveList.remove(inactiveIdx);
		} else {
			inactiveEdges.set(inactiveIdx, inactiveEdges.get(inactiveIdx).parentNode.getParentEdges().get(0));
			segsToRemoveList.set(inactiveIdx, segsToRemoveParent);
		}
	}

	/**
	 * Clean up empty reassortment edges.
	 * Uses the same pattern as EmptyEdgesNetworkOperator.RemoveAllEmptyNetworkSegments()
	 *
	 * @return log probability contribution for reverse operation
	 */
	protected double cleanEmptyEdgesTopDown() {
		double logHR = 0.0;

		// Create a local list of network edges to track during removal
		List<NetworkEdge> localNetworkEdges = new ArrayList<>(network.getEdges());

		List<Integer> removableEdges = new ArrayList<>();
		for (int i = 0; i < localNetworkEdges.size(); i++) {
			NetworkEdge edge = localNetworkEdges.get(i);
			if (!edge.isRootEdge() && edge.childNode.isReassortment()
					&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() == 0) {
				removableEdges.add(i);
			}
		}

		// Remove all empty edges iteratively
		while (removableEdges.size() > 0) {
			// TODO: For proper Hastings ratio, would need:
			// logHR -= Math.log(1.0 / removableEdges.size());

			int edgeInd = 0; // For now, deterministically remove first one (can randomize later)

			logHR += removeEmptyReassortmentEdgeAdapted(localNetworkEdges, localNetworkEdges.get(removableEdges.get(edgeInd)));

			if (logHR == Double.NEGATIVE_INFINITY)
				return Double.NEGATIVE_INFINITY;

			// Rebuild the list of removable edges
			removableEdges = new ArrayList<>();
			for (int i = 0; i < localNetworkEdges.size(); i++) {
				NetworkEdge edge = localNetworkEdges.get(i);
				if (!edge.isRootEdge() && edge.childNode.isReassortment()
						&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() == 0) {
					removableEdges.add(i);
				}
			}
		}

		// TODO: Reverse probability - Poisson probability of adding n empty edges
		// logHR += Math.log(Math.pow(lambda, nEmptyEdgesRemoved)) - lambda - logFactorial(nEmptyEdgesRemoved);

		return logHR;
	}

	/**
	 * Remove an empty reassortment edge, following the exact pattern from EmptyEdgesNetworkOperator.
	 * Updates the localNetworkEdges list to track removed edges.
	 */
	private double removeEmptyReassortmentEdgeAdapted(List<NetworkEdge> localNetworkEdges, NetworkEdge edgeToRemove) {
		double logHR = 0.0;

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

		if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
			network.setRootEdge(secondEdgeToExtend);
		} else {
			NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
			NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
			secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
			secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);
			secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
		}

		return logHR;
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

