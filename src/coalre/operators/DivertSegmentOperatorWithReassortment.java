package coalre.operators;

import beast.base.core.Function;
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

	private Function reassortmentRate;
	private CoalescentWithReassortment coalescentDistr;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		if (coalescentDistrInput.get() != null) {
			coalescentDistr = coalescentDistrInput.get();
			reassortmentRate = coalescentDistr.reassortmentRateInput.get();
		}
	}

	@Override
	public double networkProposal() {
		double logHR = 0.0;
		
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
		for (int i = 0; i < networkEdges.size(); i++) {
			if (networkEdges.get(i).childNode.isReassortment() && networkEdges.get(i).hasSegments.cardinality() > 0) {
				reverseSourceEdgeCount++;
			}
		}

		logHR += Math.log(1.0 / reverseSourceEdgeCount);

		return logHR;
	}

	protected double divertSegments(NetworkEdge destEdge, NetworkEdge sourceEdge, BitSet segsToDivert) {
		double logHR = 0.0;

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		
		getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);

		logHR -= addSegmentsToAncestors(destEdge, segsToDivert);
		logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert);
		
		if (reconnectSegmentTrees(treeChildNodeList, destEdge, segsToDivert))
			return Double.NEGATIVE_INFINITY;

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
	 * Remove segments from this edge and ancestors.
	 *
	 * @param edge         edge at which to start removal
	 * @param segsToRemove segments to remove from edge and ancestors
	 * @return log probability of reverse operation
	 */
	protected double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
		double logP = 0.0;

		segsToRemove = (BitSet) segsToRemove.clone();
		segsToRemove.and(edge.hasSegments);

		if (segsToRemove.isEmpty())
			return logP;

		edge.hasSegments.andNot(segsToRemove);

		if (edge.isRootEdge())
			return logP;

		if (edge.parentNode.isReassortment()) {

			logP += Math.log(0.5) * segsToRemove.cardinality();

			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);
			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), segsToRemove);

		} else {
			segsToRemove.andNot(getSisterEdge(edge).hasSegments);

			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);

		}

		return logP;
	}

	/**
	 * Add segments to this edge and ancestors, potentially sampling new reassortment events.
	 *
	 * @param edge      edge at which to start addition
	 * @param segsToAdd segments to add to the edge and ancestors
	 * @return log probability of operation
	 */
	protected double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd) {
		double logP = 0.0;

		segsToAdd = (BitSet) segsToAdd.clone();
		segsToAdd.andNot(edge.hasSegments);

		if (segsToAdd.isEmpty())
			return logP;

		// do this in the loop, so we can keep track of the segments that are added to the edge
		// edge.hasSegments.or(segsToAdd);

		if (edge.isRootEdge())
			return logP;

		// startwith keeping a list of active edges and the segs to Add
		List<NetworkEdge> activeEdges = new ArrayList<>();
		List<BitSet> segsToAddList = new ArrayList<>();
		List<Boolean> existingEdges = new ArrayList<>(); // edge that can't coalesce, because we just follow it's history
		activeEdges.add(edge);
		segsToAddList.add(segsToAdd);
		existingEdges.add(true);

		// Track the current time (start at the bottom of the lowest edge)
		double currentTime = edge.childNode.getHeight();


		// now, sample the time to the next reassortment or coalescent event on the active edges,
		// i.e. essentially simulate the history of these edges and segments until there is nothing left. The
		// goal is to avoid the creation of empty
		while (!activeEdges.isEmpty()) {
			int coalLines = 0;
			double[] reassortmentObsProb = new double[activeEdges.size()];
			double totalReassortmentProb = 0.0;
			double timeToNextEdgeEvent = Double.POSITIVE_INFINITY;
			for (int i = 0; i < activeEdges.size(); i++) {
				if (existingEdges.get(i)) {
					int nSegs = activeEdges.get(i).hasSegments.cardinality() + segsToAddList.get(i).cardinality();
					// now, calculate the observation probability of a reassortment event, that does not involve any of the segments that are already on the edge
					double obsProb = 1.0 - 2.0 * Math.pow(0.5, nSegs);
					obsProb -= 
					reassortmentObsProb[i] = nSegs * reassortmentRate.getArrayValue() * obsProb;
					totalReassortmentProb += reassortmentObsProb[i];
					if (timeToNextEdgeEvent > activeEdges.get(i).parentNode.getHeight()) {
						timeToNextEdgeEvent = activeEdges.get(i).parentNode.getHeight();
					}
				}else{
					// For newly created edges (not existing), they can coalesce
					coalLines++;
					// Calculate reassortment obs prob for the segments being added
					int nSegs = segsToAddList.get(i).cardinality();
					double obsProb = 1.0 - 2.0 * Math.pow(0.5, nSegs);
					reassortmentObsProb[i] = nSegs * reassortmentRate.getArrayValue() * obsProb;
					totalReassortmentProb += reassortmentObsProb[i];
				}
			}

			
			// Sample time to next reassortment
			double timeToNextReassortment = Double.POSITIVE_INFINITY;
			if (totalReassortmentProb > 0) {
				timeToNextReassortment = Randomizer.nextExponential(totalReassortmentProb);
			}
			
			// Sample time to next coalescence using network intervals approach
			// Like AddRemoveReassortmentCoalescent, use the network event list to track lineages
			double timeToNextCoalescence = Double.POSITIVE_INFINITY;
			double coalRate = 0.0;
			if (coalLines >= 1) {
				// Get network events to determine lineages at different times
				List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();
				
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
						int totalLineages = prevEvent.lineages + coalLines;
						coalRate = 0.5 * totalLineages * (totalLineages - 1);
						
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
					int totalLineages = coalLines;  // Only new lineages left
					coalRate = 0.5 * totalLineages * (totalLineages - 1);
					
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
			
			// Update current time
			currentTime += timeUntilNextEvent;
			
			// Determine what type of event occurred
			if (timeUntilNextEvent == timeToNextCoalescence) {
				// Coalescence event - TODO: implement
				break; // Placeholder
			} else if (timeUntilNextEvent == timeToNextReassortment) {
				// Reassortment event - TODO: implement
				break; // Placeholder
			} else {
				// Edge event - reached parent node of an existing edge
				// TODO: implement
				break; // Placeholder
			}
		}

		if (edge.parentNode.isReassortment()) {

			BitSet segsToAddLeft = new BitSet();
			BitSet segsToAddRight = new BitSet();

			for (int segIdx = segsToAdd.nextSetBit(0); segIdx != -1; segIdx = segsToAdd.nextSetBit(segIdx + 1)) {
				if (Randomizer.nextBoolean())
					segsToAddLeft.set(segIdx);
				else
					segsToAddRight.set(segIdx);

				logP += Math.log(0.5);
			}

			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);

		} else {
			// Sample potential reassortment event on this edge
			if (reassortmentRate != null && segsToAdd.cardinality() > 0) {
				double edgeLength = edge.parentNode.getHeight() - edge.childNode.getHeight();
				
				// Sample time to reassortment for the segments being added
				// Rate is proportional to the number of segments being added
				double rate = segsToAdd.cardinality() * reassortmentRate.getArrayValue();
				double timeToReassortment = Randomizer.nextExponential(rate);
				
				// Check if reassortment occurs on this edge
				if (timeToReassortment < edgeLength) {
					// Create reassortment event
					double reassortmentTime = edge.childNode.getHeight() + timeToReassortment;
					
					// Split segments at reassortment
					BitSet segsToAddLeft = new BitSet();
					BitSet segsToAddRight = new BitSet();
					
					for (int segIdx = segsToAdd.nextSetBit(0); segIdx != -1; segIdx = segsToAdd.nextSetBit(segIdx + 1)) {
						if (Randomizer.nextBoolean()) {
							segsToAddLeft.set(segIdx);
						} else {
							segsToAddRight.set(segIdx);
						}
						logP += Math.log(0.5);
					}
					
					// Only create observable reassortment (both sides have segments)
					if (segsToAddLeft.cardinality() > 0 && segsToAddRight.cardinality() > 0) {
						// Log probability of sampling this reassortment time
						logP += -rate * timeToReassortment + Math.log(rate);
						
						// Create new reassortment node
						NetworkNode reassortmentNode = new NetworkNode();
						reassortmentNode.setHeight(reassortmentTime);
						
						// Split the edge at the reassortment point
						NetworkNode parentNode = edge.parentNode;
						parentNode.removeChildEdge(edge);
						reassortmentNode.addChildEdge(edge);
						
						// Create two new edges going up from the reassortment
						NetworkEdge leftEdge = new NetworkEdge();
						NetworkEdge rightEdge = new NetworkEdge();
						
						reassortmentNode.addParentEdge(leftEdge);
						reassortmentNode.addParentEdge(rightEdge);
						parentNode.addChildEdge(leftEdge);
						parentNode.addChildEdge(rightEdge);
						
						// Set segments on the new edges
						leftEdge.hasSegments = (BitSet) edge.hasSegments.clone();
						leftEdge.hasSegments.or(segsToAddLeft);
						
						rightEdge.hasSegments = (BitSet) edge.hasSegments.clone();
						rightEdge.hasSegments.or(segsToAddRight);
						
						// Add edges to network edge list
						networkEdges.add(leftEdge);
						networkEdges.add(rightEdge);
						
						// Continue recursively on both branches
						logP += addSegmentsToAncestors(leftEdge, segsToAddLeft);
						logP += addSegmentsToAncestors(rightEdge, segsToAddRight);
						
						return logP;
					} else {
						// Unobservable reassortment, proceed without creating it
						// Log probability of not having observable reassortment in this interval
						logP += -rate * edgeLength;
					}
				} else {
					// No reassortment on this edge
					// Log probability of not having reassortment in this interval
					logP += -rate * edgeLength;
				}
			}
			
			// No reassortment created, continue normally
			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
		}

		return logP;
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

