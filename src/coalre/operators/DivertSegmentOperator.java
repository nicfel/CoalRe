package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.distribution.CoalescentWithReassortment;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertSegmentOperator extends EmptyEdgesNetworkOperator {

	public Input<Boolean> divertOneSegmentInput = new Input<>("divertOneSegment",
			"If true, only one segment is diverted", false);
	
    public Input<CoalescentWithReassortment> coalescentDistrInput = new Input<>("coalescentWithReassortment",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);

	
	
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
//
//		List<NetworkEdge> sourceEdges = networkEdges.stream().filter(e -> e.childNode.isReassortment())
//				.filter(e -> e.hasSegments.cardinality() > 0).collect(Collectors.toList());

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

//
//		int reverseSourceEdgeCount = (int) (networkEdges.stream().filter(e -> e.childNode.isReassortment())
//				.filter(e -> e.hasSegments.cardinality() > 0).count());

		logHR += Math.log(1.0 / reverseSourceEdgeCount);

		return logHR;
	}

	protected double divertSegments(NetworkEdge destEdge, NetworkEdge sourceEdge, BitSet segsToDivert) {
		double logHR = 0.0;

//		System.out.println(network.getExtendedNewickVerbose());
		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		
//		System.out.println(network.getExtendedNewickVerbose(5));
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
						
						
//						System.out.println(network.getExtendedNewickVerbose());
//						System.out.println(newGrandParent.getHeight());
//						System.out.println(i);
//						System.out.println(segmentTrees.get(i) +";");

						newGrandParent.removeChild(newSibling);
						oldParent.removeChild(oldSibling);
						newGrandParent.addChild(oldSibling);
						
						oldParent.addChild(newSibling);
						oldParent.setHeight(newNodeHeights[i]);
						
//						System.out.println(segmentTrees.get(i) +";");
						
//						System.exit(0);
						// TODO
//						System.out.println("do nothing "  + oldParent.getHeight() + " " + newNodeHeights[i]);
						
//						return true;
//						throw new IllegalArgumentException("do nothing " + oldParent.getHeight() + " " + newNodeHeights[i]);
						
					}else {
						
	//					System.out.println(oldGrandParent.getHeight());
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
							
//							
	
						}else {
//							System.out.println("basicmove " + newGrandParent.getHeight() + " " + oldGrandParent.getHeight() + " " + oldParent.getHeight() + " " + newNodeHeights[i]);
							newGrandParent.removeChild(newSibling);
							oldParent.removeChild(oldSibling);
							newGrandParent.addChild(oldParent);
							oldParent.addChild(newSibling);
							oldParent.setHeight(newNodeHeights[i]);
						}
					}
				}
//
//				
//				
//				
//				System.out.println(segmentTrees.get(i) +";");
//				System.out.println("Updated segment tree " + i);
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
						
//						System.out.println(segIdx + " " + edge.childNode.getHeight());

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
//        	if (segmentTrees.size() > 0) {
//	        	BitSet segsToRemoveLeft = (BitSet) segsToRemove.clone();
//	        	segsToRemoveLeft.and(getSisterEdge(edge).hasSegments);
//	        	// loop over the segsToRemove Left and update the Tree Nodes
//	        	while (segsToRemoveLeft.nextSetBit(0) != -1) {
//					int segIdx = segsToRemoveLeft.nextSetBit(0);
//					// get the corresponding tree node
//					Node treeNode = segmentTrees.get(segIdx).getNode(edge.parentNode.segmentIndices[segIdx]);
//					treeNodes[segIdx] = treeNode;
//					// get the child node index that is does not have number treeNodesIn[segIdx]
//					Node otherChild = treeNode.getChild(0).getNr;
//							
//				}	
//        	}
			segsToRemove.andNot(getSisterEdge(edge).hasSegments);

			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);

		}

		return logP;
	}
	
	protected double addSegmentsToAncestorsThroughSimulation(NetworkEdge edge, BitSet segsToAdd) {
		double logP = 0.0;

		List<NetworkEdge> currentEdges = new ArrayList<>();
		List<Boolean> alreadyExists = new ArrayList<>();
		List<BitSet> segmentsToAddList = new ArrayList<>();
		currentEdges.add(edge);
		alreadyExists.add(true);
		segmentsToAddList.add(segsToAdd);

		double currentTime = edge.childNode.getHeight();

		// now, simulate the reassortment and coalescence of new edges.
		// existing edges can only have reassortment but not coalescent events
		while (currentEdges.size() > 0) {

			// Process each lineage to calculate rates and observation probabilities
			int coalLins = 0; // Free lineages that can coalesce
			double[] observableReassortmentProbs = new double[currentEdges.size()];
			double totalReassortmentObsProb = 0.0; // Total observable reassortment probability
			double timeToNextNetworkEvent = Double.POSITIVE_INFINITY;

			for (int i = 0; i < currentEdges.size(); i++) {
				if (alreadyExists.get(i)) {
					// Existing edge: find time to next network event
					NetworkEdge e = currentEdges.get(i);
					double timeToParent = e.parentNode.getHeight() - currentTime;
					if (timeToParent < timeToNextNetworkEvent) {
						timeToNextNetworkEvent = timeToParent;
					}

					// For existing edges, only count observable reassortments involving NEW segments
					// A reassortment is observable if it splits the new segments from existing segments
					// OR if it splits the new segments among themselves
					// For now, we only care about the new segments
					BitSet newSegs = segmentsToAddList.get(i);
					int numNewSegs = newSegs.cardinality();
					// Observation probability = 1 - 2*(0.5)^n
					// This is the prob that reassortment creates two non-empty groups from new segments
					observableReassortmentProbs[i] = 1.0 - 2.0 * Math.pow(0.5, numNewSegs);
					// correct for any event taht includes an existing segment
					// TODO correct for existing segments
				} else {
					// Free lineage: can coalesce and can have observable reassortment
					coalLins++;

					// For free lineages, all segments can reassort
					BitSet allSegs = segmentsToAddList.get(i);
					int numSegs = allSegs.cardinality();
					observableReassortmentProbs[i] = 1.0 - 2.0 * Math.pow(0.5, numSegs);
					totalReassortmentObsProb += observableReassortmentProbs[i];
				}

			// Calculate time to next coalescent event (only for free lineages)
			double timeToNextCoal = Double.POSITIVE_INFINITY;
			if (coalLins >= 1) {
				// TODO: include the other lineages existing in the network for coalescence
				
				double currentTransformedTime = coalescentDistrInput.get().populationFunction.getIntensity(currentTime);
				double transformedTimeToNextCoal = Randomizer.nextExponential(0.5 * coalLins);
				timeToNextCoal = coalescentDistrInput.get().populationFunction.getInverseIntensity(
						transformedTimeToNextCoal + currentTransformedTime) - currentTime;
			}

			// Calculate time to next reassortment event (using total observable reassortment probability)
			double timeToNextReassortment = Double.POSITIVE_INFINITY;
			if (totalReassortmentObsProb > 0) {
				if (coalescentDistrInput.get().timeVaryingReassortmentRatesInput.get() != null) {
					double currentTransformedReaTime = coalescentDistrInput.get().timeVaryingReassortmentRatesInput.get()
							.getIntensity(currentTime);
					double transformedTimeToNextRea = Randomizer.nextExponential(totalReassortmentObsProb);
					timeToNextReassortment = coalescentDistrInput.get().timeVaryingReassortmentRatesInput.get()
							.getInverseIntensity(transformedTimeToNextRea + currentTransformedReaTime) - currentTime;
				} else {
					timeToNextReassortment = Randomizer.nextExponential(
							totalReassortmentObsProb * coalescentDistrInput.get().reassortmentRateInput.get().getArrayValue());
				}
			}

			// get which of the three events happens next
			double minTime = Math.min(timeToNextCoal, Math.min(timeToNextReassortment, timeToNextNetworkEvent));	
			currentTime += minTime;
			// now process the event
			if (minTime == timeToNextCoal) {
				// coalescent event between two free lineages	
			} else if (minTime == timeToNextReassortment) {
				// reassortment event
			} else {
				// network event
			}	

			break; // Temporary - remove when implementing event handling
		}

		return logP;
	}



	/**
	 * Add segments to this edge and ancestors.
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

		edge.hasSegments.or(segsToAdd);

		if (edge.isRootEdge())
			return logP;

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

			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
		}

		return logP;
	}

	
	protected double addSegmentToNewAncestor(NetworkEdge edge, BitSet segsToAdd) {
		// sample time to next coalescent or reassortment event
        double currentTransformedTime = populationFunction.getIntensity(currentTime);
        double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
        double timeToNextCoal = populationFunction.getInverseIntensity(
                transformedTimeToNextCoal + currentTransformedTime) - currentTime;

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

        // Let n be the total number of possible segments in sourceSegments,
        // and k be the number actually chosen in the subset.
        int n = sourceSegments.cardinality();
        int k = chosenSubset.cardinality();


        // Fully correct conditional log-probability:
        //    log( p^k * (1-p)^(n-k) )  -  log( 1 - p^n - (1-p)^n )
        // =  k*log(p) + (n-k)*log(1-p) - log( 1 - p^n - (1-p)^n )
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

                // Now pick each segment with probability 0.1
                if (Randomizer.nextDouble() < binomProb) {
                    destSegments.set(segIdx);
                }
            }

        } while (destSegments.cardinality() == 0
                || destSegments.cardinality() == sourceSegments.cardinality());
//        System.out.println(destSegments + " " + sourceSegments);
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
