package coalre.operators;

import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Operator that removes a segment from parental edges and re-adds it starting from the same edge.
 */
public class RemoveSegmentAndResimulate extends DivertSegmentAndResimulate {


	@Override
	public void initAndValidate() {
		super.initAndValidate();
		if (coalescentDistrInput.get() != null) {
			coalescentDistr = coalescentDistrInput.get();
		}
	}

	@Override
	public double networkProposal() {
		double logHR = 0.0;
		// pick a segment at random
		int segmentIdx = Randomizer.nextInt(network.getSegmentCount());
		System.out.println(network.getExtendedNewickVerbose(segmentIdx));
		// get the corresponding tree node
		Node n = segmentTrees.get(segmentIdx).getRoot();
		// get the two children of the tree node
		Node leftChild = n.getLeft();
		Node rightChild = n.getRight();
		// now, find th two edges in the network that correspond to the children
		NetworkEdge edge1 = null;
		NetworkEdge edge2 = null;
		for (NetworkEdge e : networkEdges) {
			if ((e.childNode.isCoalescence() 
					&& e.childNode.getChildEdges().get(0).hasSegments.get(segmentIdx) 
					&& e.childNode.getChildEdges().get(1).hasSegments.get(segmentIdx)) ||
					(e.childNode.isLeaf() && e.hasSegments.get(segmentIdx))){
				int nodeIdx = e.childNode.segmentIndices[segmentIdx];
				if (nodeIdx == leftChild.getNr()) {
					edge1 = e;
				} else if (nodeIdx == rightChild.getNr()) {
					edge2 = e;
				}
			}			
		}
		
		
		
		
		// now, set edge1 to be the one with the higher childNode height
		if (edge1.childNode.getHeight() < edge2.childNode.getHeight()) {
			NetworkEdge temp = edge1;
			edge1 = edge2;
			edge2 = temp;
		}
		// now, trace edge2 up until the edge crosses the height of edge1
		while (edge2.parentNode.getHeight() > edge1.parentNode.getHeight()) {
			if (edge2.parentNode.isCoalescence()) {
				edge2 = edge2.parentNode.getParentEdges().get(0);
			} else {
				if (edge2.parentNode.getParentEdges().get(0).hasSegments.get(segmentIdx)) {
					edge2 = edge2.parentNode.getParentEdges().get(0);
				} else {
					edge2 = edge2.parentNode.getParentEdges().get(1);
				}
			}
		}

		// now, start the segs to divert, the active and inactive edges
		List<NetworkEdge> activeEdges = new ArrayList<>();
		List<BitSet> segsToAddList = new ArrayList<>();
		List<NetworkEdge> inactiveEdges = new ArrayList<>();
		List<BitSet> segsToRemoveList = new ArrayList<>();
		activeEdges.add(edge1);
		activeEdges.add(edge2);
		inactiveEdges.add(edge1);
		inactiveEdges.add(edge2);

		BitSet segsToDivert = new BitSet();
		segsToDivert.set(segmentIdx);
		segsToAddList.add(segsToDivert);
		segsToAddList.add(segsToDivert);
		segsToRemoveList.add(segsToDivert);
		segsToRemoveList.add(segsToDivert);

		coalescentDistr.intervals.update();
		List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();

		double currentTime = edge1.childNode.getHeight();
		logHR += divertSegmentsToAncestors(new ArrayList<>(), inactiveEdges, new ArrayList<>(), segsToRemoveList, 
				currentTime, networkEventList, false, false);
		System.out.println(network.getExtendedNewickVerbose(segmentIdx) + "\n");

		logHR += divertSegmentsToAncestors(activeEdges, new ArrayList<>(), segsToAddList, new ArrayList<>(), 
				currentTime, networkEventList, false, false);
		
		System.out.println(network.getExtendedNewickVerbose(segmentIdx) + "\n");
		System.exit(0);

		return logHR;
	}

	protected double divertSegmentsToAncestors(List<NetworkEdge> activeEdges, 
			List<NetworkEdge> inactiveEdges, 
			List<BitSet> segsToAddList, 
			List<BitSet> segsToRemoveList, 
			double currentTime, List<NetworkEvent> networkEventList, 
			boolean zeroActiveReassortment, boolean zeroInactiveReassortment) {
		double logHR = 0.0;
		
		NetworkEdge rootEdge = network.getRootEdge();
		List<NetworkEdge> edgesAdded = new ArrayList<>();
		
		boolean exitafternext = false;
		double timeToRoot = network.getRootEdge().childNode.getHeight();
		
		int round = 0;

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

					reassortmentObsProb[i] = 2.0 * (Math.pow(0.5, m)*(1-Math.pow(0.5, k)));

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
	                double transformedTimeToNextRea = Randomizer.nextExponential(reassortmentObsProb.length);
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
			boolean coalAboveRoot = false;
			

			if (coalLines >= 1 || negLines >= 1) {
				
				double checkTime = currentTime;
				NetworkEvent prevEvent = null;
				
				// Calculate the earliest stopping time (reassortment or edge event)
				double stopTime = Math.min(Math.min(currentTime + timeToNextReassortment, timeToNextEdgeEvent), timeToRoot);
				
				// Find the appropriate interval and sample coalescence time
				for (NetworkEvent event : networkEventList) {
					if (event.time > checkTime) {
						// Check if we've reached the stop time (reassortment or edge event)
						if (checkTime >= stopTime) {
							break;
						}
						
						// prevEvent.lineages is the number of existing lineages in the network
						// Add our new lineages (coalLines) to get total
						coalRate =  0.5*coalLines * (coalLines - 1) + coalLines * (prevEvent.lineages - negLines);
						
						reverseCoalRate =  0.5*negLines * (negLines - 1) + negLines * (prevEvent.lineages - negLines);

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
						
						logCoaHRInterval -= coalSurv;
						logCoaHRInterval += reverseCoalSurv;
												
						// Check if coalescence happens before next event AND before stop time
						if (coalescenceTime == minNextTime) {
							// add the event contribution
							double coalDens = Math.log(1.0/coalescentDistr.populationFunction.getPopSize(coalescenceTime));
							logHR -= coalDens;
							logCoaHR -= coalDens;
							timeToNextCoalescence = coalescenceTime - currentTime;
							break;
						}
						
						// Move to next interval, but don't go past stop time
						checkTime = Math.min(event.time, stopTime);
					}
					prevEvent = event;
				}
				// If we've gone past all events and haven't reached stop time, coalescence can still happen above the root
				if (stopTime > checkTime && checkTime >= network.getRootEdge().childNode.getHeight() && coalLines > 0) {
					coalRate =  0.5*(coalLines * (coalLines-1));
					if (zeroActiveReassortment)
						coalRate = 0.0;

					
					
					double startTime = Math.max(checkTime, network.getRootEdge().childNode.getHeight());
					
					double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(startTime);
					double transformedTimeToNextCoal = Randomizer.nextExponential(coalRate);
					double coalescenceTime = coalescentDistr.populationFunction.getInverseIntensity(
							transformedTimeToNextCoal + currentTransformedTime);
					
					double nextEvent = Math.min(coalescenceTime,  stopTime);
					
					double integral = coalescentDistr.populationFunction.getIntegral(startTime, nextEvent);
					
					double coalSurv2 = -coalRate * integral;
					logHR -= coalSurv2;
					logCoaHRInterval -= coalSurv2;
					
					

					coalAboveRoot = true;

					// Reverse survival: negLines lineages among themselves
					if (negLines > 0) {
						throw new IllegalArgumentException("Negative lineages above root not allowed");
					}

					
					// Only use this if it happens before stop time
					if (coalescenceTime < stopTime) {
						timeToNextCoalescence = coalescenceTime - currentTime;
						double coalDens2 = Math.log(1.0/coalescentDistr.populationFunction.getPopSize(coalescenceTime));
						logHR -= coalDens2;
						logCoaHR -= coalDens2;
					}
				}
			}
			
			
			
			// Determine which event happens first
			double treeEventDiff = Math.min(timeToNextEdgeEvent, timeToRoot) - currentTime;
			double timeUntilNextEvent = Math.min(Math.min(timeToNextCoalescence, timeToNextReassortment),
					treeEventDiff);
			
			double timeUntilNextEdgeEvent = timeToNextEdgeEvent - currentTime;
			
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
					logReaHRInterval -= reassortSurvFwd;
					logReaHRInterval += reassortSurvRev;
				} else {
					reassortSurvFwd = -totalReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue() * timeUntilNextEvent;
					reassortSurvRev = -totalReverseReassortmentProb * coalescentDistr.reassortmentRateInput.get().getArrayValue() * timeUntilNextEvent;
					logHR -= reassortSurvFwd;
					logHR += reassortSurvRev;
				}
			}
			
			if (activeEdges.size()>1|| inactiveEdges.size()>1)
				exitOnError=true;

			
			// Determine what type of event occurred and handle it
			if (timeUntilNextEvent == timeToNextCoalescence) {
				// add the event contribution
				double coalEventHR = handleCoalescenceEvent(coalLines, activeEdges, excludedEdges, segsToAddList, inactiveEdges, segsToRemoveList, currentTime, rootEdge, edgesAdded, coalAboveRoot);
				logHR += coalEventHR;
				if(exitafternext) {
					if (activeEdges.size()==1) {
						// put all segments on the root edge
						activeEdges.get(0).hasSegments.or(segsToAddList.get(0));
						return logHR;
					}
				}

				if (currentTime==Double.POSITIVE_INFINITY) {
					
					throw new IllegalArgumentException("Current time infinite");
				}
			} else if (timeUntilNextEvent == timeToNextReassortment) {
				// add the event contribution
				double reassortDens = 0.0;
				logHR -= reassortDens;
				logReaHR -= reassortDens;
					
				double reassortEventHR = handleReassortmentEvent(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, reassortmentObsProb, totalReassortmentProb, currentTime, edgesAdded);
				logHR += reassortEventHR;
			} else if (timeUntilNextEvent == timeUntilNextEdgeEvent) {
				double edgeTraversalHR = handleEdgeTraversalEvent(nextInactiveEventIndex, nextEventIndex, inactiveEdges, segsToRemoveList, activeEdges, 
						segsToAddList, reverseCoalRate, totalReverseReassortmentProb, reassortmentObsProbInactive, negLines, excludedEdgesReverse);
				logHR += edgeTraversalHR;
			}
			
			// add the root edge as an active edge if we have reached the root
			if (timeToRoot < Double.POSITIVE_INFINITY && currentTime == timeToRoot && !activeEdges.isEmpty()) {

				if (!activeEdges.contains(rootEdge)) {
					activeEdges.add(rootEdge);
					segsToAddList.add((BitSet) rootEdge.hasSegments.clone());
					rootEdge.hasSegments.clear();
					if (inactiveEdges.contains(rootEdge)) {
						segsToAddList.get(activeEdges.indexOf(rootEdge)).andNot(segsToRemoveList.get(inactiveEdges.indexOf(rootEdge)));
						segsToRemoveList.remove(inactiveEdges.indexOf(rootEdge));
						inactiveEdges.remove(rootEdge);
					}
				}else {
					int idx = activeEdges.indexOf(rootEdge);
					segsToAddList.get(idx).or(rootEdge.hasSegments);
					rootEdge.hasSegments.clear();
					throw new IllegalArgumentException("Root edge already active when reaching root");
				}
				exitafternext = true;
				timeToRoot = Double.POSITIVE_INFINITY;
			}
			round++;
			if (round>100000) {
				System.err.println("Too many rounds in divertSegmentsToAncestors stuck at time " + currentTime + " " + network.getRootEdge().childNode.getHeight());
				System.err.println("likely happens because the reassortment rate is >> than 1/Ne, meaning that the network resampled would go on to infinite time, error message here for debugging");				
				return Double.NEGATIVE_INFINITY;
			}
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
	@Override
	protected double handleCoalescenceEvent(int coalLines, List<NetworkEdge> activeEdges, List<NetworkEdge> excludedEdges, List<BitSet> segsToAddList,
			List<NetworkEdge> inactiveEdges, List<BitSet> segsToRemoveList,
			double currentTime, NetworkEdge rootEdge, List<NetworkEdge> edgesAdded, boolean coalAboveRoot) {
		double logHR = 0.0;
		
		// Step 1: sample which of the coalLines will coalesce
		int coalLineIdx1 = Randomizer.nextInt(coalLines);
		// Forward proposal probability: choosing which coal line
		
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

		if (coexistingLineages.contains(edge1)) {
			throw new IllegalArgumentException("Error: coexisting lineages contains edge1: this should not happen");
		}
		
		
		double[] sampleProb = new double[coexistingLineages.size()];
		
		
		for (int i = 0; i < coexistingLineages.size(); i++) {
			if (activeEdges.contains(coexistingLineages.get(i)) && coexistingLineages.get(i).hasSegments.isEmpty())
				sampleProb[i] = 0.5;
			else
				sampleProb[i] = 1;
		}
		
		double sumProb = 0.0;
		for (double p : sampleProb) {
			sumProb += p;
		}
		

		double cumsum = 0.0;
		double rand = Randomizer.nextDouble() * sumProb;
		int coalLineIdx2 = -1;
		for (int i = 0; i < sampleProb.length; i++) {
			cumsum += sampleProb[i];
			if (rand <= cumsum) {
				coalLineIdx2 = i;
				break;
			}
		}
		if (coalAboveRoot)
			coalLineIdx2 = Randomizer.nextInt(coexistingLineages.size());

		// Step 3: create a new coalescent node and add it to the network
		NetworkNode coalescentNode = new NetworkNode();
		coalescentNode.setHeight(currentTime);
				
		NetworkEdge edge2 = coexistingLineages.get(coalLineIdx2);
		
		int iidx2 = inactiveEdges.indexOf(edge2);
		
		boolean edge2IsRoot = edge2==network.getRootEdge() ? true : false;
		boolean edge1IsRoot = edge1==network.getRootEdge() ? true : false;
		
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
			
			
			if (iidx2 != -1) {
				edge2.hasSegments.andNot(segsToRemoveList.get(iidx2));				
				segsToRemoveList.get(iidx2).andNot(edge1.hasSegments);				
				if (!segsToRemoveList.get(iidx2).isEmpty()) {
					inactiveEdges.set(iidx2, parentEdge);
				}
			}
			
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
			}

			segsToAddList.get(idx).andNot(edge2.hasSegments);

			if (segsToAddList.get(idx).isEmpty()) {
				activeEdges.remove(idx);
				segsToAddList.remove(idx);
			}else {
				activeEdges.set(idx, parentEdge);
			}

			networkEdges.add(parentEdge);
			if (iidx2 != -1) {
				edge2.hasSegments.andNot(segsToRemoveList.get(iidx2));				
				segsToRemoveList.get(iidx2).andNot(edge1.hasSegments);				
				if (!segsToRemoveList.get(iidx2).isEmpty()) {
					inactiveEdges.set(iidx2, parentEdge);
				}
			}
		}
		
		if (edge1IsRoot || edge2IsRoot) {
			network.setRootEdge(coalescentNode.getParentEdges().get(0));
		}
				
		if (iidx2 != -1 && segsToRemoveList.get(iidx2).isEmpty()) {
			inactiveEdges.remove(iidx2);
			segsToRemoveList.remove(iidx2);
		}
		
		return logHR;
	}

	@Override
	protected double handleReassortmentEvent(List<NetworkEdge> activeEdges, List<NetworkEdge> inactiveEdges,
			List<BitSet> segsToAddList, List<BitSet> segsToRemoveList,
			double[] reassortmentObsProb, double totalReassortmentProb, double currentTime, List<NetworkEdge> edgesAdded) {
		double logHR = 0.0;

		
		int reassortEdgeIdx = Randomizer.nextInt(activeEdges.size());

		// Sample which segments will reassort (each with 0.5 probability, conditional on observable reassortment)
		BitSet segsToReassort = segsToAddList.get(reassortEdgeIdx);
		BitSet segsToStay = new BitSet();
		BitSet segsToGo = new BitSet();
		
		segsToStay.clear();
		segsToGo.clear();
		for (int segIdx = segsToReassort.nextSetBit(0); segIdx != -1; segIdx = segsToReassort.nextSetBit(segIdx + 1)) {
			if (Randomizer.nextBoolean()) {
				segsToStay.set(segIdx);
			} else {
				segsToGo.set(segIdx);
			}
		}

		if (segsToGo.isEmpty())
			return 0.0;
		
		
		NetworkEdge parentToStay = new NetworkEdge();
		parentToStay.hasSegments = (BitSet) activeEdges.get(reassortEdgeIdx).hasSegments.clone();
		
		if (parentToStay.hasSegments.isEmpty() && segsToStay.isEmpty())
			return 0.0;
				
		int segsLeft = activeEdges.get(reassortEdgeIdx).hasSegments.cardinality();
		int idx =-1;
		if (inactiveEdges.contains(activeEdges.get(reassortEdgeIdx))){
			idx = inactiveEdges.indexOf(activeEdges.get(reassortEdgeIdx));
			segsLeft -= segsToRemoveList.get(idx).cardinality();
		}

		// prob that previously, the event was unobserved
		double probNoReassortVisibleBefore = 2*Math.pow(0.5, activeEdges.get(reassortEdgeIdx).hasSegments.cardinality());

		if (Randomizer.nextDouble() > probNoReassortVisibleBefore)
			return 0.0;
		
		if (idx!=-1){
			activeEdges.get(reassortEdgeIdx).hasSegments.andNot(segsToRemoveList.get(idx));
			inactiveEdges.set(idx, parentToStay);
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
//		System.out.println(activeEdges.get(reassortEdgeIdx).hasSegments);

		activeEdges.get(reassortEdgeIdx).hasSegments.or(segsToReassort);

		
		
		double reassortDens = 0.0;
		
		boolean unobservableEvent = false;
		
		if (segsToStay.isEmpty() && segsLeft==0) {
//			System.out.println(parentToStay.hasSegments + "dsaf" + segsToStay + " " + segsLeft + " " + segsToGo + " " + activeEdges.get(reassortEdgeIdx).hasSegments);
		}else {
			logHR -= (activeEdges.get(reassortEdgeIdx).hasSegments.cardinality() -1)*Math.log(0.5); // Reverse proposal: choosing which segments go where (conditioned on at least one on each side)
			if (coalescentDistr.timeVaryingReassortmentRates != null) {
				reassortDens = Math.log(1*coalescentDistr.timeVaryingReassortmentRates.getPopSize(currentTime));
			}else {
				reassortDens = Math.log(1*coalescentDistr.reassortmentRateInput.get().getArrayValue());
			}
		}// else, the event would not actually be observed, and we ignore its contribution to the HR
		
		logHR -= reassortDens;

		boolean wasRoot = activeEdges.get(reassortEdgeIdx)==network.getRootEdge() ? true : false;
		
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
		
		if (wasRoot) {
			network.setRootEdge(parentToStay);
		}
		
		logReaHR += logHR;
		
		return logHR;
	}

	@Override
	protected double handleEdgeTraversalEvent(int nextInactiveEventIndex, int nextEventIndex,
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

	@Override
	protected double handleInactiveEdgeTraversal(int idx, List<NetworkEdge> inactiveEdges, 
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

	@Override
	protected double handleActiveEdgeTraversal(int idx, List<NetworkEdge> activeEdges, 
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

	@Override
	protected double handleReassortment(List<NetworkEdge> activeEdges, List<BitSet> segsToAddList,
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
		
		// now, check what to remove
		BitSet segsToRemoveLeft = new BitSet();
		BitSet segsToRemoveRight = new BitSet();
		int segs = edge.hasSegments.cardinality();
		
		boolean createdEmptyEdge = false;

		// remove from inactive edges
		if (iidx != -1) {
			int segsCard = segsToRemoveList.get(iidx).cardinality();
			segsToRemoveLeft.or(segsToRemoveList.get(iidx));
			segsToRemoveRight.or(segsToRemoveList.get(iidx));
			
			segsToRemoveLeft.and(edge.parentNode.getParentEdges().get(0).hasSegments);
			segsToRemoveRight.and(edge.parentNode.getParentEdges().get(1).hasSegments);
			// Reverse proposal probability: For each segment, choose a parent at random
			if (inactiveEdges.get(iidx).hasSegments.cardinality() == segsCard && segsToGoLeft.isEmpty() && segsToGoRight.isEmpty()){// move creates two empty parents
				createdEmptyEdge = true;
			}else if (segsToRemoveLeft.cardinality()==edge.parentNode.getParentEdges().get(0).hasSegments.cardinality() && segsToGoLeft.isEmpty()){ // move creates one empty parent
				createdEmptyEdge = true;
			}else if (segsToRemoveRight.cardinality()==edge.parentNode.getParentEdges().get(1).hasSegments.cardinality() && segsToGoRight.isEmpty()){  // move creates one empty parent
				createdEmptyEdge = true;
			}
			// the totalReverseReassortmentProb condition is to ignore the first event
			if (totalReverseReassortmentProb>0.0) {
				if (createdEmptyEdge) {
					
					double logBinomval = (segs-1)*Math.log(0.5);
					logHR += logBinomval;

					// Reverse proposal probability: choosing which edge					
					if (coalescentDistr.timeVaryingReassortmentRates != null) {
						double popSize = coalescentDistr.timeVaryingReassortmentRates.getPopSize(edge.parentNode.getHeight());
						logHR += Math.log(popSize);
					} else {
						logHR += Math.log(coalescentDistr.reassortmentRateInput.get().getArrayValue());
					}			
				}else {
					logHR += Math.log(0.5) * segsCard;
				}
			}
			
			edge.hasSegments.andNot(segsToRemoveList.get(iidx));
			inactiveEdges.remove(iidx);
			segsToRemoveList.remove(iidx);
		}
		
		// randomly assign each segment to go left or right
		if (idx != -1) {
			if (createdEmptyEdge) {
				// follow the non empty parent with all segments
				if (segsToRemoveLeft.cardinality()==edge.parentNode.getParentEdges().get(0).hasSegments.cardinality()){
					segsToGoRight.or(segsToAddList.get(idx));
				}else {
					segsToGoLeft.or(segsToAddList.get(idx));
				}

			}else {
				for (int segIdx = segsToAddList.get(idx).nextSetBit(0); segIdx != -1; segIdx = segsToAddList.get(idx).nextSetBit(segIdx + 1)) {
					if (Randomizer.nextBoolean()) {
						segsToGoLeft.set(segIdx);
					} else {
						segsToGoRight.set(segIdx);
					}
				}
				// Forward proposal probability: For each segment, choose a parent at random
				logHR -= Math.log(0.5) * segsToAddList.get(idx).cardinality();
			}
			
			edge.hasSegments.or(segsToAddList.get(idx));
			activeEdges.remove(idx);
			segsToAddList.remove(idx);

		}

		
		
		if (createdEmptyEdge) {
			// check if either parent is now empty
			boolean oneIsEmpty = false;
			BitSet leftParent = (BitSet) edge.parentNode.getParentEdges().get(0).hasSegments.clone();
			leftParent.or(segsToGoLeft);
			leftParent.andNot(segsToRemoveLeft);
			
			BitSet rightParent = (BitSet) edge.parentNode.getParentEdges().get(1).hasSegments.clone();
			rightParent.or(segsToGoRight);
			rightParent.andNot(segsToRemoveRight);
			
			if (!leftParent.isEmpty() && !rightParent.isEmpty()) {
				System.out.println("error : neither parent is empty after creating empty edge");
			}					
		}else {
			boolean oneIsEmpty = false;
			BitSet leftParent = (BitSet) edge.parentNode.getParentEdges().get(0).hasSegments.clone();
			leftParent.or(segsToGoLeft);
			leftParent.andNot(segsToRemoveLeft);
			
			BitSet rightParent = (BitSet) edge.parentNode.getParentEdges().get(1).hasSegments.clone();
			rightParent.or(segsToGoRight);
			rightParent.andNot(segsToRemoveRight);
			if (leftParent.isEmpty() || rightParent.isEmpty()) {
				System.out.println("error : one parent is empty after reassortment traversal without creating empty edge");
			}
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
		
		logReaHR += logHR;
		return logHR;
	}
	
	@Override
	protected double handleCoalescence(List<NetworkEdge> activeEdges, List<BitSet> segsToAddList,
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
			if (iidx1 != -1 && segsToRemoveList.get(iidx1).cardinality()==inactiveEdges.get(iidx1).hasSegments.cardinality() && segsToAdd1.isEmpty()) {
				calculateReverseProb = true;
			}
			if (iidx2 != -1 && segsToRemoveList.get(iidx2).cardinality() == inactiveEdges.get(iidx2).hasSegments.cardinality() && segsToAdd2.isEmpty()) {
				calculateReverseProb = true;
			}
//			System.out.println("rev coal" + calculateReverseProb);
			if (calculateReverseProb) {

				logHR += Math.log(1.0/(coalescentDistr.populationFunction.getPopSize(edge.parentNode.getHeight())));				
//				System.out.println("remove coal at time " + edge.parentNode.getHeight());
//				System.out.println(segsToRemoveList.get(iidx1).cardinality() + " " + segsToRemoveList.get(iidx2).cardinality());
//				System.out.println(segsToAdd1.cardinality() + " " + segsToAdd2.cardinality());
//				System.out.println(edge.hasSegments.cardinality() + " " + sisterEdge.hasSegments.cardinality());
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
//				    logHR += Math.log(1.0 / reverseCoalLines);
//				    int partnerCount = sisterAlsoEmpty ? reverseCoexistingLineages.size() : reverseCoexistingLineages.size() + 1;

//				    int partnerCount = reverseCoexistingLineages.size() + 1;
//				    partnerCount = Math.max(1, partnerCount);
                    
//				    logHR += Math.log(1.0 / partnerCount);
//				    System.out.println("c " + partnerCount);
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
//		System.out.println("After coale traversal: " + edge.parentNode.getHeight() + " " + inactiveEdges.size() + " " + activeEdges.size() + " " + logHR);
		logCoaHR += logHR;
		return logHR;
	}


}
