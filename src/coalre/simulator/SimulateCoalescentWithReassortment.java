package coalre.simulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.Network;
import coalre.network.NetworkNode;

public class SimulateCoalescentWithReassortment extends Network implements StateNodeInitialiser {

    final public Input<Double> rRateInput = new Input<>("rRate",
            "reassortmentRate", Validate.REQUIRED);
    final public Input<Double> coalRateInput = new Input<>("coalRate",
            "coalescent rate", Validate.REQUIRED);
    final public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "dummy segment tree", new ArrayList<>());

    final public Input<Network> networkInput = new Input<>("network",
            "dummy network", Validate.REQUIRED);

    final public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set specifying leaf ages.");

    final public Input<TaxonSet> taxonSetInput = new Input<>("taxonset",
            "set of taxa that correspond to the leafs in the network",
            Validate.XOR, traitSetInput);


    ArrayList<Double> samplingTimes;
    ArrayList<NetworkEdge> extantLineages;
    int highestNetworkNodeNr;
    int reassortmentNumber;
    int[] highestNodeNr;
    TaxonSet taxonset;

    int nSegments;

    public void initAndValidate() {
        // get the sampling times in order
        samplingTimes = new ArrayList<>();
        extantLineages = new ArrayList<>();

        TraitSet traitSet = traitSetInput.get();
        if (traitSet != null && traitSet.isDateTrait()) {
            for (String taxonName : traitSet.taxaInput.get().getTaxaNames()) {
                samplingTimes.add(traitSet.getValue(taxonName));
            }
        } else {
            for (int i=0; i<taxonSetInput.get().getTaxonCount(); i++)
                samplingTimes.add(0.0);
        }

        nSegments = segmentTreesInput.get().size();

        if (nSegments==0) {
            throw new IllegalArgumentException("Need at least one segment tree!");
        }

        simulateNetwork();

        super.initAndValidate();
    }

    public void simulateNetwork() {
        int samplingInterval = 0;
        highestNetworkNodeNr = 0;
        reassortmentNumber = 0;
        double currTime = 0;
        double nextSamplingTime;
        do {
            // get the timing of the next sampling event
            if (samplingInterval < samplingTimes.size()) {
                nextSamplingTime = Collections.min(samplingTimes) - currTime;
            } else {
                nextSamplingTime = Double.POSITIVE_INFINITY;
            }

            // get the current propensities
            int k = extantLineages.size();
            double coalProp = k * (k - 1) / 2 * coalRateInput.get();
            double reasProp = k * rRateInput.get();
            double totalProp = coalProp + reasProp;

            // next event time
            double nextNonSamplingEvent = Randomizer.nextExponential(totalProp);
            if (nextNonSamplingEvent < nextSamplingTime) {
                currTime += nextNonSamplingEvent;
                if (Randomizer.nextDouble() * totalProp < coalProp) {
//        			System.out.println("c");
                    coalesce(currTime);
                } else {
//        			System.out.println("r");
//        			System.out.println(coalProb/(coalProb+reasProb));
                    reassort(currTime);
                }
            } else {
//    			System.out.println("s");
                currTime += nextSamplingTime;
                sample(currTime, samplingInterval);
                samplingInterval++;
            }

        }
        while (extantLineages.size() > 1 || nextSamplingTime < Double.POSITIVE_INFINITY);

        if (extantLineages.size() != 1) {
            throw new IllegalArgumentException("simulation turned awkward");
        }

        setRoot(extantLineages.get(0).node);
    }

    private void sample(double samplingTime, int samplingInterval) {
        HashMap<Integer, Node> segs = new HashMap<>();

        double minSampleTime = Collections.min(samplingTimes);
        int minIndex = -1;

        // check for the index of the minSampleTime
        for (int j = 0; j < samplingTimes.size(); j++) {
            if (samplingTimes.get(j) == minSampleTime) {
                minIndex = j;
                break;
            }
        }

        Boolean[] hasSegs = new Boolean[nSegments];

        for (int i = 0; i < nSegments; i++) {
            hasSegs[i] = true;
        }

        // sample the network node
        NetworkNode n = new NetworkNode();
        n.setTaxonLabel(taxonset.getTaxonId(minIndex));
        n.setHeight(minSampleTime);

        n.setHasSegments(hasSegs);

        samplingTimes.set(minIndex, Double.POSITIVE_INFINITY);
        extantLineages.add(n);
    }

    private void coalesce(double coalescentTime) {
        // sample the first "bunch" of lineages to coalesce
        int lineages1 = Randomizer.nextInt(activeLineages.size());
        int lineages2 = Randomizer.nextInt(activeLineages.size());
        while (lineages1 == lineages2)
            lineages2 = Randomizer.nextInt(activeLineages.size());

        HashMap<Integer, Node> h1 = new HashMap<>(activeLineages.get(lineages1));
        HashMap<Integer, Node> h2 = new HashMap<>(activeLineages.get(lineages2));

        // check if there is any overlap between the bunch of lineages
        Integer[] segs1 = new Integer[h1.size()];
        Integer[] segs2 = new Integer[h2.size()];
        h1.keySet().toArray(segs1);
        h2.keySet().toArray(segs2);


        Set<Integer> segs11 = h1.keySet();
        Set<Integer> segs22 = h2.keySet();

        // get all segments in 1 that are in 2
        segs11.retainAll(segs22);

        // remove the daugher lineages from being active lineages

        HashMap<Integer, Node> parents = new HashMap<>();

        Integer[] overlapSegs = new Integer[segs11.size()];
        segs11.toArray(overlapSegs);

        Boolean[] hasSegs = new Boolean[segmentsTreeInput.get().size()];


        // coalesce overlapping segments
        for (int i = 0; i < overlapSegs.length; i++) {
            Node p = new Node();
            p.setHeight(coalescentTime);
            p.setNr(highestNodeNr[overlapSegs[i]] + samplingTimes.size());
//    		p.setParent(null);
            p.setLeft(h1.get(overlapSegs[i]));
            p.setRight(h2.get(overlapSegs[i]));
            h1.get(overlapSegs[i]).setParent(p);
            h2.get(overlapSegs[i]).setParent(p);

            highestNodeNr[overlapSegs[i]]++;
//    		p.getLength()

            hasSegs[overlapSegs[i]] = true;

            parents.put(overlapSegs[i], p);
        }


//    	System.out.println("....");    	
        for (Integer aSegs1 : segs1) {
            // check if already added
            if (parents.get(aSegs1) == null) {
                parents.put(aSegs1, activeLineages.get(lineages1).get(aSegs1));
            }
        }

        for (Integer aSegs2 : segs2) {
            // check if already added
            if (parents.get(aSegs2) == null)
                parents.put(aSegs2, activeLineages.get(lineages2).get(aSegs2));
        }

        NetworkNode p = new NetworkNode();
        p.setHeight(coalescentTime);
        p.setNr(highestNetworkNodeNr + samplingTimes.size());

        p.setLeft(extantLineages.get(lineages1));
        p.setRight(extantLineages.get(lineages2));

        p.setHasSegments(hasSegs);

        extantLineages.get(lineages1).setParent(p);
        extantLineages.get(lineages2).setParent(p);

        highestNetworkNodeNr++;


        activeLineages.remove(Math.max(lineages1, lineages2));
        activeLineages.remove(Math.min(lineages1, lineages2));

        extantLineages.remove(Math.max(lineages1, lineages2));
        extantLineages.remove(Math.min(lineages1, lineages2));

        activeLineages.add(parents);


        extantLineages.add(p);

    }

    private void reassort(double reassortmentTime) {
        // get which part is reassorting
        int lineages = Randomizer.nextInt(activeLineages.size());


        // get the nodes
        HashMap<Integer, Node> h = activeLineages.get(lineages);

        HashMap<Integer, Node> goesLeft = new HashMap<>();
        HashMap<Integer, Node> goesRight = new HashMap<>();

        Boolean[] hasSegs_left = new Boolean[segmentsTreeInput.get().size()];
        Boolean[] hasSegs_right = new Boolean[segmentsTreeInput.get().size()];

        //keys
        Integer[] keys = new Integer[h.size()];
        h.keySet().toArray(keys);

        for (Integer key : keys) {
            if (Randomizer.nextBoolean()) {
                goesLeft.put(key, h.get(key));
                hasSegs_left[key] = true;
            } else {
                goesRight.put(key, h.get(key));
                hasSegs_right[key] = true;
            }
        }

        // check if the reassortment event is actually observable
        if (goesLeft.size() > 0 && goesRight.size() > 0) {
            // remove old
            activeLineages.remove(lineages);

            if (goesLeft.size() > 0)
                activeLineages.add(goesLeft);
            if (goesRight.size() > 0)
                activeLineages.add(goesRight);

            // the ones that go left have an actual parent
            NetworkNode p_left = new NetworkNode();
            NetworkNode p_right = new NetworkNode();

            p_left.setHeight(reassortmentTime);
            p_right.setHeight(reassortmentTime);

            p_left.setNr(highestNetworkNodeNr + samplingTimes.size());
            highestNetworkNodeNr++;
            p_right.setNr(highestNetworkNodeNr + samplingTimes.size());
            p_left.setReassortmentNumber(reassortmentNumber);
            p_right.setReassortmentNumber(reassortmentNumber);

            reassortmentNumber++;
            highestNetworkNodeNr++;

            p_left.setLeft(extantLineages.get(lineages));
//    		p_right.setLeft(extantLineages.get(lineages));

            p_left.setHasSegments(hasSegs_left);
            p_right.setHasSegments(hasSegs_right);

            extantLineages.get(lineages).setParent(p_left);
            extantLineages.get(lineages).setSecondParent(p_right);

            extantLineages.remove(lineages);
            extantLineages.add(p_left);
            extantLineages.add(p_right);
        }
    }

    @Override
    public void initStateNodes() {
        // get the sampling times in order
        samplingTimes = new ArrayList<>();
        activeLineages = new ArrayList<>();

        taxonset = segmentsTreeInput.get().get(0).getTaxonset();
        TraitSet datetrait = segmentsTreeInput.get().get(0).getDateTrait();
        // get the sampling times of the taxa
        for (int i = 0; i < taxonset.getNrTaxa(); i++)
            samplingTimes.add(datetrait.getValue(taxonset.getTaxonId(i)));


        highestNodeNr = new int[segmentsTreeInput.get().size()];
        for (int i = 0; i < highestNodeNr.length; i++)
            highestNodeNr[i] = 0;

        simulateNetwork();

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.addAll(segmentsTreeInput.get());
    }


}
