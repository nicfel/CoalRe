package coalre.simulator;

import java.util.*;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.Binomial;
import beast.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.Network;
import coalre.network.NetworkNode;

public class SimulateCoalescentWithReassortment extends Network implements StateNodeInitialiser {

    final public Input<Double> rRateInput = new Input<>("rRate",
            "reassortmentRate", Validate.REQUIRED);

    final public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    final public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment trees to initialize.", new ArrayList<>());

    final public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set specifying leaf ages.");

    final public Input<TaxonSet> taxonSetInput = new Input<>("taxonset",
            "set of taxa that correspond to the leafs in the network",
            Validate.XOR, traitSetInput);


    private ArrayList<NetworkEdge> extantLineages;
    private ArrayList<NetworkNode> remainingSampleNodes;

    private PopulationFunction populationFunction;

    private int nSegments;

    public void initAndValidate() {
        // get the sampling times in order
        remainingSampleNodes = new ArrayList<>();
        extantLineages = new ArrayList<>();

        TraitSet traitSet = traitSetInput.get();
        if (traitSet != null && traitSet.isDateTrait()) {

            for (String taxonName : traitSet.taxaInput.get().getTaxaNames()) {
                NetworkNode sampleNode = new NetworkNode();
                sampleNode.setTaxonLabel(taxonName);
                sampleNode.setHeight(traitSet.getValue(taxonName));
                remainingSampleNodes.add(sampleNode);
            }
            remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

        } else {
            for (String taxonName : taxonSetInput.get().getTaxaNames()) {
                NetworkNode sampleNode = new NetworkNode();
                sampleNode.setTaxonLabel(taxonName);
                sampleNode.setHeight(0.0);
                remainingSampleNodes.add(sampleNode);
            }
        }

        populationFunction = populationFunctionInput.get();

        nSegments = segmentTreesInput.get().size();

        if (nSegments==0) {
            throw new IllegalArgumentException("Need at least one segment tree!");
        }

        simulateNetwork();

        super.initAndValidate();
    }

    public void simulateNetwork() {
        double currentTime = 0;
        double timeUntilNextSample;
        do {
            // get the timing of the next sampling event
            if (!remainingSampleNodes.isEmpty()) {
                timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
            } else {
                timeUntilNextSample = Double.POSITIVE_INFINITY;
            }

            // get the current propensities
            int k = extantLineages.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = Randomizer.nextExponential(0.5*k*(k-1));
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = Randomizer.nextExponential(k*rRateInput.get());

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal)
                    coalesce(currentTime);
                else
                    reassort(currentTime);
            } else {
                currentTime += timeUntilNextSample;
                sample();
            }

        }
        while (extantLineages.size() > 1 || !remainingSampleNodes.isEmpty());

        setRoot(extantLineages.get(0).getChildNode());
    }

    private void sample() {
        // sample the network node
        NetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        hasSegs.set(0, nSegments-1);
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        extantLineages.add(lineage);

        remainingSampleNodes.remove(0);
    }

    private void coalesce(double coalescentTime) {
        // Sample the pair of lineages that are coalescing:
        NetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        NetworkEdge lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        while (lineage1 == lineage2)
            lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

        // Create coalescent node
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        lineage1.setParentNode(coalescentNode);
        lineage2.setParentNode(coalescentNode);

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineage1.hasSegments);
        hasSegments.or(lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);
    }

    private void reassort(double reassortmentTime) {
        NetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

        BitSet hasSegs_left = new BitSet();
        BitSet hasSegs_right = new BitSet();

        for (int segIdx = lineage.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineage.hasSegments.nextSetBit(segIdx+1)) {
            if (Randomizer.nextBoolean()) {
                hasSegs_left.set(segIdx);
            } else {
                hasSegs_right.set(segIdx);
            }
        }

        // Stop here if reassortment event is unobservable
        if (hasSegs_left.cardinality() == 0 || hasSegs_right.cardinality() == 0)
            return;

        // Create reassortment node

        NetworkNode node = new NetworkNode();
        node.setHeight(reassortmentTime)
                .addChildEdge(lineage);
        lineage.setParentNode(node);

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);
    }

    @Override
    public void initStateNodes() {
        initAndValidate();
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.addAll(segmentTreesInput.get());
    }


}
