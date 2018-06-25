package coalre.simulator;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

public class SimulatedCoalescentNetwork extends Network implements StateNodeInitialiser {

    final public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRate",
            "Rate of reassortment (per lineage per unit time)", Validate.REQUIRED);

    final public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    final public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "One or more segment trees to initialize.", new ArrayList<>());

    private ArrayList<NetworkEdge> extantLineages;
    private ArrayList<NetworkNode> remainingSampleNodes;

    private PopulationFunction populationFunction;
    private RealParameter reassortmentRate;

    private int nSegments;

    public void initAndValidate() {
        remainingSampleNodes = new ArrayList<>();
        extantLineages = new ArrayList<>();

        nSegments = segmentTreesInput.get().size();

        if (nSegments==0) {
            throw new IllegalArgumentException("Need at least one segment tree!");
        }

        TraitSet traitSet = segmentTreesInput.get().get(0).getDateTrait();
        if (traitSet != null) {

            for (String taxonName : traitSet.taxaInput.get().getTaxaNames()) {
                NetworkNode sampleNode = new NetworkNode();
                sampleNode.setLabel(taxonName);
                sampleNode.setHeight(traitSet.getValue(taxonName));
                remainingSampleNodes.add(sampleNode);
            }
            remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

        } else {
            TaxonSet taxonSet = segmentTreesInput.get().get(0).getTaxonset();

            if (taxonSet == null)
                throw new IllegalArgumentException("Segment trees must define" +
                        " either a trait set or a taxon set.");

            for (String taxonName : taxonSet.getTaxaNames()) {
                NetworkNode sampleNode = new NetworkNode();
                sampleNode.setLabel(taxonName);
                sampleNode.setHeight(0.0);
                remainingSampleNodes.add(sampleNode);
            }
        }

        populationFunction = populationFunctionInput.get();

        reassortmentRate = reassortmentRateInput.get();

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
            double transformedTimeToNextCoal = k>0 ? Randomizer.nextExponential(0.5*k*(k-1)) : 0.0;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = k>0 ? Randomizer.nextExponential(k*reassortmentRate.getValue()) : 0.0;

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

        setRootEdge(extantLineages.get(0));
    }

    private void sample() {
        // sample the network node
        NetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        hasSegs.set(0, nSegments);
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        extantLineages.add(lineage);
        n.addParentEdge(lineage);

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
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineage1.hasSegments);
        hasSegments.or(lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

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
        lineage.parentNode = node;

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

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
