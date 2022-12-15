package coalre.simulator;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

public class SimulatedCoalescentNetwork extends Network {

    public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRate",
            "Rate of reassortment (per lineage per unit time)", Validate.REQUIRED);

    public Input<Function> binomialProbInput = new Input<>("binomialProb",
            "Probability parameter in binomial reassortment distribution.");

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "One or more segment trees to initialize.", new ArrayList<>());

    public Input<Integer> nSegmentsInput = new Input<>("nSegments",
            "Number of segments. Used if no segment trees are supplied.");

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set used to assign leaf ages.");

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set used to define leaves");

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
            "If false, segment tree objects won't be updated to agree with simulated " +
                    "network. (Default true.)", true);

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");

    private PopulationFunction populationFunction;
    private RealParameter reassortmentRate;
    private Function binomialProb;
    

    private int nSegments;

    public void initAndValidate() {

        if (segmentTreesInput.get().isEmpty()) {
            nSegments = nSegmentsInput.get();
        }else {
            nSegments = segmentTreesInput.get().size();
            segmentNames = new String[nSegments];
            // initialize names of segments
            baseName = "Tree.t:";
//            baseName = segmentTreesInput.get().get(0).getID();
//            for (int segIdx=1; segIdx<nSegments; segIdx++) {
//        		String newBase = "";
//            	for (int i = 0; i < baseName.length(); i++) {
//    				if (baseName.substring(0,i+1).contentEquals(segmentTreesInput.get().get(segIdx).getID().substring(0, i+1))) {
//    					newBase = baseName.substring(0,i+1);
//    				}
//    				
//            	}
//            	baseName = newBase;
//            }
            
            for (int segIdx=0; segIdx<nSegments; segIdx++) {
            	segmentNames[segIdx] = segmentTreesInput.get().get(segIdx).getID().replace(baseName, "");
            }
            segmentTreesInput.get().clear();
        }

        populationFunction = populationFunctionInput.get();
        reassortmentRate = reassortmentRateInput.get();
        binomialProb = binomialProbInput.get();

        if (nSegments==0) {
            throw new IllegalArgumentException("Need at least one segment!");
        }

        // Set up sample nodes:

        List<NetworkNode> sampleNodes = new ArrayList<>();
    

        TaxonSet taxonSet = null;
        if (traitSetInput.get() != null)
            taxonSet = traitSetInput.get().taxaInput.get();
        else if (taxonSetInput.get() != null)
            taxonSet = taxonSetInput.get();
        else if (!segmentTreesInput.get().isEmpty())
        	taxonSet = segmentTreesInput.get().get(0).getTaxonset();
        else
            throw new IllegalArgumentException("Taxon set must be specified " +
                    "using either taxonSet, traitSet or provided by a segmentTree input.");

        TraitSet traitSet = null;
        if (traitSetInput.get() != null)
            traitSet = traitSetInput.get();
        else if (!segmentTreesInput.get().isEmpty())
            traitSet = segmentTreesInput.get().get(0).getDateTrait();

        for (int taxonIndex=0; taxonIndex<taxonSet.getTaxonCount(); taxonIndex++) {
            String taxonName = taxonSet.getTaxonId(taxonIndex);

            NetworkNode sampleNode = new NetworkNode();
            sampleNode.setTaxonLabel(taxonName);
            sampleNode.setTaxonIndex(taxonIndex);

            if (traitSet != null)
                sampleNode.setHeight(traitSet.getValue(taxonName));
            else if (!segmentTreesInput.get().isEmpty())
                sampleNode.setHeight(segmentTreesInput.get().get(0).getNode(taxonIndex).getHeight());
            else
                sampleNode.setHeight(0.0);

            sampleNodes.add(sampleNode);
        }

        // Perform network simulation:
        simulateNetwork(sampleNodes);

        // Update segment trees:
        if (enableSegTreeUpdateInput.get()) {
            for (int segIdx = 0; segIdx < nSegments; segIdx++) {
                Tree segmentTree = segmentTreesInput.get().get(segIdx);
                updateSegmentTree(segmentTree, segIdx);
                segmentTree.setEverythingDirty(false);
            }
        }

        // Write simulated network to file if requested
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {

                ps.println(toString());

            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to output file '"
                        + fileNameInput.get() + "'.");
            }
        }

        super.initAndValidate();
    }

    private double getBinomialProb() {
        return binomialProb != null
                ? binomialProb.getArrayValue()
                : 0.5;
    }

    /**
     * Simulate network under coalescent with reassortment model.
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateNetwork(List<NetworkNode> sampleNodes) {
    	System.out.println("resample");

        List<NetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);
        List<NetworkEdge> extantLineages = new ArrayList<>();

        remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

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
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*reassortmentRate.getValue()) : Double.POSITIVE_INFINITY;

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal)
                    coalesce(currentTime, extantLineages);
                else
                    reassort(currentTime, extantLineages);
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }

        }
        while (extantLineages.size() > 1 || !remainingSampleNodes.isEmpty());

        setRootEdge(extantLineages.get(0));
    }

    private void sample(List<NetworkNode> remainingSampleNodes, List<NetworkEdge> extantLineages) {
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

    private void coalesce(double coalescentTime, List<NetworkEdge> extantLineages) {
        // Sample the pair of lineages that are coalescing:
        NetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        NetworkEdge lineage2;
        do {
            lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        } while (lineage1 == lineage2);

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

    private void reassort(double reassortmentTime, List<NetworkEdge> extantLineages) {
        NetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

        BitSet hasSegs_left = new BitSet();
        BitSet hasSegs_right = new BitSet();

        for (int segIdx = lineage.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineage.hasSegments.nextSetBit(segIdx+1)) {
            if (Randomizer.nextDouble() < getBinomialProb()) {
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
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);
    }

}
