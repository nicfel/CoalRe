package coalre;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.simulator.SimulatedCoalescentNetwork;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

public abstract class CoalReTestClass {

    protected TaxonSet getTaxonSet(int nTaxa) {

        List<Taxon> taxa = new ArrayList<>();

        for (int i=0; i<nTaxa; i++) {
            taxa.add(new Taxon("t" + i));
        }

        return new TaxonSet(taxa);
    }


    protected TraitSet getDateTraitSet(TaxonSet taxonSet, Supplier<Double> timeFunc) {

        TraitSet traitSet = new TraitSet();

        StringBuilder traitSetValue = new StringBuilder();

        traitSet.initByName("traitname", "date",
                "taxa", taxonSet,
                "value", taxonSet.getTaxaNames().stream()
                        .map(n -> n + "=" + timeFunc.get())
                        .collect(Collectors.joining(",")));

        return traitSet;
    }

    protected TraitSet getSerialDateTraitSet(TaxonSet taxonSet, double timeWindow) {
        return getDateTraitSet(taxonSet, () -> Randomizer.nextDouble()*timeWindow);
    }

    protected TraitSet getContempDateTraitSet(TaxonSet taxonSet) {
        return getDateTraitSet(taxonSet, () -> 0.0);
    }

    protected List<Tree> getSegmentTreeObjects(int nSegments, TraitSet traitSet) {

        List<Tree> segmentTrees = new ArrayList<>();

        for (int seg=0; seg<nSegments; seg++) {
            Tree tree = new Tree();
            tree.initByName("trait", traitSet, "taxonset", traitSet.taxaInput.get());
            tree.setID("tree"+seg);
            segmentTrees.add(tree);
        }

        return segmentTrees;
    }

    protected Network getContempNetwork(int nTaxa, int nSegments, double reassortmentRate) {

        TaxonSet taxonSet = getTaxonSet(nTaxa);
        TraitSet traitSet = getContempDateTraitSet(taxonSet);
        List<Tree> segmentTrees = getSegmentTreeObjects(nSegments, traitSet);

        return getContempNetwork(segmentTrees, reassortmentRate, traitSet);
    }

    protected Network getContempNetwork(List<Tree> segmentTrees, double reassortmentRate, TraitSet traitSet) {
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        SimulatedCoalescentNetwork network = new SimulatedCoalescentNetwork();
        network.initByName(
                "segmentTree", segmentTrees,
                "traitSet", traitSet,
                "populationModel", popFunc,
                "reassortmentRate", new RealParameter(String.valueOf(reassortmentRate)),
                "enableSegmentTreeUpdate", false);

        return network;
    }

}
