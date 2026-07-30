package coalre.simulator;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class SuperspreadingSIRwithReassortment extends SIRwithReassortment implements Loggable {

    public Input<Double> kInput = new Input<>("k",
            "k-value for the negative binomial distribution", 0.1);
    
	public Input<Double> recoveredPercentageInput = new Input<>("recoveredPercentage",
			"percentage of recovered individuals at the start of the simulation", 0.0);


    @Override
    public void initAndValidate(){

        List<String> IDs = new ArrayList<>();
        if (segmentTreesInput.get().isEmpty()) {
            nSegments = (int) nSegmentsInput.get();
        }else {
            nSegments = segmentTreesInput.get().size();
            segmentNames = new String[nSegments];
            // initialize names of segments
            baseName = "";
            for (int segIdx=0; segIdx<nSegments; segIdx++) {
                IDs.add(segmentTreesInput.get().get(segIdx).getID());
                segmentNames[segIdx] = segmentTreesInput.get().get(segIdx).getID().replace(baseName, "");
            }
//            segmentTreesInput.get().clear();
        }

        simulateNetwork();
        // Create the taxon sets from the sampled individuals
        List<Taxon> taxonList = new ArrayList<>();
        for (NetworkNode leaf : getLeafNodes()) {
            Taxon taxon = new Taxon(leaf.getTaxonLabel());
            taxonList.add(taxon);
        }
        // make a TaxonSet from the taxa
        TaxonSet taxonset = new TaxonSet();
        taxonset.initByName("taxon", taxonList);

        // create the trait set using the sampling times
        TraitSet traitSet = new TraitSet();
        // make a string of values with leaf.getTaxonLabel()+"=" +leaf.getHeight() and individual entries seprated by ,
        String values = "";
        for (NetworkNode leaf : getLeafNodes()) {
            values += leaf.getTaxonLabel()+"=" +leaf.getHeight() + ",";
        }
        traitSet.initByName("value", values.substring(0, values.length()-1), "taxa", taxonset, "traitname", "date-backward");

        // Update segment trees:
        if (enableSegTreeUpdateInput.get()) {
            for (int segIdx = 0; segIdx < nSegments; segIdx++) {
                // get the segment tree, set the taxonsets and update the segment tree based on the network
                Tree segmentTree = new Tree();
                segmentTree.initByName("taxonset", taxonset);
                updateSegmentTree(segmentTree, segIdx);
                segmentTree.setID(IDs.get(segIdx));
//                segmentTreesInput.get().add( segmentTree);
                segmentTree.setEverythingDirty(false);
                segmentTreesInput.get().set(segIdx, segmentTree);
            }
        }
        // Write simulated network to file if requested
        super.initAndValidate();

        for (int segIdx = 0; segIdx < nSegments; segIdx++) {
            System.out.println(segmentTreesInput.get().get(segIdx).getID());
            System.out.println(segmentTreesInput.get().get(segIdx).getRoot().toNewick());
            System.out.println(segmentTreesInput.get().get(segIdx).getNodeCount());

        }

        // Write simulated network to file if requested
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {

                ps.println(toString().replace("0.0)", "0.0000000001)"));
                for (int segIdx = 0; segIdx < nSegments; segIdx++) {
                    ps.println(segmentTreesInput.get().get(segIdx).getRoot().toNewick() +";");
                }
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to output file '"
                        + fileNameInput.get() + "'.");
            }
        }
        super.initAndValidate();
    }

    @Override
    protected void simulateNetwork() {
        List<Individual> sampledIndividuals = new ArrayList<>();

        // set the starting conditions for simulations
        do {
            sirEvents = new ArrayList<>();

            double R = (int) (recoveredPercentageInput.get() * populationSizeInput.get());
            double S = populationSizeInput.get() - 1 - R;
            double I = 1;
            int individualNr = 0;
            // start the simulations while keeping track of the different individuals
            List<Individual> activeIndividuals = new ArrayList<>();
            sampledIndividuals = new ArrayList<>();

            Individual root = new Individual(individualNr);
            individualNr++;
            activeIndividuals.add(root);
            double time = 0.0;
            do {
                double transmissionRate = I * recoveryRateInput.get(); // This should be recovery rate, as we model the Reff through the avg number of secondary infections below
                double recoveryRate = I * recoveryRateInput.get();
                double waningRate = R * waningImmunityRateInput.get();

                double totalRate = transmissionRate + recoveryRate + waningRate;
                double nextEventTime = Math.log(1.0 / Randomizer.nextDouble()) / totalRate;

                // pick which event happens next based on the rates, by drawing at random
                double random = Randomizer.nextDouble() * totalRate;

                if (random < transmissionRate) {
                    // pick the number of offsprings from a negative binomial distribution with R and k
                    // from a gamma and a poisson distribution
                    double secondary_infections = transmissionRateInput.get()/recoveryRateInput.get();
                    double gamma = Randomizer.nextGamma(kInput.get(), kInput.get() / (double) secondary_infections);
                    int nOffspring = (int) Randomizer.nextPoisson(gamma);
                    int isR = 0;
                    // for each offspring, randomly sample if it hits an R
                    for (int i = 0; i < nOffspring; i++) {
                        if (Randomizer.nextDouble() < (R+1) / (double) (S + I + R)) {
                            // sample the individual
                            isR++;
                        }
                    }
                    nOffspring -= isR;

                    List<Integer> involvedIndividuals = new ArrayList<>();

                    // don't do anything if there are no offspring
                    if (nOffspring > 0){
//                        System.out.println(nOffspring + " " + activeIndividuals.size());
                        // transmission event
                        // pick a random individual to transmit to
                        Individual parent = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));

                        // build the superspreading event as a series of regular infection events
                        for (int i = 0; i < nOffspring; i++) {
                            System.out.println( " time " + time);
                            // set time for parent
                            parent.setTime(time + nextEventTime);
                            // create two new individuals
                            Individual child1 = new Individual(individualNr);
                            individualNr++;
                            Individual child2 = new Individual(individualNr);
                            individualNr++;
                            // start the first new lineage
                            parent.addChild(child1);
                            child1.addParent(parent);
                            // add the child to the second parent
                            parent.addChild(child2);
                            child2.addParent(parent);

                            // add the child to the list of active individuals
                            activeIndividuals.add(child1);

                            involvedIndividuals.add(child1.number);
                            involvedIndividuals.add(child2.number);
                            involvedIndividuals.add(parent.number);


                            // pick if the individual to be infected, assuming an individual is co-infected once at most
                            System.out.println(I + " coinf " + activeIndividuals.size());
                            double probCoInf = S / (double) (S + activeIndividuals.size() - 2 - i);
                            if (Randomizer.nextDouble() > probCoInf) {
                                // offset for the number of active individuals to draw from
                                int totOptions = i==0? 1 : 2 + i;
                                System.out.println(totOptions + " " + activeIndividuals.size());

                                // co-infection event
                                // pick another active individual that is not the parent, while making sure it is not
                                // the same as any of the previous individuals in this event
                                Individual parent2 = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()-totOptions));
                                while (parent2 == parent) {
                                    parent2 = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()-totOptions));
                                }
                                parent2.setTime(time + nextEventTime);

                                // add the child to the second parent
                                child2.setTime(time + nextEventTime);

                                // create a third child
                                Individual child3 = new Individual(individualNr);

//                                System.out.println("other lineage " + parent2.number);

//                                involvedIndividuals.add(parent2.number);
//                                involvedIndividuals.add(child3.number);

                                individualNr++;
                                child3.addParent(child2);
                                child2.addChild(child3);

                                child3.addParent(parent2);
                                parent2.addChild(child3);

                                activeIndividuals.add(child3);
                                activeIndividuals.remove(parent2);
                                sirEvents.add(new SIREvents(0, time + nextEventTime));
                            } else {
                                // update the population counts
                                S--;
                                I++;
                                activeIndividuals.add(child2);
                                sirEvents.add(new SIREvents(1, time + nextEventTime));
                            }
                            // remove the parent
                            activeIndividuals.remove(parent);
                            parent = child1;
                            time += 0.0000000001;
                        }

//                        System.out.println(involvedIndividuals);
                        // check if all the involved individuals are unique
//                        if (involvedIndividuals.size() != involvedIndividuals.stream().distinct().count()) {
//                            System.out.println("not unique");
//                            System.out.println(involvedIndividuals);
//                            System.exit(0);
//                        }

                    }
                } else if (random < (transmissionRate + recoveryRate)) {
                    // pick a random individual to recover
                    Individual individual = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                    individual.setTime(time + nextEventTime);
                    // choose if that individual will be sampled
                    if (Randomizer.nextDouble() < samplingProbabilityInput.get()) {
                        // sample the individual
                        sampledIndividuals.add(individual);
                    }
                    // remove the individual from the list of active individuals
                    activeIndividuals.remove(individual);
                    // update the population counts
                    I--;
                    R++;
                    sirEvents.add(new SIREvents(2, time + nextEventTime));
                } else {
                    System.out.println(transmissionRate + recoveryRate);
                    System.out.println(transmissionRate + " " + recoveryRate + " " + waningRate);
                    System.out.println(I);
                    System.out.println((S + I) );
                    System.exit(0);
                    // waning immunity event
                    // update the population counts
                    R--;
                    S++;
                    sirEvents.add(new SIREvents(3, time + nextEventTime));
                }
                time += nextEventTime;
//                System.out.println(time + " S: " + S + " I: " + I + " R: " + R + " " + activeIndividuals.size() + " " + sampledIndividuals.size());
            } while (time < simulationTimeInput.get() && I > 0);
            System.out.println("start building network from " + sampledIndividuals.size() + " sampled individuals" + " simulation time was " + time);
            // build the network from the sampled individuals
        }while (sampledIndividuals.size()<minSamplesInput.get());
        buildNetwork(sampledIndividuals);
//        System.out.println(lineages);
    }

}

