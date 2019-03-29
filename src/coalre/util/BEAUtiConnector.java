package coalre.util;

import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.MCMC;
import beast.core.Operator;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.UpDownOperator;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import coalre.distribution.CoalescentWithReassortment;
import coalre.simulator.SimulatedCoalescentNetwork;

/**
 * Class containing a static method used in the BEAUti template to tidy
 * up some loose ends.
 */
public class BEAUtiConnector {

    public static boolean customConnector(BeautiDoc doc) {

        int segTreeCount = 0;

        TraitSet traitSet = null;

        for (BEASTInterface p : doc.getPartitions("tree")) {
            String pId = BeautiDoc.parsePartition(p.getID());

            System.out.println(pId);

            if (doc.pluginmap.get("CoalescentWithReassortmentDummy.t:" + pId) == null)
                continue;

            segTreeCount += 1;

            GenericTreeLikelihood likelihood = (GenericTreeLikelihood) p;
            Tree segmentTree = (Tree) likelihood.treeInput.get();

           // Remove segment trees from up/down operators.

            MCMC mcmc = (MCMC) doc.mcmc.get();
            for (Operator operator : mcmc.operatorsInput.get()) {
                if (!(operator instanceof UpDownOperator))
                    continue;

                UpDownOperator upDown = (UpDownOperator) operator;

                upDown.upInput.get().remove(segmentTree);
                upDown.downInput.get().remove(segmentTree);
            }

            // Extract trait set from one of the trees to use for network.

            if (traitSet == null && segmentTree.hasDateTrait())
                traitSet = segmentTree.getDateTrait();
        }


        if (doc.pluginmap.containsKey("networkCwR.alltrees")) {
            SimulatedCoalescentNetwork network = (SimulatedCoalescentNetwork) doc.pluginmap.get("networkCwR.alltrees");

            // Update number of segments for initializer.
            network.nSegmentsInput.setValue(segTreeCount, network);

            // Provide trait set from first segment tree to network initializer:
            if (traitSet != null)
                network.traitSetInput.setValue(traitSet, network);
        }



        return false;
    }
}
