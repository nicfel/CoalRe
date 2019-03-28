package coalre.util;

import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.MCMC;
import beast.core.Operator;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.UpDownOperator;
import beast.evolution.tree.Tree;
import coalre.distribution.CoalescentWithReassortment;
import coalre.simulator.SimulatedCoalescentNetwork;

public class UpDownRemover {
    public static boolean customConnector(BeautiDoc doc) {

        int segTreeCount = 0;

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
        }

        // Update number of segments for initializer.

        if (doc.pluginmap.containsKey("networkCwR")) {
            SimulatedCoalescentNetwork network = (SimulatedCoalescentNetwork) doc.pluginmap.get("networkCwR");
            network.nSegmentsInput.setValue(segTreeCount, network);
        }

        return false;
    }
}
