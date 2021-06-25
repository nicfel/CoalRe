package coalre.util;

import beast.app.beauti.BeautiDoc;
import beast.core.*;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.UpDownOperator;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import coalre.distribution.CoalescentWithReassortment;
import coalre.network.SegmentTreeInitializer;
import coalre.operators.NetworkScaleOperator;
import coalre.simulator.SimulatedCoalescentNetwork;

import java.util.*;

/**
 * Class containing a static method used in the BEAUti template to tidy
 * up some loose ends.
 */
public class BEAUtiConnector {

    public static boolean customConnector(BeautiDoc doc) {

        int segTreeCount = 0;

        TraitSet traitSet = null;

        MCMC mcmc = (MCMC) doc.mcmc.get();

        Set<RealParameter> parametersToScaleUp = new HashSet<>();
        Set<RealParameter> parametersToScaleDown = new HashSet<>();

        for (BEASTInterface p : doc.getPartitions("tree")) {
            String pId = BeautiDoc.parsePartition(p.getID());

            System.out.println(pId);

            DummyTreeDistribution dummy = (DummyTreeDistribution)doc.pluginmap.get("CoalescentWithReassortmentDummy.t:" + pId);

            if (dummy == null || !dummy.getOutputs().contains(doc.pluginmap.get("prior")))
                continue;

            segTreeCount += 1;

            GenericTreeLikelihood likelihood = (GenericTreeLikelihood) p;
            Tree segmentTree = (Tree) likelihood.treeInput.get();


            // Tell each segment tree initializer which segment index it's
            // initializing.  (Better way to do this?)
            SegmentTreeInitializer segmentTreeInitializer =
                    (SegmentTreeInitializer) doc.pluginmap.get("segmentTreeInitializerCwR.t:" + pId);
            segmentTreeInitializer.segmentIndexInput.setValue(segTreeCount-1, segmentTreeInitializer);

            // Ensure segment tree initializers come first in list:
            // (This is a hack to ensure that the RandomTree initializers
            // are the ones removed by StateNodeInitializerListInputEditor.customConnector().)
            mcmc.initialisersInput.get().remove(segmentTreeInitializer);
            mcmc.initialisersInput.get().add(0, segmentTreeInitializer);


            // Remove segment trees from standard up/down operators.

            for (Operator operator : mcmc.operatorsInput.get()) {
                if (!(operator instanceof UpDownOperator))
                    continue;

                UpDownOperator upDown = (UpDownOperator) operator;

                boolean segmentTreeScaler = upDown.upInput.get().contains(segmentTree)
                        || upDown.downInput.get().contains(segmentTree);

                if (segmentTreeScaler) {
                    upDown.upInput.get().remove(segmentTree);
                    upDown.downInput.get().remove(segmentTree);
                }

            }

            // Clock rates to add to network up/down:
            for (Operator operator : mcmc.operatorsInput.get()) {
                if (!(operator instanceof UpDownOperator))
                    continue;

                UpDownOperator upDown = (UpDownOperator) operator;

                // Note: built-in up/down operators scale trees _down_ while
                // ours scales trees _up_, hence the up/down reversal.
                for (BEASTObject o : upDown.upInput.get()) {
                    if (o instanceof RealParameter) {
                        if (o.getID().contains("clock"))
                            parametersToScaleDown.add((RealParameter)o);
                    }
                }

                for (BEASTObject o : upDown.downInput.get()) {
                    if (o instanceof RealParameter) {
                        if (o.getID().contains("clock"))
                            parametersToScaleUp.add((RealParameter)o);
                    }
                }
            }

            // Extract trait set from one of the trees to use for network.

            if (traitSet == null && segmentTree.hasDateTrait())
                traitSet = segmentTree.getDateTrait();
        }


        // Add clock rates to network up/down operator.

        NetworkScaleOperator networkUpDown = (NetworkScaleOperator)doc.pluginmap.get("networkUpDownCwR.alltrees");
        if (networkUpDown != null) {
            List<RealParameter> paramsCurrent;

            paramsCurrent = networkUpDown.upParametersInput.get();
            for (RealParameter param : paramsCurrent) {
                if (param.getID().toLowerCase().contains("clock"))
                    networkUpDown.upParametersInput.get().remove(param);
            }

            paramsCurrent = networkUpDown.downParametersInput.get();
            for (RealParameter param : paramsCurrent) {
                if (param.getID().toLowerCase().contains("clock"))
                    networkUpDown.downParametersInput.get().remove(param);
            }

            networkUpDown.upParametersInput.get().addAll(parametersToScaleUp);
            networkUpDown.downParametersInput.get().addAll(parametersToScaleDown);
        }

        // Update network initializer:

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
