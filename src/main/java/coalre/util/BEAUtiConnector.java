package coalre.util;

import beastfx.app.inputeditor.BeautiDoc;
import beast.base.inference.Scalable;
import beast.base.core.BEASTInterface;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
// BEAUti 2.8 builds the spec-typed twins of these, so the legacy classes never match.
import beast.base.spec.evolution.likelihood.GenericTreeLikelihood;
import beast.base.spec.evolution.operator.AdaptableOperatorSampler;
import beast.base.spec.evolution.operator.AdaptableVarianceMultivariateNormalOperator;
import beast.base.spec.evolution.operator.UpDownOperator;
import beast.base.spec.inference.operator.Transform;
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

        Set<Scalable> parametersToScaleUp = new HashSet<>();
        Set<Scalable> parametersToScaleDown = new HashSet<>();
        Set<Tree> segmentTrees = new HashSet<>();

        for (BEASTInterface p : doc.getPartitions("tree")) {
            String pId = BeautiDoc.parsePartition(p.getID());

            System.out.println(pId);

            DummyTreeDistribution dummy = (DummyTreeDistribution)doc.pluginmap.get("CoalescentWithReassortmentDummy.t:" + pId);

            if (dummy == null || !dummy.getOutputs().contains(doc.pluginmap.get("prior")))
                continue;

            segTreeCount += 1;

            GenericTreeLikelihood likelihood = (GenericTreeLikelihood) p;
            Tree segmentTree = (Tree) likelihood.treeInput.get();
            segmentTrees.add(segmentTree);


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

            for (Operator aos : mcmc.operatorsInput.get()) {
            	if (aos instanceof AdaptableOperatorSampler) {           		
            		
                    for (Operator operator : ((AdaptableOperatorSampler) aos).operatorsInput.get()) {
		                if (operator instanceof UpDownOperator) {

			                UpDownOperator upDown = (UpDownOperator) operator;

			                boolean segmentTreeScaler = upDown.upInput.get().contains(segmentTree)
			                        || upDown.downInput.get().contains(segmentTree);
			
			                if (segmentTreeScaler) {
			                    upDown.upInput.get().remove(segmentTree);
			                    upDown.downInput.get().remove(segmentTree);
			                }
		                }else if (operator instanceof AdaptableVarianceMultivariateNormalOperator) {
		                	
		                	AdaptableVarianceMultivariateNormalOperator amvn = (AdaptableVarianceMultivariateNormalOperator) operator;
		                	for (Transform t : amvn.transformationsInput.get())
		                		if (t instanceof Transform.UnivariableTransform)
		                			((Transform.UnivariableTransform) t).functionInput.get().remove(segmentTree);
		                	
		                }
                    }
            	}

            }

            // Clock rates to add to network up/down:
            List<Operator> removeOps  = new ArrayList<>();
            for (Operator aos : mcmc.operatorsInput.get()) {
            	if (aos instanceof AdaptableOperatorSampler) {           		
                    for (Operator operator : ((AdaptableOperatorSampler) aos).operatorsInput.get()) {
		                if (!(operator instanceof UpDownOperator))
		                    continue;

		                UpDownOperator upDown = (UpDownOperator) operator;

		                // Note: built-in up/down operators scale trees _down_ while
		                // ours scales trees _up_, hence the up/down reversal.
		                for (Scalable s : upDown.upInput.get()) {
		                    if (isClockParameter(s))
		                        parametersToScaleDown.add(s);
		                }

		                for (Scalable s : upDown.downInput.get()) {
		                    if (isClockParameter(s))
		                        parametersToScaleUp.add(s);
		                }
                    }
                    if (((AdaptableOperatorSampler) aos).treeInput.get().size()!=0)
                    	removeOps.add(aos);                    
            	}
            }
            
            mcmc.operatorsInput.get().removeAll(removeOps);

            // Extract trait set from one of the trees to use for network.

            if (traitSet == null && segmentTree.hasDateTrait())
                traitSet = segmentTree.getDateTrait();
        }


        // Add clock rates to network up/down operator.

        NetworkScaleOperator networkUpDown = (NetworkScaleOperator)doc.pluginmap.get("networkUpDownCwR.alltrees");
        if (networkUpDown != null) {
            // removeIf, not a for-each with remove(): the list being iterated is the very
            // list the Input holds, so removing during iteration threw ConcurrentModificationException.
            networkUpDown.upParametersInput.get().removeIf(BEAUtiConnector::isClockParameter);
            networkUpDown.downParametersInput.get().removeIf(BEAUtiConnector::isClockParameter);

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
          
            network.segmentTreesInput.get().clear();
            network.segmentTreesInput.get().addAll(segmentTrees);

        }
        
        return false;
    }

    /**
     * Scalable is not a BEASTInterface, so the ID has to be read off the underlying object.
     */
    private static boolean isClockParameter(Scalable s) {
        return s instanceof BEASTInterface b
                && b.getID() != null
                && b.getID().toLowerCase().contains("clock");
    }
}
