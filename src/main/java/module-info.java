open module coalre {
    requires beast.pkgmgmt;
    requires beast.base;
    // coalre.util.InitFromTree uses beastfx.app.inputeditor.BeautiDoc
    requires beast.fx;

    requires java.desktop;               // Swing/AWT in coalre.networkannotator
    requires java.xml;                   // org.w3c.dom in Network.getDOM()

    requires org.antlr.antlr4.runtime;   // checked-in Network.g4 parser
    requires colt;                       // network operators and simulators
    requires commons.math3;              // linear algebra in coalre.dynamics

    exports coalre.distribution;
    exports coalre.dynamics;
    exports coalre.network;
    exports coalre.network.parser;
    exports coalre.networkannotator;
    exports coalre.operators;
    exports coalre.simulator;
    exports coalre.statistics;
    exports coalre.util;

    provides beast.base.core.BEASTInterface with
        coalre.distribution.CoalescentWithReassortment,
        coalre.distribution.emptyNetworkEdgesPrior,
        coalre.distribution.NetworkDistribution,
        coalre.distribution.NetworkIntervals,
        coalre.distribution.TipPrior,
        coalre.network.Network,
        coalre.network.SegmentTreeInitializer,
        coalre.operators.AddRemoveReassortment,
        coalre.operators.AddRemoveReassortmentCoalescent,
        coalre.operators.AddRemoveReassortmentExponential,
        coalre.operators.BranchUniform,
        coalre.operators.DivertSegmentOperator,
        coalre.operators.DivertSegmentAndResimulate,
        coalre.operators.RemoveSegmentAndResimulate,
        coalre.operators.AddRemoveAndResimulate,
        coalre.operators.GibbsOperatorAboveSegmentRoots,
        coalre.operators.MultiTipDatesRandomWalker,
        coalre.operators.NetworkExchange,
        coalre.operators.NetworkExchangeAndResimulate,
        coalre.operators.NetworkScaleOperator,
        coalre.operators.SubNetworkSlide,
        coalre.operators.SubNetworkLeap,
        coalre.operators.SubNetworkLeapAndResimulate,
        coalre.operators.NetworkSPR,
        coalre.operators.SubTreeSlideOnNetwork,
        coalre.operators.TipReheight,
        coalre.operators.UniformNetworkNodeHeightOperator,
        coalre.operators.UniformReassortmentReheight,
        coalre.operators.ChangePredictorOperator,
        coalre.operators.EffectSizePredictorOperator,
        coalre.simulator.InitNetworkFromExtendedNewick,
        coalre.simulator.SimulatedCoalescentNetwork,
        coalre.simulator.SIRwithReassortment,
        coalre.simulator.SISwithReassortment,
        coalre.simulator.SuperspreadingSIRwithReassortment,
        coalre.simulator.SuperspreadingStructuredSIRwithReassortment,
        coalre.statistics.NetworkStatsLogger,
        coalre.statistics.ReassortmentEventsLogger,
        coalre.statistics.ReassortmentStatsLogger,
        coalre.util.DateOffsetInitializer,
        coalre.util.DummyTreeDistribution,
        coalre.util.SpikeSlabParameter,
        coalre.util.WeightedSumDistribution,
        coalre.util.InitFromTree,
        coalre.dynamics.RecombinationDynamicsFromSpline,
        coalre.dynamics.Spline,
        coalre.dynamics.NeDynamicsFromSpline,
        coalre.dynamics.SplineTransmissionDifference,
        coalre.dynamics.PiecewiseConstantReassortmentRates,
        coalre.dynamics.PiecewiseConstantReassortmentRateScalers,
        coalre.dynamics.Difference,
        coalre.dynamics.LogDifference,
        coalre.dynamics.SkygrowthNeDynamics,
        coalre.dynamics.SkygrowthReassortmentRatesFromSkygrowthNe,
        coalre.dynamics.GLMReassortmentRates;
}
