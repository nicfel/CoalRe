package coalre.operators;

import beast.util.Randomizer;

public class MergeSplitOperator extends NetworkOperator {

    @Override
    protected double networkProposal() {

        if (Randomizer.nextBoolean())
            return split();
        else
            return merge();
    }

    protected double split() {
        return 0.0;
    }

    protected double merge() {
        return 0.0;
    }
}
