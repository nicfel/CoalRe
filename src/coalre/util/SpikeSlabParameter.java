package coalre.util;

import beast.base.core.BEASTObject;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;

public class SpikeSlabParameter extends BEASTObject implements Function {

    public Input<BooleanParameter> indicatorInput = new Input<>("indicator",
            "Boolean parameter indicating which elements are spike and which are slab.",
            Input.Validate.REQUIRED);

    public Input<Function> spikeValuesInput = new Input<>("spikeValues",
            "Value of each element when in spike mode.",
            Input.Validate.REQUIRED);

    public Input<Function> slabValuesInput = new Input<>("slabValues",
            "Value of each element when in slab mode.",
            Input.Validate.REQUIRED);

    Function spikeValues, slabValues;
    BooleanParameter indicators;

    SpikeSlabParameter() { }

    @Override
    public void initAndValidate() {

        indicators = indicatorInput.get();
        spikeValues = spikeValuesInput.get();
        slabValues = slabValuesInput.get();

        if (indicators.getDimension() != spikeValues.getDimension()
                || indicators.getDimension() != slabValues.getDimension()) {
            throw new IllegalArgumentException("Dimensions of all inputs to " +
                    "SpikeSlabParameter must match.");
        }

    }

    @Override
    public int getDimension() {
        return indicators.getDimension();
    }

    @Override
    public double getArrayValue(int i) {
        return indicators.getValue(i)
                ? spikeValues.getArrayValue(i)
                : slabValues.getArrayValue(i);
    }

    @Override
    public double getArrayValue() {
        return getArrayValue(0);
    }
}
