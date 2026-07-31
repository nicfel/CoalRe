package coalre.util;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.spec.domain.Real;
import beast.base.spec.type.BoolVector;
import beast.base.spec.type.RealVector;

import java.util.AbstractList;
import java.util.List;

/**
 * A derived real vector whose ith element is taken from either the spike or the
 * slab vector, selected by the ith indicator.
 */
public class SpikeSlabParameter extends BEASTObject implements RealVector<Real> {

    public Input<BoolVector> indicatorInput = new Input<>("indicator",
            "Boolean parameter indicating which elements are spike and which are slab.",
            Input.Validate.REQUIRED);

    public Input<RealVector<Real>> spikeValuesInput = new Input<>("spikeValues",
            "Value of each element when in spike mode.",
            Input.Validate.REQUIRED);

    public Input<RealVector<Real>> slabValuesInput = new Input<>("slabValues",
            "Value of each element when in slab mode.",
            Input.Validate.REQUIRED);

    RealVector<Real> spikeValues, slabValues;
    BoolVector indicators;

    SpikeSlabParameter() { }

    @Override
    public void initAndValidate() {

        indicators = indicatorInput.get();
        spikeValues = spikeValuesInput.get();
        slabValues = slabValuesInput.get();

        if (indicators.size() != spikeValues.size()
                || indicators.size() != slabValues.size()) {
            throw new IllegalArgumentException("Dimensions of all inputs to " +
                    "SpikeSlabParameter must match.");
        }

    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }

    @Override
    public double get(int i) {
        return indicators.get(i)
                ? spikeValues.get(i)
                : slabValues.get(i);
    }

    /**
     * Elements are derived on demand rather than stored, so this is a live view
     * over the spike/slab inputs instead of a copied list.
     */
    @Override
    public List<Double> getElements() {
        return new AbstractList<>() {
            @Override
            public Double get(int i) {
                return SpikeSlabParameter.this.get(i);
            }

            @Override
            public int size() {
                return indicators.size();
            }
        };
    }
}
