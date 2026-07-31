package coalre.dynamics;

import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

import java.util.ArrayList;
import java.util.List;

public class SplineTransmissionDifference extends CalculationNode implements RealVector<Real> {

    public Input<Spline> splineInput = new Input<>("spline", "spline to use for the population function", Input.Validate.REQUIRED);
    public Input<Spline> spline2Input = new Input<>("spline2", "spline to use for the population function");

    double[] difference;
    double[] storedDifference;
    boolean needsRecompute = true;

    @Override
    public void initAndValidate() {
        if (splineInput.get().infectedIsNe) {
            difference = new double[splineInput.get().splineCoeffs.length];
            storedDifference = new double[splineInput.get().splineCoeffs.length];
        }else {
            difference = new double[splineInput.get().splineCoeffs.length - 1];
            storedDifference = new double[splineInput.get().splineCoeffs.length - 1];
        }
    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }

    @Override
    public int size() {
        return difference.length;
    }

    @Override
    public List<Double> getElements() {
        if (needsRecompute) {
            compute();
        }
        List<Double> elements = new ArrayList<>(difference.length);
        for (double v : difference) elements.add(v);
        return elements;
    }

    @Override
    public double get(int dim) {
        if (needsRecompute) {
            compute();
        }
        return difference[dim];
    }

    void compute() {

        if (splineInput.get().infectedIsNe){
            double[] value = new double[splineInput.get().splineCoeffs.length+1];
            for (int i = 0; i <= splineInput.get().splineCoeffs.length; i++) {
                value[i] = splineInput.get().InfectedInput.get().get(i);
				if (spline2Input.get() != null) {
					value[i] += spline2Input.get().InfectedInput.get().get(i);
				}
            }

            for (int i = 1; i <= splineInput.get().splineCoeffs.length; i++) {
                difference[i - 1] = value[i - 1] - value[i];
            }

        }else {
            double[] transmissionRates = new double[splineInput.get().splineCoeffs.length];
            for (int i = 0; i < splineInput.get().splineCoeffs.length; i++) {
                transmissionRates[i] = splineInput.get().uninfectiousRate.get() -
                        splineInput.get().splineCoeffs[i][2];
            }

            for (int i = 1; i < splineInput.get().splineCoeffs.length; i++) {
                difference[i - 1] = transmissionRates[i - 1] - transmissionRates[i];
            }
        }

        needsRecompute = false;
    }

    @Override
    public void store() {
        System.arraycopy(difference, 0, storedDifference, 0, difference.length);
        super.store();
    }

    @Override
    public void restore() {
        double [] tmp = storedDifference;
        storedDifference = difference;
        difference = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }


}
