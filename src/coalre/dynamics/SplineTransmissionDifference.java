package coalre.dynamics;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;

public class SplineTransmissionDifference extends CalculationNode implements Function {

    public Input<Spline> splineInput = new Input<>("spline", "spline to use for the population function", Input.Validate.REQUIRED);

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
    public int getDimension() {
        return difference.length;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return difference[0];
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return difference[dim];
    }

    void compute() {

        if (splineInput.get().infectedIsNe){
            double[] value = new double[splineInput.get().splineCoeffs.length+1];
            for (int i = 0; i <= splineInput.get().splineCoeffs.length; i++) {
                value[i] = splineInput.get().InfectedInput.get().getArrayValue(i);
            }

            for (int i = 1; i <= splineInput.get().splineCoeffs.length; i++) {
                difference[i - 1] = value[i - 1] - value[i];
            }

        }else {
            double[] transmissionRates = new double[splineInput.get().splineCoeffs.length];
            for (int i = 0; i < splineInput.get().splineCoeffs.length; i++) {
                transmissionRates[i] = splineInput.get().uninfectiousRate.getValue() -
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
