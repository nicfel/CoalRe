package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.util.List;


/**
 * @author Nicola F. Mueller
 */
@Description("Computes time varying recombination rates as a population"+
        " function from spline interpolation of the number of infected over time")
public class RecombinationDynamicsFromSpline extends PopulationFunction.Abstract implements Loggable {
    final public Input<RealParameter> InfectedToRhoInput = new Input<>("InfectedToRho",
            "the value that maps the number of infected or the Ne to the reassortment rate ", Input.Validate.REQUIRED);
    final public Input<Spline> splineInput = new Input<>("spline",
            "Spline to use for the population function", Input.Validate.REQUIRED);

    boolean NesKnown = false;

    Spline spline;

    RealParameter InfectedToRho;

    double rateRatio;

    @Override
    public void initAndValidate() {
        InfectedToRho = InfectedToRhoInput.get();
        spline = splineInput.get();
    }


    @Override
    public List<String> getParameterIds() {
        return null;
    }

    @Override
    public double getPopSize(double t) {
        if (!spline.update())
            return Double.NaN;

        // check which time t is, if it is larger than the last time, return the last Ne
        int interval = spline.gridPoints-1;
        for (int i = 0; i < spline.gridPoints; i++){
            if (t < spline.time[i]){
                interval = i-1;
                break;
            }
        }
        return spline.I[interval]*InfectedToRho.getValue();
    }

    public double getIntegral(double from, double to) {
        if (!spline.update())
            return Double.NaN;

        // compute the integral of Ne's between from an to
        double NeIntegral = 0;
        int intervalFrom = spline.gridPoints-1;
        for (int i = 0; i < spline.gridPoints; i++) {
            // get the first time larger than from
            if (from < spline.time[i]) {
                intervalFrom = i-1;
                break;
            }
        }
        for (int i = intervalFrom; i < spline.gridPoints; i++) {
            double rate = spline.I[i]*InfectedToRho.getValue();
            // if i==intervalFrom, we have to start compute the diff from there
            if (spline.time[i+1] > to) {
                if (i == intervalFrom) {
                    NeIntegral += rate * (to - from);
                }else{
                    NeIntegral += rate * (to - spline.time[i]);
                }
                break;
            }else if (i == intervalFrom) {
                NeIntegral += rate * (spline.time[i+1] - from);
            }else{
                NeIntegral += rate * (spline.time[i + 1] - spline.time[i]);
            }
        }
        return NeIntegral;
    }

    @Override
    public double getIntensity(double v){
        return getIntegral(0, v);
    }

    @Override
    public double getInverseIntensity(double v) {
        for (int i = 0; i < spline.gridPoints; i++) {
            double rate = spline.I[i]*InfectedToRho.getValue();
            v -= rate * (spline.time[i + 1] - spline.time[i]);
            if (v<0){
                v += rate * (spline.time[i + 1] - spline.time[i]);
                // solve for the final time
                return spline.time[i] + v/rate;
            }
        }
        return spline.time[spline.gridPoints-1] + v/(spline.I[spline.gridPoints-1]*InfectedToRho.getValue());
    }

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    public void init(PrintStream printStream) {
        for (int i = 0; i < spline.gridPoints; i+=20) {
            printStream.print("reassortment" + i + "\t");
        }
    }

    @Override
    public void log(long l, PrintStream printStream) {
        for (int i = 0; i < spline.gridPoints; i+=20) {
            printStream.print(spline.I[i]*InfectedToRho.getValue() + "\t");
        }
    }

    @Override
    public void close(PrintStream printStream) {

    }


}