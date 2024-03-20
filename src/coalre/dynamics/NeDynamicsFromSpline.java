package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;

import java.io.PrintStream;
import java.util.List;


/**
 * @author Nicola F. Mueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class NeDynamicsFromSpline extends PopulationFunction.Abstract implements Loggable {

    final public Input<Spline> splineInput = new Input<>("spline",
            "Spline to use for the population function", Input.Validate.REQUIRED);

    boolean NesKnown = false;

    Spline spline;

    boolean returnNaN = false;

    @Override
    public void initAndValidate() {
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
        return spline.I[interval]/(spline.transmissionRate[interval]*2);
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
            // if i==intervalFrom, we have to start compute the diff from there
            if (spline.time[i+1] > to) {
                if (i == intervalFrom) {
                    NeIntegral += (spline.transmissionRate[i]/spline.I[i]) * (to - from);
                }else{
                    NeIntegral += (spline.transmissionRate[i]/spline.I[i]) * (to - spline.time[i]);
                }
                break;
            }else if (i == intervalFrom) {
                NeIntegral += (spline.transmissionRate[i]/spline.I[i]) * (spline.time[i+1] - from);
            }else{
                NeIntegral += (spline.transmissionRate[i] / spline.I[i]) * (spline.time[i + 1] - spline.time[i]);
            }
        }
        return 2*NeIntegral;
    }

    @Override
    public double getIntensity(double v) {
        return getIntegral(0,v);
    }

    @Override
    public double getInverseIntensity(double v) {
        // divide by 2 to avoid a division in every step
        v/=2;
        for (int i = 0; i < spline.gridPoints; i++) {
            v -= (spline.transmissionRate[i] / spline.I[i]) * (spline.time[i + 1] - spline.time[i]);
            if (v<0){
                v += (spline.transmissionRate[i] / spline.I[i]) * (spline.time[i + 1] - spline.time[i]);
                // solve for the final time
                return spline.time[i] + v/(spline.transmissionRate[i] / spline.I[i]);
            }
        }
        return spline.time[spline.gridPoints-1] + v/(spline.transmissionRate[spline.gridPoints-1] / spline.I[spline.gridPoints-1]);

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
            printStream.print("logNe_" + i + "\t");
        }
        for (int i = 0; i < spline.gridPoints; i+=20) {
            printStream.print("logI_" + i + "\t");
        }
        for (int i = 0; i < spline.gridPoints; i+=1) {
            printStream.print("transmissionRate" + i + "\t");
        }

    }

    @Override
    public void log(long l, PrintStream printStream) {
        for (int i = 0; i < spline.gridPoints; i+=20) {
            printStream.print(Math.log(spline.I[i]/spline.transmissionRate[i]) + "\t");
        }
        for (int i = 0; i < spline.gridPoints; i+=20) {
            printStream.print(spline.I[i] + "\t");
        }
        for (int i = 0; i < spline.gridPoints; i+=1) {
            printStream.print(spline.transmissionRate[i] + "\t");
        }

    }

    @Override
    public void close(PrintStream printStream) {

    }
}