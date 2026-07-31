package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;

import java.util.Arrays;
import java.util.List;

/**
 * @author Nicola F. Mueller
 */
@Description("Population function with defines Ne's at points in time and interpolated between them. Parameter has to be in log space. The Ne's are used to compute the transmission rates, which are used for the coalescent process. ")
public class SkygrowthNeDynamics extends PopulationFunction.Abstract {

	// written to (setDimension) in initAndValidate: needs the concrete param type
	final public Input<RealVectorParam<Real>> logNeInput = new Input<>("logNe", "Nes over time in log space",
			Input.Validate.REQUIRED);
	final public Input<beast.base.spec.type.RealVector<NonNegativeReal>> rateShiftsInput = new Input<>("rateShifts",
			"When to switch between elements of Ne", Input.Validate.REQUIRED);

	RealVectorParam<Real> Ne;
	beast.base.spec.type.RealVector<NonNegativeReal> rateShifts;

	boolean NesKnown = false;
	double[] growth;
	double[] growth_stored;
	double[] invRatio;
	double[] invRatio_stored;

	@Override
	public void initAndValidate() {
		Ne = logNeInput.get();
		rateShifts = rateShiftsInput.get();
		Ne.setDimension(rateShifts.size());
		growth = new double[rateShifts.size()];
		invRatio = new double[rateShifts.size()];
		recalculateNe();
	}

	@Override
	public List<String> getParameterIds() {
		return null;
	}

	@Override
	public double getPopSize(double t) {
		int i = getIntervalNr(t);
		double timediff = t;
		timediff -= rateShifts.get(i);

		return Math.exp(Ne.get(i) - growth[i] * timediff);
	}

	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.size()-1; i++)
			if (t < rateShifts.get(i+1))
				return i;

		// after the last interval, just keep using the last element
		return rateShifts.size()-1;
	}

	private int getLaterIntervalNr(double t, int startPoint) {
		// check which interval t + offset is in
		for (int i = startPoint; i < rateShifts.size()-1; i++)
			if (t < rateShifts.get(i+1))
				return i;
		// after the last interval, just keep using the last element
		return rateShifts.size()-1;
	}

	@Override
	public double getIntegral(double start, double finish) {
		if (start == finish)
			return 0.0;

		// get the interval "start" is in
		int first_int = getIntervalNr(start);
		// get the interval "finish" is in
		int last_int = getLaterIntervalNr(finish, first_int);

		double weighted = 0.0;
		double curr_time = start;	

		for (int i = first_int; i <= last_int; i++) {
			if (i > rateShifts.size()) {
				throw new IllegalArgumentException("rate shifts out of bounds");
			}

			double next_time = Math.min(getNextTime(i), finish);
			double r = growth[i];

			double rateShift = rateShifts.get(i);
			double timediff1 = curr_time - rateShift;
			double timediff2 = next_time - rateShift;
			
			if (r == 0.0) {
				weighted += (next_time - curr_time) / Math.exp(Ne.get(i));
			} else {
				weighted += (Math.exp(timediff2 * r) - Math.exp(timediff1 * r)) / invRatio[i];
			}

			curr_time = next_time;
		}
		return weighted;

	}
	
	
	@Override
	public double getIntensity(double t) {
		return getIntegral(0,t);
	}


	@Override
	public double getInverseIntensity(double x) {
		
		if (x == Double.POSITIVE_INFINITY) {
			return Double.POSITIVE_INFINITY;
		}

		int i = 0;
		double curr_time = 0;
		double integral = 0;
		
		do
		{
			double next_time = getNextTime(i);
			double r = growth[i];
	
			double timediff1 = curr_time-rateShifts.get(i);
			double timediff2 = next_time-rateShifts.get(i);
	
			double old_diff = x - integral;
	
			if (r == 0.0) {
				integral += (next_time - curr_time) / Math.exp(Ne.get(i));
			} else {
				integral += (Math.exp(timediff2 * r) - Math.exp(timediff1 * r)) / invRatio[i];
			}
	
			double diff = x - integral;
	
			if (diff < 0 || i == rateShifts.size()) {
				
				if (r == 0.0) {
					return Math.exp(Ne.get(i)) * old_diff + curr_time;
				} else {
					return Math.log(Math.exp(Ne.get(i)) * old_diff * r + Math.exp(timediff1 * r)) / r + curr_time;
				}
			
			}
			
			curr_time = next_time;
			i++;
			
		}while(i<=rateShifts.size());
		
		
	
		return Double.POSITIVE_INFINITY;

	}


	private double getNextTime(int i) {
		if (i < rateShifts.size()-1)
			return rateShifts.get(i+1);
		else
			return Double.POSITIVE_INFINITY;
	}


	@Override
	public boolean requiresRecalculation() {
		recalculateNe();
		return super.requiresRecalculation();
	}

	@Override
	public void store() {
		growth_stored = new double[growth.length];
		System.arraycopy(growth, 0, growth_stored, 0, growth.length);
		invRatio_stored = new double[invRatio.length];
		System.arraycopy(invRatio, 0, invRatio_stored, 0, invRatio.length);
		super.store();
	}

	@Override
	public void restore() {
		System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
		System.arraycopy(invRatio_stored, 0, invRatio, 0, invRatio_stored.length);
		super.restore();
	}

	// computes the Ne's at the break points
	private void recalculateNe() {
		growth = new double[rateShifts.size()];
		double curr_time = 0.0;
		for (int i = 0; i < Ne.size()-1; i++) {
			growth[i] = (Ne.get(i) - Ne.get(i+1))/(rateShifts.get(i+1)-curr_time);
			curr_time = rateShifts.get(i+1);
			invRatio[i] = Math.exp(Ne.get(i)) * growth[i];
			
		}
		NesKnown = true;
	}
}