package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;
import java.util.List;

/**
 * @author Nicola F. Mueller
 */
@Description("Population function with defines Ne's at points in time and interpolated between them. Parameter has to be in log space. The Ne's are used to compute the transmission rates, which are used for the coalescent process. ")
public class SkygrowthReassortmentRatesFromSkygrowthNe extends PopulationFunction.Abstract {

	final public Input<RealParameter> logNeInput = new Input<>("logNe", "Nes over time in log space");
	final public Input<RealParameter> NeToReassortmentInput = new Input<>("neToReassortment",
			"the value that maps the number of infected or the Ne to the reassortment rate ");
	final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
			"When to switch between elements of Ne", Input.Validate.REQUIRED);

	RealParameter Ne;
	RealParameter NeToReassortment;
	RealParameter rateShifts;

	boolean NesKnown = false;
	double[] growth;
	double[] growth_stored;

	double[] rates;
	double[] stored_rates;
	
	@Override
	public void initAndValidate() {
		
		NeToReassortment = NeToReassortmentInput.get();
		rateShifts = rateShiftsInput.get();
		if (logNeInput.get()!=null) {
			Ne = logNeInput.get();
			Ne.setDimension(rateShifts.getDimension());
		}
		NeToReassortment.setDimension(rateShifts.getDimension());
		growth = new double[rateShifts.getDimension()];
		
		rates = new double[rateShifts.getDimension()];
		stored_rates = new double[rateShifts.getDimension()];
		
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
		timediff -= rateShifts.getValue(i);
		
		return Math.exp(rates[i] - growth[i] * timediff);
	}

	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.getDimension()-1; i++)
			if (t < rateShifts.getValue(i+1))
				return i;

		// after the last interval, just keep using the last element
		return rateShifts.getDimension()-1;
	}

	@Override
	public double getIntegral(double start, double finish) {
		if (start == finish)
			return 0.0;

		// get the interval "start" is in
		int first_int = getIntervalNr(start);
		// get the interval "finish" is in
		int last_int = getIntervalNr(finish);

		double weighted = 0.0;
		double curr_time = start;	

		for (int i = first_int; i <= last_int; i++) {
			if (i > rateShifts.getDimension()) {
				throw new IllegalArgumentException("rate shifts out of bounds");
			}

			double next_time = Math.min(getNextTime(i), finish);
			double r = growth[i];

			double timediff1 = curr_time;
			double timediff2 = next_time;
			timediff1 -= rateShifts.getArrayValue(i);
			timediff2 -= rateShifts.getArrayValue(i);
			
			if (r == 0.0) {
				weighted += (next_time - curr_time) * Math.exp(rates[i]);
			} else {
				weighted += Math.exp(rates[i]) / -r * ( Math.exp(-timediff2 * r) - Math.exp(-timediff1 * r));
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

		int i = 0;
		double curr_time = 0;
		double integral = 0;
		
		do
		{
			double next_time = getNextTime(i);
			double r = growth[i];
	
//			double timediff1 = curr_time-rateShifts.getArrayValue(i);
			double timediff2 = next_time-rateShifts.getArrayValue(i);
	
			double old_diff = x - integral;
	
			if (r == 0.0) {
				integral += (next_time - curr_time) * Math.exp(rates[i]);
			} else {
				integral +=  Math.exp(rates[i]) / -r * (Math.exp(-timediff2 * r)-1.0);
			}
			
	
			double diff = x - integral;
	
			if (diff < 0 || i == rateShifts.getDimension()) {
				
				if (r == 0.0) {
					return old_diff/Math.exp(rates[i]) + curr_time;
				} else {
					return Math.log((old_diff * -r)/ Math.exp(rates[i]) +1.0)/-r  + curr_time;
				}
			
			}
			
			curr_time = next_time;
			i++;
			
		}while(i<=rateShifts.getDimension());
		
		
	
		return Double.POSITIVE_INFINITY;

	}


	private double getNextTime(int i) {
		if (i < rateShifts.getDimension()-1)
			return rateShifts.getArrayValue(i+1);
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
		stored_rates = new double[rates.length];
		System.arraycopy(rates, 0, stored_rates, 0, rates.length);
		super.store();
	}

	@Override
	public void restore() {
		System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
		System.arraycopy(stored_rates, 0, rates, 0, stored_rates.length);
		super.restore();
	}

	// computes the Ne's at the break points
	private void recalculateNe() {
		growth = new double[rateShifts.getDimension()];
		rates = new double[rateShifts.getDimension()];
		double curr_time = 0.0;
		if (logNeInput.get()!=null) {			
			for (int i = 0; i < Ne.getDimension(); i++) {
				rates[i] = Ne.getArrayValue(i) + NeToReassortment.getArrayValue(i);
			}
		}else {
			for (int i = 0; i < NeToReassortment.getDimension()-1; i++) {
				rates[i] = NeToReassortment.getArrayValue(i);
			}
		}
		for (int i = 0; i < NeToReassortment.getDimension()-1; i++) {
			growth[i] = (rates[i] - rates[i+1])
					/(rateShifts.getValue(i+1)-curr_time);
			curr_time = rateShifts.getValue(i+1);
		}
		NesKnown = true;
	}
}