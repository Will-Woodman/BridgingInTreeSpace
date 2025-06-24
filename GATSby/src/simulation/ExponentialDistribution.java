/**
 *  ExponentialDistribution
    Copyright (C) 2012  Sarah E. Heaps

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact the author at:  <sarah.heaps@ncl.ac.uk>
 */

package simulation;

/**
 * Exponential distribution.
 * 
 */

import cern.jet.random.tdouble.Exponential;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class ExponentialDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    /* Instance variables */
    private double rate;
    private Exponential coltObject;
    
    /** Constructor */
    public ExponentialDistribution(double lambda) {
        rate = lambda;
        coltObject = new Exponential(rate, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public ExponentialDistribution(double lambda, int seed) {
        rate = lambda;
        coltObject = new Exponential(rate, new DoubleMersenneTwister(seed));
    }
    public ExponentialDistribution(double lambda, DoubleMersenneTwister randomEngine) {
        rate = lambda;
        coltObject = new Exponential(rate, randomEngine);
    }
    
    public void resetRandomEngineSeed() {
        coltObject = new Exponential(rate, Random.getEngine());
    }
    
    public void setParams(double lambda)  {
        rate = lambda;
        coltObject.setState(rate);
    }

    public double getRate() {
        return rate;
    }
    
    /** log probability density function of the exponential distribution */
    public double logpdf(double x)
    {
        return Math.log(rate) - rate*x;
    }

    /** probability density function of the exponential distribution */
    public double pdf(double x)
    {
        return Math.exp(this.logpdf(x));
    }
    
    /** cumulative distribution function of the exponential distribution */
    public double cdf(double x)
    {
        return 1.0 - Math.exp(-rate*x);
    }
    
    /** log probability density function of the exponential distribution */
    public static double logpdf(double x, double lambda)
    {
        return Math.log(lambda) - lambda*x;
    }

    /** probability density function of the exponential distribution */
    public static double pdf(double x, double lambda)
    {
        return Math.exp(logpdf(x,lambda));
    }
    
    /** cumulative density function of the exponential distribution */
    public static double cdf(double x, double lambda)
    {
        return 1.0 - Math.exp(-lambda*x);
    }

    /** Random number generation */
    public double sample() {
        return coltObject.nextDouble();
    }
    /** Random number generation: bypass internal state values */
    public double sample(double lambda) {
        return coltObject.nextDouble(lambda);
    }
    
    /** Quantiles */
    public static double quantile(double prob, double lambda){
        return -Math.log(1.0 - prob)/lambda;
    }
    public double quantile(double prob){
        return -Math.log(1.0 - prob)/rate;
    }
    

    
}
