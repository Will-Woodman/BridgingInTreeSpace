/**
 *  LogNormalDistribution
    Copyright (C) 2012  Sarah E Heaps

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
 * Log-normal distribution.
 * 
 */

import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class LogNormalDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    /* Data */
    private double mu, sigma;
    private Normal coltObject;
    
    /** Constructor */
    public LogNormalDistribution(double mlog, double sdlog) {
        mu = mlog;
        sigma = Math.abs(sdlog);
        coltObject = new Normal(mu, sigma, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public LogNormalDistribution(double mlog, double sdlog, int seed) {
        mu = mlog;
        sigma = Math.abs(sdlog);
        coltObject = new Normal(mu, sigma, new DoubleMersenneTwister(seed));
    }
    public LogNormalDistribution(double mlog, double sdlog, DoubleMersenneTwister randomEngine) {
        mu = mlog;
        sigma = Math.abs(sdlog);
        coltObject = new Normal(mu, sigma, randomEngine);
    }
    
    public void resetRandomEngineSeed() {
        coltObject = new Normal(mu, sigma, Random.getEngine());
    }
    
    public void setParams(double mlog, double sdlog)  {
        mu = mlog;
        sigma = Math.abs(sdlog);
        coltObject.setState(mu, sigma);
    }
    
    public double getMu() {
        return mu;
    }
    
    public double getSigma() {
        return sigma;
    }
    
    /** log probability density function of the Log-normal distribution */
    public double logpdf(double x)
    {
        return -Math.log(2.0*Math.PI)/2.0 - Math.log(sigma) - (Math.log(x)-mu)*(Math.log(x)-mu)/(2.0*sigma*sigma) - Math.log(x);
    }

    /** probability density function of the Log-normal distribution */
    public double pdf(double x)
    {
        return Math.exp(this.logpdf(x));
    }
    
    /** cumulative distribution function of the Log-normal distribution */
    public double cdf(double x)
    {
        return coltObject.cdf(Math.log(x));
    }
    
    /** log probability density function of the Log-normal distribution */
    public static double logpdf(double x, double mean, double sd)
    {
        return NormalDistribution.logpdf(Math.log(x),mean,sd) - Math.log(x);
    }

    /** probability density function of the Log-normal distribution */
    public static double pdf(double x, double mean, double sd)
    {
        return Math.exp(logpdf(x,mean,sd));
    }
    
    /** cumulative density function of the Log-normal distribution */
    public static double cdf(double x, double mean, double sd)
    {
        return NormalDistribution.cdf(Math.log(x),mean,sd);
    }

    /** Random number generation */
    public double sample() {
        return Math.exp(coltObject.nextDouble());
    }
    /** Random number generation: bypass internal state values */
    public double sample(double mean, double sd) {
        return Math.exp(coltObject.nextDouble(mean, sd));
    }
    
    /** Quantiles */
    public static double quantile(double prob, double mean, double sd){
        return Math.exp(NormalDistribution.quantile(prob,mean,sd));
    }
    public double quantile(double prob){
        return Math.exp(mu+sigma*cern.jet.stat.tdouble.Probability.normalInverse(prob));
    }
    
    
}
