/**
 * NormalDistribution
    Copyright (C) 2012  Tom M. W. Nye

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

    Contact the author at:  <tom.nye@ncl.ac.uk>
                            <http://www.mas.ncl.ac.uk/~ntmwn/>
 */

package simulation;

/**
 * Normal distribution.
 * Version of colt normal distrib with quantiles added.
 */

import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class NormalDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Data */
    private double mu, sigma;
    private Normal coltObject;

    public NormalDistribution(double mean, double sd) {
        // Set up data
        mu = mean;
        sigma = Math.abs(sd);
        coltObject = new Normal(mu, sigma, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public NormalDistribution(double mean, double sd, int seed) {
        // Set up data
        mu = mean;
        sigma = Math.abs(sd);
        coltObject = new Normal(mu, sigma, new DoubleMersenneTwister(seed));
    }
    public NormalDistribution(double mean, double sd, DoubleMersenneTwister randomEngine) {
        // Set up data
        mu = mean;
        sigma = Math.abs(sd);
        coltObject = new Normal(mu, sigma, randomEngine);
    }
    
    public void resetRandomEngineSeed() {
        coltObject = new Normal(mu, sigma, Random.getEngine());
    }

    public void setParams(double mean, double sd)  {
        mu = mean;
        sigma = Math.abs(sd);
        coltObject.setState(mu, sigma);
    }
    
    public double getMu() {
        return mu;
    }
    
    public double getSigma() {
        return sigma;
    }

    /** log probability density function of the Normal distribution */
    public double logpdf(double x)
    {
        return -Math.log(2.0*Math.PI)/2.0 - Math.log(sigma) - (x-mu)*(x-mu)/(2.0*sigma*sigma);
    }

    /** probability density function of the Normal distribution */
    public double pdf(double x)
    {
        return Math.exp(this.logpdf(x));
    }
    
    /** cumulative distribution function of the Normal distribution */
    public double cdf(double x)
    {
        return coltObject.cdf(x);
    }
    
    /** log probability density function of the Normal distribution */
    public static double logpdf(double x, double mean, double sd)
    {
        return -Math.log(2.0*Math.PI)/2.0 - Math.log(sd) - (x-mean)*(x-mean)/(2.0*sd*sd);
    }

    
    /** probability density function of the Normal distribution */
    public static double pdf(double x, double mean, double sd)
    {
        return Math.exp(logpdf(x,mean,sd));
    }
    /** cumulative density function of the Normal distribution */
    public static double cdf(double x, double mean, double sd)
    {
        return cern.jet.stat.tdouble.Probability.normal((x-mean)/sd);
    }

    /** Random number generation */
    public double sample() {
        return coltObject.nextDouble();
    }
    /** Random number generation: bypass internal state values */
    public double sample(double mean, double sd) {
        return coltObject.nextDouble(mean, sd);
    }

    /** Quantiles */
    public static double quantile(double prob, double mean, double sd){
        return mean+sd*cern.jet.stat.tdouble.Probability.normalInverse(prob);
    }
    public double quantile(double prob){
        return mu+sigma*cern.jet.stat.tdouble.Probability.normalInverse(prob);
    }

    /* ---------------------------------------------------------------------- */
   
        
    }
