/**
 *  TruncatedNormalDistribution
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
 * Truncated Normal distribution.
 * Based on R's functions in the msm package. The sampling methods are, in turn,
 * based on the algorithms from Robert, C. P. Simulation of truncated normal 
 * variables. Statistics and Computing (1995) 5, 121-125
 */

import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.Exponential;
import cern.jet.random.tdouble.DoubleUniform;
import treebase.AlgorithmError;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class TruncatedNormalDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Data */
    private double mu, sigma;
    private Normal coltObject;
    private DoubleUniform uniformDistrib;
    private Exponential expDistrib;
    private double upper,lower;
    

    public TruncatedNormalDistribution(double mean, double sd, double low, double up) throws AlgorithmError {
        // Set up data
        mu = mean;
        sigma = Math.abs(sd);
        coltObject = new Normal(mu, sigma, Random.getEngine());
        if(low>=up) {
            throw new AlgorithmError("Lower truncation point must be less than upper trucation point in truncated normal distribution.");
        }
        lower = low;
        upper = up;
        expDistrib = new Exponential(1.0, Random.getEngine());
        uniformDistrib = new DoubleUniform(Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public TruncatedNormalDistribution(double mean, double sd, double low, double up, int seed) throws AlgorithmError {
        // Set up data
        mu = mean;
        sigma = Math.abs(sd);
        coltObject = new Normal(mu, sigma, new DoubleMersenneTwister(seed));
        if(low>=up) {
            throw new AlgorithmError("Lower truncation point must be less than upper trucation point in truncated normal distribution.");
        }
        lower = low;
        upper = up;
        expDistrib = new Exponential(1.0, new DoubleMersenneTwister(seed));
        uniformDistrib = new DoubleUniform(new DoubleMersenneTwister(seed));
    }
    public TruncatedNormalDistribution(double mean, double sd, double low, double up, DoubleMersenneTwister randomEngine) throws AlgorithmError {
        // Set up data
        mu = mean;
        sigma = Math.abs(sd);
        coltObject = new Normal(mu, sigma, randomEngine);
        if(low>=up) {
            throw new AlgorithmError("Lower truncation point must be less than upper trucation point in truncated normal distribution.");
        }
        lower = low;
        upper = up;
        expDistrib = new Exponential(1.0, randomEngine);
        uniformDistrib = new DoubleUniform(randomEngine);
    }
    
    public void resetRandomEngineSeed() {
        coltObject = new Normal(mu, sigma, Random.getEngine());
        expDistrib = new Exponential(1.0, Random.getEngine());
        uniformDistrib = new DoubleUniform(Random.getEngine());
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

    public void setTruncationPoints(double low, double up) throws AlgorithmError {
        if(low>=up) {
            throw new AlgorithmError("Lower truncation point must be less than upper trucation point in truncated normal distribution.");
        }
        lower = low;
        upper = up;
    }
    
    public double getLower() {
        return lower;
    }
    
    public double getUpper() {
        return upper;
    }

    /** log probability density function of the truncated Normal distribution */
    public double logpdf(double x)
    {
        if( (x>=lower) & (x<=upper) )
            return -Math.log(2.0*Math.PI)/2.0 - Math.log(sigma) - (x-mu)*(x-mu)/(2.0*sigma*sigma) - Math.log(coltObject.cdf(upper) - coltObject.cdf(lower));
        else return Double.NEGATIVE_INFINITY;
    }
    
    /** probability density function of the truncated Normal distribution */
    public double pdf(double x)
    {
        return Math.exp(this.logpdf(x));
    }
    
    /** cumulative distribution function of the truncated Normal distribution */
    public double cdf(double x)
    {
        if( (x>=lower) & (x<=upper) ) {
            double denom = coltObject.cdf(upper) - coltObject.cdf(lower);
            double numer = coltObject.cdf(x) - coltObject.cdf(lower);
            return numer/denom;
        } else if(x < lower) {
            return 0.0;
        } else {
            return 1.0;
        }
    }
    
    /** log probability density function of the truncated Normal distribution */
    public static double logpdf(double x, double mean, double sd, double low, double up) throws AlgorithmError
    {
        if(low>=up) {
            throw new AlgorithmError("Lower truncation point must be less than upper trucation point in truncated normal distribution.");
        }
        if( (x>=low) & (x<=up) )
           return NormalDistribution.logpdf(x,mean,sd) - Math.log(NormalDistribution.cdf(up,mean,sd) - NormalDistribution.cdf(low,mean,sd));
        else return Double.NEGATIVE_INFINITY;
    }

    /** probability density function of the truncated Normal distribution */
    public static double pdf(double x, double mean, double sd, double low, double up) throws AlgorithmError
    {
        return Math.exp(logpdf(x,mean,sd,low,up));
    }
    /** cumulative density function of the truncated Normal distribution */
    public static double cdf(double x, double mean, double sd, double low, double up) throws AlgorithmError
    {
        if(low>=up) {
            throw new AlgorithmError("Lower truncation point must be less than upper trucation point in truncated normal distribution.");
        }
        if( (x>=low) & (x<=up) ) {
            double denom = NormalDistribution.cdf(up,mean,sd) - NormalDistribution.cdf(low,mean,sd);
            double numer = NormalDistribution.cdf(x,mean,sd) - NormalDistribution.cdf(low,mean,sd);
            return numer/denom;
        } else if(x < low) {
            return 0.0;
        } else {
            return 1.0;
        }
    }

    /** Random number generation */
    public double sample() {
        double mean = getMu();
        double sd = getSigma();
        double low = (lower - mean)/sd;
        double up = (upper - mean)/sd;
        double z;
        if( (low < 0 & Double.isInfinite(up)) | (up > 0 & Double.isInfinite(low)) | (!Double.isInfinite(low) & !Double.isInfinite(up) & low < 0 & up > 0 & (up - low) > Math.sqrt(2*Math.PI)) ) {
            z=sampleNormal(mean,sd,low,up);
        } else if( low >= 0 & up > (low + 2.0*Math.sqrt(Math.E)*Math.exp((low*low-low*Math.sqrt(low*low+4.0))/4.0)/(low+Math.sqrt(low*low+4.0))) ) {
            z=sampleLowExponential(mean,sd,low,up);
        } else if( up <= 0 & -low > (-up + 2.0*Math.sqrt(Math.E)*Math.exp((up*up-(-up)*Math.sqrt(up*up+4.0))/4.0)/(-up+Math.sqrt(up*up+4.0))) ) {
            z=sampleUpExponential(mean,sd,low,up);
        } else {
            z=sampleUniform(mean,sd,low,up);
        }
        return z*sd+mean;
    }
    
    public double sampleNormal(double mean, double sd, double low, double up) {
        double y;
        do {
            y=coltObject.nextDouble(0.0,1.0);
        } while( y < low | y > up );
        return y;
    }
    
    public double sampleLowExponential(double mean, double sd, double low, double up) {
        double a = (low+Math.sqrt(low*low+4.0))/2.0;
        double y,u;
        do {
            y=expDistrib.nextDouble(a)+low;
            u=uniformDistrib.nextDouble();
        } while( u > Math.exp(-(y-a)*(y-a)/2.0) | y > up );
        return y;
    }
    
    public double sampleUpExponential(double mean, double sd, double low, double up) {
        double a = (-up+Math.sqrt(up*up+4.0))/2.0;
        double y,u;
        do {
            y=expDistrib.nextDouble(a)-up;
            u=uniformDistrib.nextDouble();
        } while( u > Math.exp(-(y-a)*(y-a)/2.0) | y > -low );
        return y;
    }
    
    public double sampleUniform(double mean, double sd, double low, double up) {
        double rho,u,y;
        do {
            y=uniformDistrib.nextDoubleFromTo(low,up);
            if(low > 0) {
                rho=Math.exp((low*low-y*y)/2.0);
            } else if(up < 0) {
                rho=Math.exp((up*up-y*y)/2.0);
            } else {
                rho=Math.exp(-y*y/2.0);
            }
            u=uniformDistrib.nextDouble();
        } while( u > rho );
        return y;
    }
    
    /** Random number generation: bypass internal state values */
    public double sample(double mean, double sd, double low, double up) {
        low = (low - mean)/sd;
        up = (up - mean)/sd;
        double z;
        if( (low < 0 & Double.isInfinite(up)) | (up > 0 & Double.isInfinite(low)) | (!Double.isInfinite(low) & !Double.isInfinite(up) & low < 0 & up > 0 & (up - low) > Math.sqrt(2*Math.PI)) ) {
            z=sampleNormal(mean,sd,low,up);
        } else if( low >= 0 & up > (low + 2.0*Math.sqrt(Math.E)*Math.exp((low*low-low*Math.sqrt(low*low+4.0))/4.0)/(low+Math.sqrt(low*low+4.0))) ) {
            z=sampleLowExponential(mean,sd,low,up);
        } else if( up <= 0 & -low > (-up + 2.0*Math.sqrt(Math.E)*Math.exp((up*up-(-up)*Math.sqrt(up*up+4.0))/4.0)/(-up+Math.sqrt(up*up+4.0))) ) {
            z=sampleUpExponential(mean,sd,low,up);
        } else {
            z=sampleUniform(mean,sd,low,up);
        }
        return z*sd+mean;
    }


}
