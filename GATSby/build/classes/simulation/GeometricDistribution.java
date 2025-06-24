/**
 * GeometricDistribution
    Copyright (C) 2015  Tom M. W. Nye

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
 * Simulate from a geometric distribution on 1,2,3,... with param p
 */

import cern.jet.random.tdouble.DoubleUniform;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import treebase.AlgorithmError;

public class GeometricDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Data */
    private double p;
    private DoubleUniform u;

    public GeometricDistribution(double q) {
        this(q, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public GeometricDistribution(double q, int seed) {
        this(q, new DoubleMersenneTwister(seed));
    }
    public GeometricDistribution(double q, DoubleMersenneTwister randomEngine) {
        u = new DoubleUniform(randomEngine);
        p=q;
    }
    
    public void setMean(double mu) throws AlgorithmError {
        if (mu<=1.0) {
            throw new AlgorithmError("Bad mean passed to GeometricDistribution.");
        }
        p = 1.0/mu;
    }
    
    public void setP(double q) throws AlgorithmError {
        if ((q>=1.0)||(q<=0.0)) {
            throw new AlgorithmError("Bad p passed to GeometricDistribution.");
        }
        p = q;
    }
    
    public void resetRandomEngineSeed() {
        u = new DoubleUniform(Random.getEngine());
    }

    public int sample() {
        double x = u.nextDouble();
        return (int) Math.ceil(Math.log(1-x) / Math.log(1-p));  
    }
    
    public double logpmf(int k) {
        return (k-1)*Math.log(1-p)+Math.log(p);
    }
    
    public double pmf(int k) {
        return Math.exp(logpmf(k));
    }
    
    public double cdf(int k) {
        return 1.0-Math.pow((1-p), k);
    }
    

}
