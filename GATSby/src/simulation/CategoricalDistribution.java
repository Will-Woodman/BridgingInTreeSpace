/**
 * CategoricalDistribution
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
 * Simulate from a distribution on 0,..,k-1 with probs p1,...,pk
 */

import cern.jet.random.tdouble.DoubleUniform;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class CategoricalDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /* Data */
    private double[] cumsum;
    private double[] probabilities;
    private DoubleUniform u;

    public CategoricalDistribution(double[] probs) {
        this(probs, Random.getEngine());
    }
    
    /* Constructor which initialise the probability vector to the normalised
     unit vector */
    public CategoricalDistribution(int n) {
        this(n, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public CategoricalDistribution(double[] probs, int seed) {
        this(probs,new DoubleMersenneTwister(seed));
    }
    public CategoricalDistribution(double[] probs, DoubleMersenneTwister randomEngine) {
        u = new DoubleUniform(randomEngine);
        cumsum = new double[probs.length];
        double s = 0.0;
        for (int i=0; i<probs.length; i++) {
            s += probs[i];
            cumsum[i] = s;
        }
        probabilities = new double[probs.length];
        System.arraycopy(probs, 0, probabilities, 0, probs.length);
    }
    public CategoricalDistribution(int n, int seed) {
        this(n, new DoubleMersenneTwister(seed));
    }
    public CategoricalDistribution(int n, DoubleMersenneTwister randomEngine) {
        u = new DoubleUniform(randomEngine);
        cumsum = new double[n];
        probabilities = new double[n];
        for (int i=0; i<n; i++) {
            cumsum[i] = (i+1.0)/n;
            probabilities[i] = 1.0/n;
        }
    }
    
    public void resetRandomEngineSeed() {
        u = new DoubleUniform(Random.getEngine());
    }

    public int sample() {
        double x = u.nextDouble();
         for (int i=0; i<cumsum.length; i++) {
            if (x<=cumsum[i]) return i;
        }
        return cumsum.length-1;
    }

    public void changeProbabilities(double[] newProbs) throws treebase.AlgorithmError {
        if (newProbs.length!=probabilities.length) throw new treebase.AlgorithmError("Probabilities wrong length in CategoricalDistribution.");
        double s = 0.0;
        for (int i=0; i<newProbs.length; i++) {
            s += newProbs[i];
            cumsum[i] = s;
        }
        System.arraycopy(newProbs, 0, probabilities, 0, newProbs.length);
    }

    public double logpmf(int i) {
        if(i<0 || i>=probabilities.length) return Double.NEGATIVE_INFINITY;
        else return Math.log(probabilities[i]);
    }
    
    public double logcdf(int i) {
        if(i<0) return Double.NEGATIVE_INFINITY;
        else if(i>=probabilities.length) return 0.0;
        else return Math.log(cumsum[i]);
    }
    
    public double pmf(int i) {
        return Math.exp(logpmf(i));
    }
    
    public double cdf(int i) {
        return Math.exp(logcdf(i));
    }
    
    /* Useful sampler utils! */

    public static Object sampleFromSetWithoutReplacement(java.util.HashSet theSet, DoubleUniform theRnd) {
        java.util.ArrayList temp = new java.util.ArrayList();
        temp.addAll(theSet);
        int n = theRnd.nextIntFromTo(0,temp.size()-1);
        Object o = temp.get(n);
        theSet.remove(o);
        return o;
    }

    public static Object sampleFromSetWithReplacement(java.util.HashSet theSet, DoubleUniform theRnd) {
        java.util.ArrayList temp = new java.util.ArrayList();
        temp.addAll(theSet);
        int n = theRnd.nextIntFromTo(0,temp.size()-1);
        Object o = temp.get(n);
        return o;
    }

    public static Object sampleFromListWithoutReplacement(java.util.List theList, DoubleUniform theRnd) {
        int n = theRnd.nextIntFromTo(0,theList.size()-1);
        Object o = theList.get(n);
        theList.remove(o);
        return o;
    }

    public static Object sampleFromListWithReplacement(java.util.List theList, DoubleUniform theRnd) {
        int n = theRnd.nextIntFromTo(0,theList.size()-1);
        Object o = theList.get(n);
        return o;
    }
    

}
