/**
 *  DirichletDistribution
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
 * Dirichlet distribution.
 *
 */

import java.io.*;
import cern.jet.random.tdouble.Gamma;
import treebase.AlgorithmException;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class DirichletDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Instance variables */
    private double[] mean;
    private double concentration;
    private Gamma gammaDistrib;

    /** Constructor: assumes a probability vector (sum 1) */
    public DirichletDistribution(double[] probVector, double concParameter) throws AlgorithmException {
         //Error checking:
        if(concParameter <= 0.0 || !isSimplex(probVector))
            throw new AlgorithmException("probVector does not lie in the simplex"
                    + " or concentration parameter <= 0 in DirichletDistribution.");
        
        int dimension = probVector.length;
        mean = new double[dimension];
   
        // Copy values from the input param. 
        System.arraycopy(probVector,0,mean,0,dimension);
        concentration = concParameter;
        //A Gamma object is needed to generate samples from the Dirichlet distribution:
        gammaDistrib = new Gamma(1.0, 1.0, Random.getEngine());
        
    }
    
    /** Constructor from general vector */
    public DirichletDistribution(double[] argVector) throws AlgorithmException {
        int i, dimension = argVector.length;
        
        //Error checking:
        for(i=0;i<dimension;i++){
            if(argVector[i] <= 0)
                throw new AlgorithmException("argVector has entries <= 0"
                    + " in DirichletDistribution.");
        }
        
        mean = new double[dimension];
        concentration=0.0;
        
        for(i=0;i<dimension;i++)
            concentration+=argVector[i];
        for(i=0;i<dimension;i++)
            mean[i]=argVector[i]/concentration;
        gammaDistrib = new Gamma(1.0, 1.0, Random.getEngine());

    }
    
    /** Default constructor: creates Dirichlet(1,1,...,1) */
    public DirichletDistribution(int n) {
        int i;
        
        mean = new double[n];
        concentration=1.0;

        for(i=0;i<n;i++)
            mean[i]=1.0/n;
        gammaDistrib = new Gamma(1.0, 1.0, Random.getEngine());

    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public DirichletDistribution(double[] probVector, double concParameter, int seed) throws AlgorithmException {
         //Error checking:
        if(concParameter <= 0.0 || !isSimplex(probVector))
            throw new AlgorithmException("probVector does not lie in the simplex"
                    + " or concentration parameter <= 0 in DirichletDistribution.");
        
        int dimension = probVector.length;
        mean = new double[dimension];
   
        // Copy values from the input param. 
        System.arraycopy(probVector,0,mean,0,dimension);
        concentration = concParameter;
        //A Gamma object is needed to generate samples from the Dirichlet distribution:
        gammaDistrib = new Gamma(1.0, 1.0, new DoubleMersenneTwister(seed));
        
    }
    public DirichletDistribution(double[] probVector, double concParameter, DoubleMersenneTwister randomEngine) throws AlgorithmException {
         //Error checking:
        if(concParameter <= 0.0 || !isSimplex(probVector))
            throw new AlgorithmException("probVector does not lie in the simplex"
                    + " or concentration parameter <= 0 in DirichletDistribution.");
        
        int dimension = probVector.length;
        mean = new double[dimension];
   
        // Copy values from the input param. 
        System.arraycopy(probVector,0,mean,0,dimension);
        concentration = concParameter;
        //A Gamma object is needed to generate samples from the Dirichlet distribution:
        gammaDistrib = new Gamma(1.0, 1.0, randomEngine);
        
    }
    public DirichletDistribution(double[] argVector, int seed) throws AlgorithmException {
        int i, dimension = argVector.length;
        
        //Error checking:
        for(i=0;i<dimension;i++){
            if(argVector[i] <= 0)
                throw new AlgorithmException("argVector has entries <= 0"
                    + " in DirichletDistribution.");
        }
        
        mean = new double[dimension];
        concentration=0.0;
        
        for(i=0;i<dimension;i++)
            concentration+=argVector[i];
        for(i=0;i<dimension;i++)
            mean[i]=argVector[i]/concentration;
        gammaDistrib = new Gamma(1.0, 1.0, new DoubleMersenneTwister(seed));

    }
    public DirichletDistribution(double[] argVector, DoubleMersenneTwister randomEngine) throws AlgorithmException {
        int i, dimension = argVector.length;
        
        //Error checking:
        for(i=0;i<dimension;i++){
            if(argVector[i] <= 0)
                throw new AlgorithmException("argVector has entries <= 0"
                    + " in DirichletDistribution.");
        }
        
        mean = new double[dimension];
        concentration=0.0;
        
        for(i=0;i<dimension;i++)
            concentration+=argVector[i];
        for(i=0;i<dimension;i++)
            mean[i]=argVector[i]/concentration;
        gammaDistrib = new Gamma(1.0, 1.0, randomEngine);

    }
    public DirichletDistribution(int n, int seed) {
        int i;
        
        mean = new double[n];
        concentration=1.0;

        for(i=0;i<n;i++)
            mean[i]=1.0/n;
        gammaDistrib = new Gamma(1.0, 1.0, new DoubleMersenneTwister(seed));

    }
    public DirichletDistribution(int n, DoubleMersenneTwister randomEngine) {
        int i;
        
        mean = new double[n];
        concentration=1.0;

        for(i=0;i<n;i++)
            mean[i]=1.0/n;
        gammaDistrib = new Gamma(1.0, 1.0, randomEngine);

    }
    
    public void resetRandomEngineSeed() {
        gammaDistrib = new Gamma(1.0, 1.0, Random.getEngine());
    }
    
    /** Test whether a vector lies on the simplex (sums to 1 and positive entries) */
    private static boolean isSimplex(double[] x){
        double sum = 0.0;
        final double TOLERANCE = 1.0E-7;
        
        for(int i=0;i<x.length;i++){
            if(x[i] <= 0.0)
                return false;
            sum+=x[i];
        }
        if(Math.abs(sum-1.0) > TOLERANCE)
            return false;
        
        return true;
    }

    /** Set the mean (probability vector) */
    public void setMean(double[] probVector) throws AlgorithmException {
        int dimension = mean.length;
        
        //Error checking
        if(dimension!=probVector.length)
            throw new AlgorithmException("probVector does not have the"
                    + " correct dimension in DirichletDistribution.setMean.");
        
        if(!isSimplex(probVector))
            throw new AlgorithmException("probVector does not lie in the simplex"
                    + " or concentration parameter <= 0 in DirichletDistribution.setMean.");
 
        System.arraycopy(probVector,0,mean,0,dimension);
        
    }
    
    public void setConcentration(double concParameter) throws AlgorithmException {
        //Error checking:
        if(concParameter <= 0.0) throw new AlgorithmException("concentration parameter <= 0 in DirichletDistribution.setConcentration.");
        
        concentration = concParameter;
    }

    /** Provide access to parameters */
    public double getConcParam() {
        return concentration;
    }
    public double getConcParam(int i) {
        return concentration*mean[i];
    }
    public double[] getConcParamVector() {
        double[] ret = new double[mean.length];
        for(int i=0; i<mean.length; i++) {
            ret[i] = concentration*mean[i];
        }
        return ret;
    }
    public double[] getMean() {
        double[] m = new double[mean.length];
        System.arraycopy(mean, 0, m, 0, mean.length);
        return m;
    }

    /** log probability density function of the Dirichlet distribution */
    public double logpdf(double[] x) {
        //Error checking
        /* if(!isSimplex(x))
            throw new AlgorithmException("x does not lie in the simplex"
                    + " in DirichletDistribution.lnpdf/pdf.");
        */
        double sum=cern.jet.stat.tdouble.Gamma.logGamma(concentration);
        
        for(int i=0;i<mean.length;i++)
            sum+=(concentration*mean[i]-1)*Math.log(x[i])
                    -cern.jet.stat.tdouble.Gamma.logGamma(concentration*mean[i]);
        
        return sum;
    }
    
    /** probability density function of the Dirichlet distribution */
    public double pdf(double[] x) {
        return Math.exp(this.logpdf(x));
    }
    
    /** log probability density function of the Dirichlet distribution */
    public static double logpdf(double[] x, double[] probVector, double concParameter) {
        int dimension=x.length;
        
        //Error checking:
        //if(dimension!=probVector.length)  throw new AlgorithmException("x and probVector have different dimensions"
        //            + " in DirichletDistribution.lnpdf/pdf.");
        /*
        if(concParameter <= 0.0 || !isSimplex(x) || !isSimplex(probVector))
            throw new AlgorithmException("x or probVector does not lie in the simplex"
                    + " or concentration parameter <= 0 in DirichletDistribution.lnpdf/pdf.");
        */
        double sum=cern.jet.stat.tdouble.Gamma.logGamma(concParameter);
        
        for(int i=0;i<dimension;i++)
            sum+=(concParameter*probVector[i]-1)*Math.log(x[i])
                    -cern.jet.stat.tdouble.Gamma.logGamma(concParameter*probVector[i]);
        
        return sum;

    }
    
    /** probability density function of the Dirichlet distribution */
    public static double pdf(double[] x, double[] probVector, double concParameter) {
        return Math.exp(logpdf(x,probVector,concParameter));
    }

    /** Random number generation */
    public double[] sample() {
        int i,k=mean.length;
        double[] draw = new double[k];
        double sum=0.0;
        
        for(i=0;i<k;i++){
            draw[i]=gammaDistrib.nextDouble(mean[i]*concentration,1.0);
            sum+=draw[i];
        }
        for(i=0;i<k;i++){
            draw[i]=draw[i]/sum;
        }
        return draw;
    }
    /** Generate random sample into given existing vector */
    public void sample(double[] draw) {
        int i,k=mean.length;
        double sum=0.0;
        
        for(i=0;i<k;i++){
            draw[i]=gammaDistrib.nextDouble(mean[i]*concentration,1.0);
            sum+=draw[i];
        }
        for(i=0;i<k;i++){
            draw[i]=draw[i]/sum;
        }
    }

        
    /** Random number generation */
    public double[] sample(double[] probVector, double concParameter) {
        //Error checking:
        //if(concParameter <= 0.0 || !isSimplex(probVector)) throw new AlgorithmException("probVector does not lie in the simplex"
        //           + " or concentration parameter <= 0 in DirichletDistribution.sample.");
        
        int i,k=probVector.length;
        double[] draw = new double[k];
        double sum=0.0;
        
        for(i=0;i<k;i++){
            draw[i]=gammaDistrib.nextDouble(probVector[i]*concParameter,1.0);
            sum+=draw[i];
        }
        for(i=0;i<k;i++){
            draw[i]=draw[i]/sum;
        }
        return draw;
    }
    


}