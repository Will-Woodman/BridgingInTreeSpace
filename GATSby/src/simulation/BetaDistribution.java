/**
 *  BetaDistribution
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
 * Beta distribution.
 *
 */

import cern.jet.random.tdouble.Beta;
import treebase.AlgorithmException;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class BetaDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;
    
    /* Instance variables */
    private double shape1,shape2;
    private Beta coltObject;
    
    public BetaDistribution(double alpha, double beta) {
        // Set up instance variables
        shape1 = alpha;
        shape2 = beta;
        coltObject = new Beta(shape1, shape2, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public BetaDistribution(double alpha, double beta, int seed) {
        // Set up instance variables
        shape1 = alpha;
        shape2 = beta;
        coltObject = new Beta(shape1, shape2, new DoubleMersenneTwister(seed));
    }
    
    public BetaDistribution(double alpha, double beta, DoubleMersenneTwister randomEngine) {
        // Set up instance variables
        shape1 = alpha;
        shape2 = beta;
        coltObject = new Beta(shape1, shape2, randomEngine);
    }

    public void resetRandomEngineSeed() {
        coltObject = new Beta(shape1, shape2, Random.getEngine());
    }

    public void setParams(double alpha, double beta)  {
        shape1 = alpha;
        shape2 = beta;
        coltObject.setState(shape1, shape2);
    }
    
    public double getShape1(){
        return shape1;
    }
    
    public double getShape2(){
        return shape2;
    }
    
    /** log probability density function of the Beta distribution */
    public double logpdf(double x)
    {
        return cern.jet.stat.tdouble.Gamma.logGamma(shape1+shape2) - cern.jet.stat.tdouble.Gamma.logGamma(shape1) - cern.jet.stat.tdouble.Gamma.logGamma(shape2) + (shape1 - 1.0)*Math.log(x) + (shape2 - 1.0)*Math.log(1.0 -x);
    }

    /** probability density function of the Beta distribution */
    public double pdf(double x)
    {
        return Math.exp(this.logpdf(x));
    }
    
    /** cumulative distribution function of the Beta distribution */
    public double cdf(double x)
    {
        return coltObject.cdf(x);
    }
    
    /** log probability density function of the Beta distribution */
    public static double logpdf(double x, double alpha, double beta)
    {
        return cern.jet.stat.tdouble.Gamma.logGamma(alpha+beta) - cern.jet.stat.tdouble.Gamma.logGamma(alpha) - cern.jet.stat.tdouble.Gamma.logGamma(beta) + (alpha - 1.0)*Math.log(x) + (beta - 1.0)*Math.log(1.0 -x);
    }

    /** probability density function of the Beta distribution */
    public static double pdf(double x, double alpha, double beta)
    {
        return Math.exp(logpdf(x,alpha,beta));
    }
    
    /** cumulative density function of the Beta distribution */
    public static double cdf(double x, double alpha, double beta)
    {
        return cern.jet.stat.tdouble.Probability.beta(alpha,beta,x);
    }

    /** Random number generation */
    public double sample() {
        return coltObject.nextDouble();
    }
    
    /** Random number generation: bypass internal state values */
    public double sample(double alpha, double beta) {
        return coltObject.nextDouble(alpha, beta);
    }
    
    /** Quantile calculation: copied from the C code on page 273
        of Numerical Recipes: The Art of Scientific Computing, Third Edition
        by William H Press, Saul E Teukolsky, William T Vetterling and Brian
        P Flannery, Cambridge University Press*/
    public double quantile(double p) {
        return quantile(p, shape1, shape2);
    }
    
    public static double quantile(double p, double alpha, double beta) {
        double eps = 1e-8;
        double pp, t, u, err, x, al, h, w, afac;
        double a1=alpha-1, b1=beta-1;
        if (p <= 0) return 0;
        else if (p >= 1) return 1;
        else if (alpha >= 1 && beta >= 1) {
            //Set initial guess.
            pp = (p < 0.5)? p : 1 - p;
            t = Math.sqrt(-2*Math.log(pp));
            x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
            if (p < 0.5) x = -x;
            al = (x*x-3)/6;
            h = 2/(1/(2*alpha-1.)+1/(2*beta-1));
            w = (x*Math.sqrt(al+h)/h)-(1/(2*beta-1)-1/(2*alpha-1))*(al+5/6-2/(3*h));
            x = alpha/(alpha+beta*Math.exp(2*w));
        } else {
            double lna = Math.log(alpha/(alpha+beta));
            double lnb = Math.log(beta/(alpha+beta));
            t = Math.exp(alpha*lna)/alpha;
            u = Math.exp(beta*lnb)/beta;
            w = t + u;
            if (p < t/w) x = Math.pow(alpha*w*p,1/alpha);
            else x = 1 - Math.pow(beta*w*(1-p),1/beta);
        }
        afac = -cern.jet.stat.tdouble.Gamma.logGamma(alpha) -
                   cern.jet.stat.tdouble.Gamma.logGamma(beta) +
                   cern.jet.stat.tdouble.Gamma.logGamma(alpha+beta);
        for (int j=0;j<10;j++) {
            if (x == 0 || x == 1) return x;
            //a or b too small for accurate calculation.
            err = cdf(x, alpha, beta) - p;
            t = Math.exp(a1*Math.log(x)+b1*Math.log(1-x) + afac);
            u = err/t;
            //Halley:
            x -= (t = u/(1-0.5*Math.min(1,u*(a1/x - b1/(1-x)))));
            if (x <= 0) x = 0.5*(x + t);
            //Bisect if x tries to go neg or > 1.
            if (x >= 1) x = 0.5*(x + t + 1);
            if (Math.abs(t) < eps*x && j > 0) break;
        }
        return x;
    }
    
    
}
