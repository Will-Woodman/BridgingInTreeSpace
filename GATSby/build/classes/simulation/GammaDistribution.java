/**
 * GammaDistribution
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
 * Gamma distribution.
 * Version of colt gamma distrib.
 * Quantiles adapted from dynare implementation
 *
 */

import cern.jet.random.tdouble.Gamma;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class GammaDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Data */
    private double shape, scale;
    private Gamma coltObject;

    public GammaDistribution(double alpha, double theta) {
        shape = alpha;
        scale = theta;
        coltObject = new Gamma(shape, 1.0/scale, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public GammaDistribution(double alpha, double theta, int seed) {
        shape = alpha;
        scale = theta;
        coltObject = new Gamma(shape, 1.0/scale, new DoubleMersenneTwister(seed));
    }
    public GammaDistribution(double alpha, double theta, DoubleMersenneTwister randomEngine) {
        shape = alpha;
        scale = theta;
        coltObject = new Gamma(shape, 1.0/scale, randomEngine);
    }
    
    public void resetRandomEngineSeed() {
        coltObject = new Gamma(shape, 1.0/scale, Random.getEngine());
    }
    
    public void setParams(double alpha, double theta)  {
        shape = alpha;
        scale = theta;
        coltObject.setState(alpha, 1.0/theta);
    }

    public double getShape() {
        return shape;
    }

    public double getScale() {
        return scale;
    }

    /** log probability density function of the Gamma distribution */
    public double logpdf(double x)
    {
        return (shape-1.0)*Math.log(x) - x/scale - shape*Math.log(scale) - cern.jet.stat.tdouble.Gamma.logGamma(shape);
    }
    /** probability density function of the Gamma distribution */
    public double pdf(double x)
    {
        return Math.exp(this.logpdf(x));
    }
    /** cumulative distribution function of the Gamma distribution */
    public  double cdf(double x)
    {
        return cern.jet.stat.tdouble.Probability.gamma(1.0/scale,shape,x);
    }
    /** log probability density function of the Gamma distribution */
    public static double logpdf(double x, double shape, double scale)
    {
        return (shape-1.0)*Math.log(x) - x/scale - shape*Math.log(scale) - cern.jet.stat.tdouble.Gamma.logGamma(shape);
    }
    /** probability density function of the Gamma distribution */
    public static double pdf(double x, double shape, double scale)
    {
        return Math.exp(logpdf(x,shape,scale));
    }
    /** cumulative distribution function of the Gamma distribution */
    public static double cdf(double x, double shape, double scale)
    {
        return cern.jet.stat.tdouble.Probability.gamma(1.0/scale, shape, x); // NB: strange order of arguments due to colt weirdness!
    }

    /** Random number generation */
    public double sample() {
        return coltObject.nextDouble();
    }
    /** Random number generation */
    public double sample(double shape, double scale) {
        return coltObject.nextDouble(shape, 1.0/scale);
    }

    /** Calculate a quantile */
    public double quantile(double prob)
    {
        return quantile(prob, shape, scale);
    }

    /** Quantile calculation -- taken from Geogebra UNSTABLE!!! */
    public static double quantileOLD(double p, double alpha, double scale)
    {
        double a, b, c, ch, g, p1, v;
        double p2, q, s1, s2, s3, s4, s5, s6, t=0.0, x;
        int i;
        double aa = 0.6931471805;

        /* test arguments and initialise */

    /*!* #ifdef IEEE_754 /*4!*/
        if (Double.isNaN(p) || Double.isNaN(alpha) || Double.isNaN(scale))
        return p + alpha + scale;
    /*!* #endif /*4!*/

        if (p < 0 || p > 1 || alpha <= 0) {
        throw new java.lang.ArithmeticException("Math Error: DOMAIN");
        //      return Double.NaN;
        }
        if (/* 0 <= */ p < 0.000002) return 0;
        if (/* 1 >= */ p > 0.999998) return Double.POSITIVE_INFINITY;

        v = 2*alpha;

        c = alpha-1;
    /*!*     g = lgammafn(alpha);!!!COMMENT!!! *!*/
        g = cern.jet.stat.tdouble.Gamma.logGamma(alpha);/* log Gamma(v/2) */

    /*!*     if(v < (-1.24)*log(p)) { *!*/
        if(v < (-1.24)*java.lang.Math.log(p)) {
          /* starting approximation for small chi-squared */

    /*!*    ch = pow(p*alpha*exp(g+alpha*Constants.M_LN_2), 1/alpha); *!*/
        ch = java.lang.Math.pow(p*alpha*java.lang.Math.exp(g+alpha*aa), 1/alpha);
        if(ch < 5e-7) {
            throw new java.lang.ArithmeticException("Math Error: DOMAIN");
            //              return Double.NaN;
        }

        } else if(v > 0.32) {

        /* starting approximation using Wilson and Hilferty estimate */

        x = NormalDistribution.quantile(p, 0, 1);
        p1 = 0.222222/v;
    /*!*    ch = v*pow(x*sqrt(p1)+1-p1, 3); *!*/
        ch = v*java.lang.Math.pow(x*java.lang.Math.sqrt(p1)+1-p1, 3);

        /* starting approximation for p tending to 1 */

        if( ch > 2.2*v + 6 )
    /*!*        ch = -2*(log(1-p) - c*log(0.5*ch) + g); *!*/
            ch = -2*(java.lang.Math.log(1-p) - c*java.lang.Math.log(0.5*ch) + g);

        } else { /* starting approximation for v <= 0.32 */

        ch = 0.4;
    /*!*    a = log(1-p) + g + c*Constants.M_LN_2; *!*/
        a = java.lang.Math.log(1-p) + g + c*aa;
        do {
            q = ch;
            p1 = 1+ch*(4.67+ch);
            p2 = ch*(6.73+ch*(6.66+ch));
            t = -0.5 +(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
    /*!*        ch -= (1- exp(a+0.5*ch)*p2/p1)/t; *!*/
            ch -= (1- java.lang.Math.exp(a+0.5*ch)*p2/p1)/t;
    /*!*    } while(fabs(q/ch - 1) > EPS1); *!*/
        } while(java.lang.Math.abs(q/ch - 1) > 1e-2);
        }

        /* algorithm AS 239 and calculation of seven term taylor series */

        for( i=1 ; i <= 20 ; i++ ) {
        q = ch;
        p1 = 0.5*ch;
        p2 = p - cern.jet.stat.tdouble.Probability.gamma(p1, alpha, 1.0);
    /*!* #ifdef IEEE_754 /*4!*/
        if(Double.isInfinite(p2)) return Double.NaN;

    /*!*    t = p2*exp(alpha*Constants.M_LN_2+g+p1-c*log(ch)); *!*/
        t = p2*java.lang.Math.exp(alpha*aa+g+p1-c*java.lang.Math.log(ch));
        b = t/ch;
        a = 0.5*t-b*c;
        s1 = (210+a*(140+a*(105+a*(84+a*(70+60*a)))))/420;
        s2 = (420+a*(735+a*(966+a*(1141+1278*a))))/2520;
        s3 = (210+a*(462+a*(707+932*a)))/2520;
        s4 = (252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
        s5 = (84+2264*a+c*(1175+606*a))/2520;
        s6 = (120+c*(346+127*c))/5040;
        ch = ch+t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    /*!*    if(fabs(q/ch-1) > EPS2) *!*/
        if(java.lang.Math.abs(q/ch-1) > 5e-7)
            return 0.5*scale*ch;
        }
        throw new java.lang.ArithmeticException("Math Error: PRECISION");
        //        return 0.5*scale*ch;
    }
    
    /** Quantile calculation -- taken from Dynare:
        http://www.dynare.org/dynare-matlab-m2html/matlab/missing/stats/gaminv.html
    */
    public static double quantileCurrent(double p, double alpha, double scale) {
        // Check arguments
        if (Double.isNaN(p) || Double.isNaN(alpha) || Double.isNaN(scale)) return p + alpha + scale;
        
        if (p < 0 || p > 1 || alpha <= 0 || scale <=0) throw new java.lang.ArithmeticException("Math Error: DOMAIN");
        
        // Deal with trivial cases
        if (p == 1) {
            return Double.POSITIVE_INFINITY;
        }
        if(p ==0) {
            return 0.0;
        }
        // Deal with non-trivial cases
        double eps = Math.pow(2.0, -52.0);
        double y = alpha * scale;
        if (p < eps) {
            y = Math.sqrt(y);
        }
        double yOld = y;
        double yNew = 0.0;
        for (int i = 1; i <= 100; i++) {
            double h = (cdf(yOld, alpha, scale) - p) / pdf(yOld, alpha, scale);
            yNew = yOld - h;
            if (yNew <= eps) {
                yNew = yOld / 10;
                h = yOld - yNew;
            }
            if (Math.abs(h) < Math.sqrt(eps)) {
                break;
            }
            yOld = yNew;
        }
        return yNew;
    }
    
    public static double quantile(double p, double alpha, double scale) {
        return scale*qchisq(p, 2.0*alpha)/2.0;
    }
    
    public static double qchisq(double prob, double nu) {
        double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6, v = nu;
        
        if (p < 0.000002 || p > 0.999998 || v <=0.0) return -1.0; // Not sure about this
        
        g = cern.jet.stat.tdouble.Gamma.logGamma(v/2.0);
	xx = v/2.0;   
	c = xx - 1.0;
	if (v < -1.24*Math.log(p)) {
	    ch = Math.pow((p*xx*Math.exp(g+xx*aa)), 1.0/xx);
	    if (ch-e<0) 
		return (ch);
        } else {
            // goto l1;
            if (v <= 0.32) {
		ch = 0.4;   
		a = Math.log(1.0-p);
                // goto l2
                while(true) {
		    q = ch;  
		    p1 = 1.0+ch*(4.67+ch);  
		    p2 = ch*(6.73+ch*(6.66+ch));
		    t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		    ch -= (1.0-Math.exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		    if (Math.abs(q/ch-1.0)-0.01 <= 0.0) 
                        break;
                }
            } else {
                // goto l3
                x = NormalDistribution.quantile(p, 0.0, 1.0);
                p1 = 0.222222/v;   
	        ch = v*Math.pow((x*Math.sqrt(p1)+1.0-p1), 3.0);
	        if (ch > 2.2*v+6.0)  
	            ch = -2.0*(Math.log(1.0-p)-c*Math.log(0.5*ch)+g);
            }
        }
        // goto l4
        while(true) {
            q = ch;   
	    p1 = 0.5*ch;
	    if ((t = simulation.GammaDistribution.cdf (p1, xx, 1.0)) < 0.0) {
	        System.out.println("Error: Problem in qchisq");
	        return -1.0;
	    }
	    p2 = p-t;
	    t = p2*Math.exp(xx*aa+g+p1-c*Math.log(ch));   
	    b = t/ch;  
	    a = 0.5*t-b*c;
	    s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
	    s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
	    s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
	    s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
	    s5 = (84.0+264.0*a+c*(175.0+606.0*a)) / 2520.0;
	    s6 = (120.0+c*(346.0+127.0*c)) / 5040.0;
	    ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	    if (Math.abs(q/ch-1.0) <= e) 
	        break;
        }
        return ch;

	/*l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;   
		a = Math.log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-Math.exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (Math.abs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = NormalDistribution.quantile(p, 0.0, 1.0);
		p1 = 0.222222/v;   
		ch = v*Math.pow((x*Math.sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch = -2.0*(Math.log(1.0-p)-c*Math.log(0.5*ch)+g);
	l4:
		q = ch;   
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0) {
	            System.out.println("Error: Problem in qchisq");
	            return -1.0;
		}
		p2 = p-t;
		t = p2*Math.exp(xx*aa+g+p1-c*Math.log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a)) / 2520.0;
		s6 = (120.0+c*(346.0+127.0*c)) / 5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (Math.abs(q/ch-1.0) > e) 
			goto l4;
                
		return ch;*/
        
    }


}
