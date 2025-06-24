/** 
    F84Model
    Implementation of the F84 model for DNA substitution
  
    Copyright (C) 2011  Tom M. W. Nye

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

package sequencemodels;

/**
 * Implementation of the F84 model for DNA substitution
 */

public class F84Model extends TN93Model {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /** Constructor: fix rate matrix etc.
     p[] must be an array length 4 containing the stationary distribution in order A, G, C, T*/
    public F84Model(double[] p, double kappa) {
        super(p, (1.0+kappa/(p[2]+p[3])), (1.0+kappa/(p[0]+p[1])));
    }

    public void changeParameters(double[] p, double kappa) {
        setParameters(p, (1.0+kappa/(p[2]+p[3])), (1.0+kappa/(p[0]+p[1])));
    }

    /** Formulae from Ziheng Yang, page 16 */
    @Override
    public double calcDistance(int[] a, int[] b, double infiniteDist) throws treebase.AlgorithmError {
        if (a.length!=b.length) throw new treebase.AlgorithmError("Sequence length mismatch in F84 distance calculation");
        // Count up where a and b match
        double[] res = countTransitionsTransversions(a, b);
        double v = res[0];
        double s = res[1]+res[2];

        double x;

        x = 1.0-0.5*s/(piT*piC/piY+piA*piG/piR);
        x -= 0.5*v*(piT*piC*piR/piY+piA*piG*piY/piR)/(piT*piC*piR+piA*piG*piY);
        if (x<=0.0) return infiniteDist;
        double aa = -Math.log(x);

        x = 1.0-0.5*v/(piY*piR);
        if (x<=0.0) return infiniteDist;
        double bb = -Math.log(x);

        double d = 2.0*(piT*piC/piY+piA*piG/piR)*aa;
        d -= 2.0*(piT*piC/piY*piR+piA*piG/piR*piY-piY*piR)*bb;

        return d;
    }


}
