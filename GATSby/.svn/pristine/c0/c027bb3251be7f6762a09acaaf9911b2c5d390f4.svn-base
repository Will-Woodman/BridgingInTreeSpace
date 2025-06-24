/**
 * Optimisation1D
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
                            * */

package geodesics;

/**
 * Methods for minimising a 1D function
 */

import treebase.AlgorithmException;

public class Optimisation1D {
    
    protected static final double GOLDEN = 1.618;

    public interface Objective {
        public double evaluate(double x) throws AlgorithmException;
    }
    
    public static double[] goldenRatioOptimise(Objective theObjective, double[] initialValues, double[] objectiveValues, double tolerance) throws AlgorithmException {
        
        double x1, x2, snew, ynew;
        double s1 = initialValues[0]; 
        double s2 = initialValues[2];
        double s = initialValues[1];
        double y1 = objectiveValues[0];
        double y2 = objectiveValues[2];
        double y = objectiveValues[1];
        
        boolean converged = false;
        while (!converged) {
//                System.out.println(s1+" "+s+" "+s2+" "+y1+" "+y+" "+y2);

            x1 =s-s1;
            x2 = s2-s;
            if (x2>x1) {
                // Upper interval is wider
                snew = (GOLDEN*s+s2)/(1+GOLDEN);
                ynew = theObjective.evaluate(snew);
                if (ynew<y) {
                    // New interval is mid - new - upper
                    s1 = s; y1 = y; 
                    s = snew; y = ynew; 
                }
                else {
                    // New interval is lower - mid - new
                    s2 = snew; y2 = ynew; 
                }
            }
            else {
                // Lower interval is wider
                snew = (GOLDEN*s+s1)/(1+GOLDEN);
                ynew = theObjective.evaluate(snew);

                if (ynew<y) {
                    // New interval is lower - new - mid
                    s2 = s; y2 = y; 
                    s = snew; y = ynew; 
                }
                else {
                    // New interval is new - mid - upper
                    s1 = snew; y1 = ynew; 
                }
            }

            // Test for convergence
            if (Math.abs(s2-s1)<tolerance) converged = true;

        } // End golden ratio search

        return new double[] {s,y};
    }

    
    /** Find an initial bracket and optimise */
    public static double[] optimise(Objective theObjective, double[] bounds, double tolerance) throws AlgorithmException {
        
        double x1 = bounds[0];
        double x2 = bounds[1];
        double x = x1;
        
        double y1 = theObjective.evaluate(x1);
        double y2 = theObjective.evaluate(x2);
        double y = y1;
        
        boolean gotBracket = false;
        
        while (!gotBracket) {
            
            x = (x2+GOLDEN*x1)/(1.0+GOLDEN);
            y = theObjective.evaluate(x);
            
            if ((y1>y)&&(y2>y)) {
                gotBracket = true;
            }
            else {
                
                // Update interval
                if (y1<y2) {
                    x2 = x;
                    y2 = y;
                }
                else {
                    x1 = x;
                    y1 = y;
                }

                // Check for convergence
                if (Math.abs(x2-x1)<tolerance) {
                    // Failed to find a bracket so return min value so far!
                    return new double[] {x,y};
                }
            }
        }
        
        // Now we have a bracket, run golden ratio
        return goldenRatioOptimise(theObjective, new double[] {x1, x, x2}, new double[] {y1, y, y2}, tolerance);
    }


    /** Optimise given an initial guess at the point */
    public static double[] optimise(Objective theObjective, double initVal, double initObj, double[] bounds, double tolerance) throws AlgorithmException {

        // Check value within bds
        if ((initVal<bounds[0])||(initVal>bounds[1])) {
            throw new AlgorithmException("Initial value for 1D optimisation lies outside given bounds.");
        }

        // Find largest interval
        double x, x1, x2;
        if ((bounds[1]-initVal)>(initVal-bounds[0])) {
            // Lower interval is smaller
            x1 = bounds[0]; x = initVal;
            x2 = x1 + GOLDEN*(x-x1);
        }
        else {
            // Upper interval is smaller
            x2 = bounds[1]; x = initVal;
            x1 = x2 - GOLDEN*(x2-x);
        }

        double y1 = theObjective.evaluate(x1);
        double y2 = theObjective.evaluate(x2);
        double y = initObj;

        if ((y1>y)&&(y2>y)) {
            // Got a bracket
            return goldenRatioOptimise(theObjective, new double[] {x1, x, x2}, new double[] {y1, y, y2}, tolerance);
        }

        // If no bracket, search on smallest interval
        if (y1>y2) {
            return optimise(theObjective, new double[] {x, x2}, tolerance);
        }
        else {
            return optimise(theObjective, new double[] {x1, x}, tolerance);
        }
    }
}