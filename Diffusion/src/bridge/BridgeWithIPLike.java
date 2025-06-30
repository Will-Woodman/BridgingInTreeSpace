/*
BridgeWithIPLike
    Copyright (C) 2025  William M Woodman

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

    Contact the author at:  <w.m.woodman2@ncl.ac.uk>
                           
 */

package bridge;

import geodesics.Geodesic;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;

/** Bridge where likelihood is the independence proposal density. Used for simulating
 * bridges from the independence proposal given that all steps are simple using MCMC.
 * This forms part of the reference distribution for Chib one block marginal likelihood estimator (for unknown t0)
 */

public class BridgeWithIPLike extends ForwardStepBridge {
    
    public BridgeWithIPLike(TreeAsSplits tA, TreeAsSplits tB, int n) {
        super(tA, tB, n);
    }
    
    @Override
    /** Make empty copy for use in MCMC proposal */
    public ForwardStepBridge makeEmptyCopy() {
        return new BridgeWithIPLike(treeA, treeB, numSteps);
    }
    @Override
    public ForwardStepBridge makeFreshBridge(TreeAsSplits tA, TreeAsSplits tB) {
        return new BridgeWithIPLike(tA, tB, numSteps);
    }
    @Override
    public ForwardStepBridge makeEmptyCopy(TreeAsSplits x, boolean whichEnd) {
        if (whichEnd) return new BridgeWithIPLike(treeA, x, numSteps);
        else return new BridgeWithIPLike(x, treeB, numSteps);
    }

    /** Compute the log likelihood associated with a particular step. 
 */
    protected double[] calcStepLogLike(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException {
        //First we calculate the number of steps remaining on the bridge which informs the next steps
        int k=1;
        int j=1;
        while(k<numSteps+1){
            if(tB.equals(theTrees[k])){
               j=k;
               k=numSteps+1;
               
            }
            else{
                k=k+1;
            }
        }
        //calculate the log likelihood which depends on j:
        double PropLogLike = calcStepLogLikeProp(tA, tB, g, t0,j)[0];
        
        //System.out.println(PostLogLike+" "+PropLogLike+" "+ll);
        return new double[] {PropLogLike, 0.0};
        
    }
    
    protected double[] calcStepLogLikeProp(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0,int j) throws AlgorithmException {
        int nPrime = tA.getNumTaxa()-3;
        double var=t0/numSteps;
        double RWsd = Math.sqrt(var);
        double sd = Math.sqrt(var*((double)(numSteps-j))/((double)(numSteps-j+1)));
        int numStepsRemaining =numSteps+1-j;
        double ll;
        if(numStepsRemaining==1) ll=0;
        //here the mixture step density is normalised because the other part of the mixture density is also normalised
        else ll=computeFixedMixtureStepDensity(tA,  theTrees[numSteps], tB, sd,RWsd, numStepsRemaining, unifDist)-nPrime*Math.log(Math.sqrt(2*Math.PI));

      
        
        return new double[] {ll, 0.0};
        
        
    }
    
  
}
