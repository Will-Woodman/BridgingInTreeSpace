/*
BridgeWithApproxMVNLike
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

package bridge;

import cern.jet.random.tdouble.DoubleUniform;
import geodesics.Geodesic;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import simulation.NormalDistribution;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.TreeAsSplits;

/** Implementation of the forward bridge using approx MVN log like obtained from a 
    simple geodesic at each step.
 */

public class BridgeWithStepStoneLike extends ForwardStepBridge {
    protected double beta;
    
    public BridgeWithStepStoneLike(TreeAsSplits tA, TreeAsSplits tB, int n, double Beta) {
        super(tA, tB, n);
        beta=Beta;
    }
    
    //constructor from array of trees used for quasistatic stepping stone sampling
    public BridgeWithStepStoneLike(TreeAsSplits[] startTrees, double Beta) throws AlgorithmError {
        super(startTrees);
        beta=Beta;
    }
    
    @Override
    /** Make empty copy for use in MCMC proposal */
    public ForwardStepBridge makeEmptyCopy() {
        return new BridgeWithStepStoneLike(treeA, treeB, numSteps,beta);
    }
    @Override
    public ForwardStepBridge makeFreshBridge(TreeAsSplits tA, TreeAsSplits tB) {
        return new BridgeWithStepStoneLike(tA, tB, numSteps,beta);
    }
    @Override
    public ForwardStepBridge makeEmptyCopy(TreeAsSplits x, boolean whichEnd) {
        if (whichEnd) return new BridgeWithStepStoneLike(treeA, x, numSteps,beta);
        else return new BridgeWithStepStoneLike(x, treeB, numSteps,beta);
    }
  
    //used to set the new log likes when switching between values of beta in the quasistatic
    //stepping stone algorithm
     public void setLogLikes(double t0) {
         totalLogLike = 0;
        for (int k=1; k<=numSteps; k++) {
            try {
                Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
                double[] step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0);
                logLike[k] = step[0];
                totalLogLike+=step[0];
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error updating a bridge for new t0 value. This should not happen: the step should be simple. \nError report --"+anErr.getMessage());
            }
             
        }
       
    }
    

    /** Compute the log likelihood associated with a particular step. 
     in a similar manner to MVN step. In this case must calculate the indep proposal likelihood and the GGF likelihood*/
    protected double[] calcStepLogLike(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException {
        //need to think of a better way of doing the below to get the number of steps remaining
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
        
        double PostLogLike = calcStepLogLikePost(tA, tB, g, t0)[0];
        double PropLogLike = calcStepLogLikeProp(tA, tB, g, t0,j)[0];
        
        double ll=0;
        
        if(j==numSteps){
           ll= beta*PostLogLike;
       }
        else{
        ll = beta*PostLogLike+(1-beta)*PropLogLike;
                   }
        
        return new double[] {ll, 0.0};
        
    }
    
    //The GGF part of the density
    protected double[] calcStepLogLikePost(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException {
        int nprime = tA.getNumTaxa()-3;
        double likeStepSD = Math.sqrt(t0/numSteps);
        if (!g.isSimple()) {
            throw new AlgorithmException("Non-simple segment in BridgeWithApproxStepStoneLike");
            //return new double[] {Double.NEGATIVE_INFINITY,0};
        }
        
        double ll = 0.0;
        HashMap<Split, Split> splitMap =  g.getSimpleSplitMap();
        
        Iterator<Split> it = tA.getSplitIterator();
        while (it.hasNext()) {
            Split p = it.next();
            if (p.isTerminal()==null) {
                double la = tA.getSplitLength(p);
                if (!tB.contains(p)) {
                    // This split contracts and expands  
                    double lb = tB.getSplitLength(splitMap.get(p));
                    //GGF likelihood
                    ll += NormalDistribution.logpdf(la+lb, 0.0, likeStepSD)-Math.log(2);
                    
                }
                else {
                    // Shared split
                    
                    double lb = tB.getSplitLength(p);
                    //GGF likelihood
                    ll += NormalDistribution.logpdf(la-lb, 0.0, likeStepSD);
                } 
            }
        }

        return new double[] {ll, 0.0};
        
    }
    

// The independence proposal part of the mixture density
protected double[] calcStepLogLikeProp(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0,int j) throws AlgorithmException {
        int nPrime = tA.getNumTaxa()-3;
        double var=t0/numSteps;
        double RWsd = Math.sqrt(var);
        double sd = Math.sqrt(var*((double)(numSteps-j))/((double)(numSteps-j+1)));
        int numStepsRemaining =numSteps+1-j;
        //double ll=computeFixedMixtureStepDensity(tA,  theTrees[numSteps], tB, sd, numStepsRemaining, unifDist)-nPrime*Math.log(Math.sqrt(2*Math.PI));
        double ll=computeFixedMixtureStepDensity(tA,  theTrees[numSteps], tB, sd, RWsd, numStepsRemaining, unifDist)-(nPrime)*Math.log(Math.sqrt(2*Math.PI));
        
        
        return new double[] {ll, 0.0};
        
        
    }



   
}
    



