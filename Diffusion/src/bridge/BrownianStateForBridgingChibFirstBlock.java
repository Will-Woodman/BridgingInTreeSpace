/*
BrownianStateForBridgingChibFirstBlock
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

import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.PositiveParameter;
import MCMC.RealParameter.RealParameterProposal;
import MCMC.State;
import cern.jet.random.tdouble.DoubleUniform;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;

/**

 *Brownian State for bridging that outputs the information needed for Chib one two block estimator from
 * the full posterior MCMC sims
 * 
 */

public class BrownianStateForBridgingChibFirstBlock extends BrownianStateForBridging{
    
    public final double t_0;
    public final DensityCalculator thePrior;
    public final RealParameterProposal theRefDist;
    public double t_0Star;

    
    DoubleUniform unifDist;
    
    //Constructor
    public BrownianStateForBridgingChibFirstBlock(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, ForwardStepBridge template,DensityCalculator t0Prior,RealParameterProposal t0RefDist,double t0Star) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        thePrior = t0Prior;
        theRefDist =t0RefDist;
        t_0Star=t0Star;
    }
    
    @Override
     public String getValueString() {
         
       String L;
       //Check that the two block estimator is being used in the correct context:
        if(Objects.isNull(thePrior)){
            L="";
            System.out.println("Attempt to use the Chib two block estimator without a prior on t0. Exiting");
                System.exit(1);
        }
        else    try {
            L=getValueStringWithDisp();
       } catch (AlgorithmException ex) {
           Logger.getLogger(BrownianStateForBridgingChibFirstBlock.class.getName()).log(Level.SEVERE, null, ex);
           L="";
       }
        
        return L;
        
     }
     
    /*
    Class to output all the info needed to calculate one part of the Chib two block estimate
    Used for the simulations for the full posterior (on t0 and bridges)
    Involves computing likelihoods for current t0 and fixed value t_0^*
    */
    public String getValueStringWithDisp() throws AlgorithmException {
         
        int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
         
        String L ="";
        //double to store the log independence prop density in
        double l = 0.0;
        double lTotal =0.0;
        //double to store the log likelihood for the bridge in
        double s = 0.0;
        double sTotal= 0.0;
        for (int i=0; i<numBridges; i++) {
            try {
                s = getBridgeState(i).theBridge.getTotalLogLike();
                l = getBridgeState(i).theBridge.getTotalGGFLogLike(t_0Star);
                lTotal+=l;
                sTotal+=s;
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("zero likelihood calculated for bridge in Chib first block");
                   //double t= getBridgeState(i).theBridge.computeIndependenceLogDensityCheck(t_0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2; 
                }
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingChibFirstBlock.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        //create a new state with just t0 to get the log densities from the prior
        BrownianStateForBridgingChibFirstBlock t0State=this;
        Set<State> t0StateList = new HashSet();
        t0StateList.add(t0State);
        GlobalState newState = new GlobalState(t0StateList);
        PositiveParameter curr = (PositiveParameter) t0State.getSubState("dispersion");
        double currv = curr.getValue();
        Double PriorLogDens = thePrior.logDensity(newState, "dispersion");
        //Now calculate densities of the proposals on t0
        //Proposal density of t0 given state t_0^*:
        Double RefLogDens = theRefDist.logProposalDensity(t_0Star,currv);
        //Proposal density of t0^* given state t_0:
        Double tStarPropDens = theRefDist.logProposalDensity(currv,t_0Star);
        
        //Form the output string
        L =Double.toString(lTotal) +" "+Double.toString(sTotal)+" "+Double.toString(RefLogDens)+" "+Double.toString(PriorLogDens);
        L+= " "+Double.toString(tStarPropDens);
        return L; 
     }
     
      @Override
      public String getHeader() {
          String L ="";
      if(Objects.isNull(thePrior)){
            System.out.println("Attempt to use the Chib two block estimator without a prior on t0. Exiting");
                System.exit(1);
        }
      else L="LogLikeStar LogLike t0RefDens t0PriorDens t0StarPropDens";
      
    return L+" "+gett0State().getName();
}
     
      
      public String getHeaderWithDisp() {
          String L ="";
      
        for (int i=0; i<numBridges; i++) {
    if(i==0) L+= "PropDens"+i+" Like"+i;
    else L+= " PropDens"+i+" Like"+i;
        
      } 
    return L+" "+gett0State().getName();
}
}
