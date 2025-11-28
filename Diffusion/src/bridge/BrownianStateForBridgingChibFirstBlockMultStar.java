/*
BrownianStateForBridgingChibFirstBlockMultStar
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
 *Brownian State for bridging that outputs the information needed for both the Chib two block method with multiple fixed values
 * of t_0^* - not included in thesis
 * 
 */
public class BrownianStateForBridgingChibFirstBlockMultStar extends BrownianStateForBridging{
    
    public final double t_0;
    public final DensityCalculator thePrior;
    public final RealParameterProposal theRefDist;
    public Double[] t_0Star;
    DoubleUniform unifDist;
  
    //The below is constructor used when integrating out t0 as well as bridges
    public BrownianStateForBridgingChibFirstBlockMultStar(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, ForwardStepBridge template,DensityCalculator t0Prior,RealParameterProposal t0RefDist,Double[] t0Star) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        thePrior = t0Prior;
        theRefDist =t0RefDist;
        t_0Star=t0Star;
    }
    
    @Override
     public String getValueString() {
         
       String L;
        if(Objects.isNull(thePrior)){
            L="";
            System.out.println("Attempt to use the Chib two block estimator without a prior on t0. Exiting");
                System.exit(1);
        }
        else    try {
            L=getValueStringWithDisp();
       } catch (AlgorithmException ex) {
           Logger.getLogger(BrownianStateForBridgingChibFirstBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
           L="";
       }
        
        return L;
        
     }
     
    public String getValueStringWithDisp() throws AlgorithmException {
         
        int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
         
        String L ="";
        int numT0Stars=t_0Star.length;
        //double to store the log independence prop density in
        double[] l = new double[numT0Stars];
        double[] lTotal =new double[numT0Stars];
        //double to store the log likelihood for thde bridge in
        double s = 0.0;
        double sTotal= 0.0;
        for (int i=0; i<numBridges; i++) {
            try {
                //l = getBridgeState(i).theBridge.computeIndependenceLogDensity(t_0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                s = getBridgeState(i).theBridge.getTotalLogLike();
                //compute the log likelihoods for each of the different fixed values of t_0^*
                for(int j=0;j<numT0Stars;j++){
                l[j] = getBridgeState(i).theBridge.getTotalGGFLogLike(t_0Star[j]);
                lTotal[j]+=l[j];
                if(l[j]==Double.NEGATIVE_INFINITY) {
                    System.out.println("zero likelihood calculated for bridge in Chib first block");
                   //double t= getBridgeState(i).theBridge.computeIndependenceLogDensityCheck(t_0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                    
                }
                }
                sTotal+=s;
                
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingChibFirstBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        //create a new state with just t0 to get the log densities from the prior
        BrownianStateForBridgingChibFirstBlockMultStar t0State=this;
        Set<State> t0StateList = new HashSet();
        t0StateList.add(t0State);
        GlobalState newState = new GlobalState(t0StateList);
        PositiveParameter curr = (PositiveParameter) t0State.getSubState("dispersion");
        double currv = curr.getValue();
        Double[] RefLogDens= new Double[numT0Stars];
        Double[] tStarPropDens= new Double[numT0Stars];
        Double PriorLogDens = thePrior.logDensity(newState, "dispersion");

                
        //calculate the proposal densities on t0
        for(int j=0;j<numT0Stars;j++){
        RefLogDens[j] = theRefDist.logProposalDensity(t_0Star[j],currv);
        tStarPropDens[j] = theRefDist.logProposalDensity(currv,t_0Star[j]);
            L+=Double.toString(lTotal[j])+" "+Double.toString(RefLogDens[j])+" "+Double.toString(tStarPropDens[j])+" ";
        }
        L +=Double.toString(sTotal)+" "+Double.toString(PriorLogDens);
        return L; 
     }
     
      @Override
      public String getHeader() {
          String L ="";
        if(Objects.isNull(thePrior)){
            L="";
            System.out.println("Attempt to use the Chib two block estimator without a prior on t0. Exiting");
                System.exit(1);
        }
      else{
          int numT0Stars=t_0Star.length;
          for(int j=0;j<numT0Stars;j++){
              L+= "LogLikeStar"+j+" t0PropDens"+j+" t0StarPropDens"+j+" ";
          }
          L+="LogLike t0PriorDens";
      }
      
    return L+" "+gett0State().getName();
      }
}
      
