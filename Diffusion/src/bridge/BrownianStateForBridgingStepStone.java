/*
BrownianStateForBridgingStepStone
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
import MCMC.Kernel;
import MCMC.PositiveParameter;
import MCMC.State;
import bridge.BridgeState;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.BrownianState;
import diffbase.TreeState;
import geodesics.Geodesic;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.CategoricalDistribution;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;

/**
 *
 * @author will
 */
public class BrownianStateForBridgingStepStone extends BrownianStateForBridging{
    
    public final double t_0;
    public final double betaMinus;
    public final double beta;
    //below will be set to null if not integrating out t0:
    public final SteppingStoneBrownianStatePrior thePrior;
    
    //constructor without existing bridges
    public BrownianStateForBridgingStepStone(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, BridgeWithStepStoneLike template,double Beta,double BetaMinus) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        beta=Beta;
        betaMinus = BetaMinus;
        thePrior=null;
    }
    
    //constructor with already existing bridges - used for the quasistatic stepping stone sampling
    public BrownianStateForBridgingStepStone(TreeAsSplits x0, double t0, int m, ArrayList<BridgeState> bridges,double Beta,double BetaMinus) {
        super(x0,t0,m,bridges);
        t_0=t0;
        beta=Beta;
        betaMinus = BetaMinus;
        thePrior=null;
        
    }
    
    //constructor when integrating out t0 as well as bridges
    public BrownianStateForBridgingStepStone(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, BridgeWithStepStoneLike template,double Beta,double BetaMinus, SteppingStoneBrownianStatePrior t0Prior) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        beta=Beta;
        betaMinus = BetaMinus;
        thePrior=t0Prior;
    }
    
    /*constructor with already existing bridges - used for the quasistatic stepping stone sampling -
    constructor when integrating out t0 as well as bridges
    */
    public BrownianStateForBridgingStepStone(TreeAsSplits x0, double t0, int m, ArrayList<BridgeState> bridges,double Beta,double BetaMinus, SteppingStoneBrownianStatePrior t0Prior) {
        super(x0,t0,m,bridges);
        t_0=t0;
        beta=Beta;
        betaMinus = BetaMinus;
        thePrior=t0Prior;
        
    }
    
    @Override
    
    public String getValueString() {
        String L;
        //if there is no distribution on t0, we are not integrating out t0:
        if(Objects.isNull(thePrior)) L= getValueStringWithoutDisp();
        else{ try {
            L=getValueStringWithDisp();
        } catch (AlgorithmException ex) {
            Logger.getLogger(BrownianStateForBridgingStepStone.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("Something wrong producing t0 dist results");
            L= "";
        }
        
     }
    return L+" "+gett0State().getValueString();
}
    
     public String getValueStringWithoutDisp() {
        String L ="";
        Double l = 0.0;
        double lTotal=0.0;
        int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
        for (int i=0; i<numBridges; i++) {
            try {
               //get the stepping stone likelihood 
               double StepLike=getBridgeState(i).theBridge.getTotalLogLike();
               //get the reference distribution likelihood
               double PropDen =getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(t_0)-Math.log(Math.sqrt(2*Math.PI))*(numSteps-1)*nPrime;
               double GGFLike;
               //get the usual bridge likelihood:
               if(betaMinus==0){
                   GGFLike=getBridgeState(i).theBridge.getTotalGGFLogLike(t_0);
               }
               else{
                   //otherwise we can work it out from the Stepping stone likelihood and proposal density:
               GGFLike= (StepLike-(1-betaMinus)*PropDen)/betaMinus;
               }
               //calculate what we need for the estimate:
                l = (beta-betaMinus)*(GGFLike-PropDen);
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingStepStone.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            //lTotal +=l;
            L+=l.toString()+" ";
        }
        return L+gett0State().getValueString();
     }
     
     public String getValueStringWithDisp() throws AlgorithmException {
        String L ="";
        double l = 0.0;
        double lTotal=0.0;
        int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
        for (int i=0; i<numBridges; i++) {
            try {
               double currt0 = gett0State().getValue();
               //get the stepping stone likelihood 
               double StepLike=getBridgeState(i).theBridge.getTotalLogLike();
               //get the reference distribution likelihood
               double PropDen =getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(currt0)-Math.log(Math.sqrt(2*Math.PI))*(numSteps-1)*nPrime;
               double GGFLike;
               //get the usual bridge likelihood:
               if(betaMinus==0){
                   GGFLike=getBridgeState(i).theBridge.getTotalGGFLogLike(currt0);
               }
               else{
                   //otherwise we can work it out from the Stepping stone likelihood and proposal density:
               GGFLike= (StepLike-(1-betaMinus)*PropDen)/betaMinus;
               }
                l= (beta-betaMinus)*(GGFLike-PropDen);
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingStepStone.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            lTotal +=l;
        }
        //create a new global state to get the log densities from the prior and the reference distribution on t0
        BrownianStateForBridgingStepStone t0State=this;
        Set<State> t0StateList = new HashSet();
        t0StateList.add(t0State);
        GlobalState newState = new GlobalState(t0StateList);
        //get the densities from both the reference distribution and prior on t0
        Double PriorLogDens = thePrior.variancePrior.logDensity(newState, "dispersion");
        Double RefLogDens = thePrior.varianceProposal.logDensity(newState, "dispersion");
        
        lTotal+=(beta-betaMinus)*(PriorLogDens-RefLogDens);
        L=Double.toString(lTotal);
        return L+" "+gett0State().getValueString()+" "+t_0;
     }
     
     @Override
      public String getHeader() { 
          String L ="";
      if(Objects.isNull(thePrior)){
          for(int i=0; i<numBridges;i++){
          L+=getBridgeState(i).getName()+"Ratio ";
      }
      }
      else {
          L= "bridgeRatio"+" "+gett0State().getName();
      }
      return L;
}
   
    //a prior specifically for stepping stone that allows just a prior on t0 and implements the prior multiplied by
    //combined with the t0 marginal reference distribution
   static public class SteppingStoneBrownianStatePrior implements DensityCalculator {

        /* Instance variables: a prior for the TreeState and a Prior for the standard deviation */
        private DensityCalculator variancePrior;
        private DensityCalculator varianceProposal;
        private double temp;

        /* Constructors */

        /** Constructor
         temp is the value beta (see paper)
         */
        public SteppingStoneBrownianStatePrior(DensityCalculator v,DensityCalculator vp,double Temp) {
            variancePrior = v;
            varianceProposal = vp;
            temp=Temp;
            System.out.println(Temp+" temp");
        }

        
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName,temp, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            State subst = x.getSubState(subStateName);
            // log density from the prior
            double l = variancePrior.logDensity(x,subst.getSubState(PositiveParameter.class).getName());
            //log density from the reference distribution
            double lp = varianceProposal.logDensity(x,subst.getSubState(PositiveParameter.class).getName());
            return lp*(1-temp)+l*temp;
        }

        @Override
        public DensityCalculator makeCopy() {
            return new SteppingStoneBrownianStatePrior(variancePrior,varianceProposal,temp);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof BrownianState)) {
                throw new AlgorithmException("SimpleBrownianStatePrior not compatible with given state.");
            }
            variancePrior.checkCompatibility(x,x.getSubState(subStateName).getSubState(PositiveParameter.class).getName());
        }


    }
     
     

  
}
