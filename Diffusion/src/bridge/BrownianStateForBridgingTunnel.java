/*
BrownianStateForBridgingTunnel
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
import MCMC.RealParameter.RealParameterProposal;
import MCMC.State;
import bridge.BridgeState;
import static bridge.ForwardStepBridge.computeMixtureStepDensity;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.BrownianState;
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
 *Brownian State for bridging that outputs the information needed for both the Chib one block method and the tunnel
 * sampling method of marginal likelihood estimation
 * 
 */
public class BrownianStateForBridgingTunnel extends BrownianStateForBridging{
    
    public final double t_0;
    public final DensityCalculator thePrior;
    public final DensityCalculator theRefDist;
    DoubleUniform unifDist;
    //set the following to be true when simulating from the unnormalised independence proposal and false when simulating from the 
    //Brownian bridge distribution
    public boolean proposal;

    
    public BrownianStateForBridgingTunnel(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, ForwardStepBridge template) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        thePrior = null;
        theRefDist =null;
        proposal =false;
    }
    
    //The below is constructor used when integrating out t0 as well as bridges
    public BrownianStateForBridgingTunnel(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, ForwardStepBridge template,DensityCalculator t0Prior,DensityCalculator t0RefDist,boolean isProposal) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        thePrior = t0Prior;
        theRefDist =t0RefDist;
        proposal =isProposal;
    }
    
    @Override
     public String getValueString() {
         
       String L;
        if(Objects.isNull(thePrior)) L= getValueStringWithoutDisp();
        else    try {
            //then we are integrating out t0 as well as the bridges
            L=getValueStringWithDisp();
       } catch (AlgorithmException ex) {
           Logger.getLogger(BrownianStateForBridgingTunnel.class.getName()).log(Level.SEVERE, null, ex);
           L="";
       }
        
        return L;
        
     }
     
     //for multiple data points and no t0, we output the log likelihoods and prop densities separately for each bridge
     //propDens and then loglike
     public String getValueStringWithoutDisp() {
         
         int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
         
        String L ="";
        //double to store the log independence prop density in
        double l = 0.0;
        double lTotal =0.0;
        //double to store the log likelihood for thde bridge in
        double s = 0.0;
        double sTotal= 0.0;
        for (int i=0; i<numBridges; i++) {
            if(proposal){
            try {
                /*if proposal is true then the likelihood stored in the bridge is the independence proposal likelihood
                in practice this part is not used because we simulate directly from the independence proposal when
                t0 is not included
                */
                l = getBridgeState(i).theBridge.getTotalLogLike();//indep prop dens
                s = getBridgeState(i).theBridge.getTotalGGFLogLike(t_0);//bridge log likelihood
                L+= Double.toString(l)+" "+Double.toString(s)+" ";
                lTotal+=l;
                sTotal+=s;
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("next one");
                   double t= getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(t_0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                    
                }
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingTunnel.class.getName()).log(Level.SEVERE, null, ex);
            }
                
            }
            else{
                try {
                l = getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(t_0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;//indep prop dens
                s = getBridgeState(i).theBridge.getTotalLogLike();//bridge loglike
                lTotal+=l;
                sTotal+=s;
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("next one");
                   double t= getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(t_0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                    
                }
                if(i==numBridges-1){
                L+= Double.toString(l)+" "+Double.toString(s);
                }
                else L+=Double.toString(l)+" "+Double.toString(s)+" ";
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingTunnel.class.getName()).log(Level.SEVERE, null, ex);
            }
            }
        }
        return L; 
     }
     
    public String getValueStringWithDisp() throws AlgorithmException {
        
        int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
        double currt0 = gett0State().getValue();
        
        String L ="";
        //double to store the log independence prop density in
        double l = 0.0;
        double lTotal =0.0;
        //double to store the log likelihood for thde bridge in
        double s = 0.0;
        double sTotal= 0.0;
        
        for (int i=0; i<numBridges; i++) {
            if(proposal){
            try {
                //System.out.println(currt0);
                l = getBridgeState(i).theBridge.getTotalLogLike();//indep prop dens
                s = getBridgeState(i).theBridge.getTotalGGFLogLike(currt0);//bridge log like
                lTotal+=l;
                sTotal+=s;
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("next one");
                   double t= getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(currt0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                    
                }
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingTunnel.class.getName()).log(Level.SEVERE, null, ex);
            }
            }
            else{
                try {
                //System.out.println(currt0);
                l = getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(currt0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;////indep prop dens
                s = getBridgeState(i).theBridge.getTotalLogLike();//bridge log like
                lTotal+=l;
                sTotal+=s;
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("next one");
                   double t= getBridgeState(i).theBridge.computeIndependenceLogDensity(currt0)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                    
                }
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingTunnel.class.getName()).log(Level.SEVERE, null, ex);
            }
            }
        }
        
        //create a new globalState to get the log densities from the prior
        BrownianStateForBridgingTunnel t0State=this;
        Set<State> t0StateList = new HashSet();
        t0StateList.add(t0State);
        GlobalState newState = new GlobalState(t0StateList);
        Double PriorLogDens = thePrior.logDensity(newState, "dispersion");
        Double RefLogDens = theRefDist.logDensity(newState, "dispersion");
        
        L =Double.toString(lTotal) +" "+Double.toString(sTotal)+" "+Double.toString(RefLogDens)+" "+Double.toString(PriorLogDens)+" "+currt0;
        return L; 
     }
     
      @Override
      public String getHeader() {
          String L ="";
      if(Objects.isNull(thePrior)){
          for(int i=0; i<numBridges;i++){
          L+=getBridgeState(i).getName()+"PropDens "+getBridgeState(i).getName()+"LogLike ";
      }
      }
      else L="PropDens LogLike t0RefDens t0PriorDens Dispersion";//don't need to output separate values for the different bridges
      
    return L;
}
      
     
      
 
//Prior that gives access to both the prior distribution and reference distribution on t0
   static public class TunnelBrownianStatePrior implements DensityCalculator {

        /* Instance variables: a prior for the TreeState and a Prior for the standard deviation */
        private DensityCalculator variancePrior;
        private DensityCalculator varianceProposal;
        //set to true if simulating from the independence proposal with simple steps and to false if
        //simulating from the prior.
        private boolean isProposal;

        /* Constructors */

        /** Construct from a prior for TreeState and a prior for PositiveParameter */
        public TunnelBrownianStatePrior(DensityCalculator v,DensityCalculator vp,boolean isProp) {
            isProposal = isProp;
            variancePrior = v;
            varianceProposal = vp;
        }
        
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName,1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            State subst = x.getSubState(subStateName);
            double l;
            if(isProposal==false) l = variancePrior.logDensity(x,subst.getSubState(PositiveParameter.class).getName());//prior density
            else l = varianceProposal.logDensity(x,subst.getSubState(PositiveParameter.class).getName());//reference dist density
            return l;
           
        }

        @Override
        public DensityCalculator makeCopy() {
            return new TunnelBrownianStatePrior(variancePrior,varianceProposal,isProposal);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof BrownianState)) {
                throw new AlgorithmException("SimpleBrownianStatePrior not compatible with given state.");
            }
            variancePrior.checkCompatibility(x,x.getSubState(subStateName).getSubState(PositiveParameter.class).getName());
        }


    }      
   
   /*Prior for Chib two block methods. Density comes from the prior on t0 but also 
   have access to the proposal on t0 given by a RealParameterProposal
   */
   static public class ChibBrownianStatePrior implements DensityCalculator {

        /* Instance variables: a prior for dispersion and a proposal for dispersion */
        private DensityCalculator variancePrior;
        private RealParameterProposal varianceProposal;


        /* Constructor */
         public ChibBrownianStatePrior(DensityCalculator v,RealParameterProposal vp) {
            variancePrior = v;//the prior on t0
            varianceProposal = vp;//the proposal on t0
        }
        
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName,1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            State subst = x.getSubState(subStateName);
            double l = variancePrior.logDensity(x,subst.getSubState(PositiveParameter.class).getName());
            return l;
           
        }

        @Override
        public DensityCalculator makeCopy() {
            return new ChibBrownianStatePrior(variancePrior,varianceProposal);
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
