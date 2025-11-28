/*
BrownianStateForBridgingChibSecondBlockMultStar
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
import MCMC.OutputStringProcessing;
import MCMC.PositiveParameter;
import MCMC.RealParameter.RealParameterProposal;
import MCMC.State;
import cern.jet.random.tdouble.DoubleUniform;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;

/**
 *
 * @author will - 
 * Used for simulations in the second part of Chib two block marginal likelihood estimations 
 * Output log likelihoods and independence proposal densities for each individual bridge for fixed t_0^* 
 *  * Then generate a new value of t_0 under the proposal at t_0^* and calculate the densities again
 * This class gives the option to generate multiple new values of t_0 under the proposal at t_0^*
 * Could also be used in the standard Chib two block estimator as it does not rely on multiple values of t_0^*
 */

public class BrownianStateForBridgingChibSecondBlockMultStar extends BrownianStateForBridging{
    
    public final double t_0;
    public final DensityCalculator thePrior;
    public final RealParameterProposal theRefDist;

    public double t_0Star;
    //to stabilise the estimate for the third distribution, simulate multiple values of t_0 from rho(t_0|t_0^*):
    int numt_0Props;
    
    DoubleUniform unifDist;
    
    //The below is constructor used when integrating out t0 as well as bridges
    public BrownianStateForBridgingChibSecondBlockMultStar(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, ForwardStepBridge template,DensityCalculator t0Prior,RealParameterProposal t0RefDist,double t0Star,int numT0Props) throws AlgorithmException {
        super(x0, t0, dataPoints, constructBridge, template);
        t_0 =t0;
        thePrior = t0Prior;
        theRefDist =t0RefDist;
        t_0Star=t0Star;
        numt_0Props=numT0Props;
    }
    
    @Override
      public String getValueString(){    
      int nPrime = getBridgeState(0).theBridge.treeA.getNumTaxa()-3;
      
      
        String L ="";
        //double to store the log independence prop density in
        double l = 0.0;
        //double to store the log likelihood for thde bridge in
        double s = 0.0;
        for (int i=0; i<numBridges; i++) {
            try {
                //Compute the log independence proposal density for all of the bridges
                l= getBridgeState(i).theBridge.computeMargLikeIndependenceLogDensity(t_0Star)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                //Compute the log likelihood for all of the bridges
                s= getBridgeState(i).theBridge.getTotalLogLike();
                L+= Double.toString(l)+" "+Double.toString(s)+" ";
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("zero likelihood calculated for bridge in Chib second block");
                   //double t= getBridgeState(i).theBridge.computeIndependenceLogDensityCheck(t_0Star)-Math.log(2*Math.PI)*(numSteps-1)*nPrime/2;
                    
                }
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        
        
        //create a new global state to get the log densities from the prior
        BrownianStateForBridgingChibSecondBlockMultStar t0State=this;
        Set<State> t0StateList = new HashSet();
        t0StateList.add(t0State);
        GlobalState newState = new GlobalState(t0StateList);
        Double PriorT0StarDens=0.0;
        Double[] PriorT0Dens=new Double[numt_0Props];
        Double[] tStarPropDens=new Double[numt_0Props];
        PositiveParameter curr2 = (PositiveParameter) t0State.getSubState("dispersion");
        double currv2 = curr2.getValue();

         //get the prior density for t0
        try {
            PriorT0StarDens = thePrior.logDensity(newState, "dispersion");
        } catch (AlgorithmException ex) {
            Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("Cannot calculate t0 prior density");
        }
        
        //update t0 in the new state using the proposal distribution the specified number of times
        Double[] RefLogDensRatio=new Double[numt_0Props];
        for(int j=0;j<numt_0Props;j++){
        try {
            RefLogDensRatio[j]=theRefDist.sample(newState, "dispersion");
            PriorT0Dens[j]=thePrior.logDensity(newState, "dispersion");
        } catch (AlgorithmException ex) {
            Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("Cannot calculate t0 proposal density");
        }   
         //get the proposed value of t0 and calculate the proposal densities
        PositiveParameter curr = (PositiveParameter) t0State.getSubState("dispersion");
        double currv = curr.getValue();
       
        try {
            tStarPropDens[j] = theRefDist.logProposalDensity(currv, t_0Star);
        } catch (AlgorithmException ex) {
            Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
        }
      
        //finally get the proposal densities and likelihoods for all of the bridges with the new t0 
        curr = (PositiveParameter) newState.getSubState("dispersion");
        currv = curr.getValue();
        
        double sTotal =0.0;
        s=0.0;
        for (int i=0; i<numBridges; i++) {
            try {

                s = getBridgeState(i).theBridge.getTotalGGFLogLike(currv);

                sTotal+=s;
                if(l==Double.NEGATIVE_INFINITY) {
                    System.out.println("zero likelihood calculated for bridge in Chib second block");
                    
                }
            } catch (AlgorithmException ex) {
                Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        
        L+=sTotal+" "+RefLogDensRatio[j]+" "+PriorT0Dens[j]+" ";
        
        //ensure the t_0 state is set back to t_0Star
        try {
            newState.getSubState("dispersion").setValueFromString(Double.toString(t_0Star));
        } catch (AlgorithmException ex) {
            Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("Error setting t0Star back");
        } catch (OutputStringProcessing.InsufficientInformationException ex) {
            Logger.getLogger(BrownianStateForBridgingChibSecondBlockMultStar.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("Error setting t0Star back");
        }
        
        }
        
        return L+PriorT0StarDens; 
        
     }

     
      @Override
      public String getHeader() {
      String L ="";
      for(int i=0; i<numBridges;i++){
          L+=getBridgeState(i).getName()+"PropDens "+getBridgeState(i).getName()+"LogLike ";
      }
      for(int j =0;j<numt_0Props;j++){
          L+="TotalLogLikeNewDisp"+j+" DispPropRatio"+j+" T0Prior"+j+" ";
      }
     
      
    return L+=" T0StarPrior";
}
      
      
}
