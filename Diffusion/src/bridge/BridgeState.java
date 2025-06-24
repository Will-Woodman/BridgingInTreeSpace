/*
   BridgeState
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

import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import MCMC.PositiveParameter;
import MCMC.State;
import cern.jet.random.tdouble.DoubleUniform;
import java.util.ArrayList;
import simulation.CategoricalDistribution;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;

/**
    MCMC state representing a bridge between source x0 and data point
 */

public class BridgeState extends State {

    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /** Instance variables */
    public ForwardStepBridge theBridge;
    
    final private static boolean DEBUG = true;

    /* CONSTRUCTORS --------------------------------------------------------- */

    /** Initialize with pre-existing bridge */
    public BridgeState(ForwardStepBridge b, String stateName) {
        super(stateName);
        theBridge = b;
    }


    /* UTILITY METHODS ------------------------------------------------------ */

    /** Return list of strings -- one for each step on RW.
     Useful for debug! */
    public String[] getTreeList() {
        ArrayList<TreeAsSplits> a = theBridge.getTreePathList();
        String[] res = new String[a.size()];
        for (int i=0; i<a.size(); i++) {
            res[i] = a.get(i).toString();
        }
        return res;
    }

    /** Output the loglike for this path
     * @return  */
    //changing this now to ouput the trees from the bridge
    /*
    public String getValueString() {
        return String.format("%7.7f", theBridge.getTotalLogLike());
    }
    */
     @Override
     public String getValueString() {
        ArrayList<TreeAsSplits> a = theBridge.getTreePathList();
        String res = "";
        for (int i=0; i<a.size(); i++) {
            res = res + a.get(i).toString();
        }
        return res;
     }


    @Override
    public void setValueFromString(String str) throws AlgorithmException {
        throw new UnsupportedOperationException("Not supported.");
    }
    @Override
    public int getLengthStringRepr() {
        return 1; // Single block of text
    }

    @Override
    public void updateRunningMean(State runningMean, int runningCount) throws AlgorithmException {
        throw new UnsupportedOperationException("Not supported.");
    }
    
  
    
    /* LIKELIHOOD CALCULATOR ------------------------------------------------ */
        
    static public class BridgeLikelihoodCalculator implements DensityCalculator {
        
        public BridgeLikelihoodCalculator() {
            // Nothing to do
        }
        
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BridgeState)) {
                throw new AlgorithmException("Likelihood not compatible with given state.");
            }
//            if (!subStateName.equals(BridgeState.this.getName())) {
//                throw new AlgorithmException("Bridge likelihood calculator not compatible with given state.");
//            }
        }
       
        @Override
        /* When would this ever be needed??? */
        public DensityCalculator makeCopy() {
            //DEBUG
            System.out.println("Suprise! Call to makeCopy in BrownianBridgeState likelihood.");
            return new BridgeLikelihoodCalculator();
        }
        
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            double l = logDensity(x,subStateName);
            return l*likelihoodTemperature;
        }

        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BridgeState)) {
                throw new AlgorithmError("Bridge likelihood calculator not compatible with given state.");
            }

            return ((BridgeState)x.getSubState(subStateName)).theBridge.getTotalLogLike();
         }


    }


    /* PROPOSALS  ----------------------------------------------------------- */

    /** Independence proposal */
    static protected class IndependenceProposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
                
        private double likelihoodIncrement; 
        private int numCalls = 0, numFails = 0;


        public IndependenceProposal() {
            name = "Brownian bridge independence proposal";
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BridgeState)) {
                throw new AlgorithmException("Bridge independence proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
        }

        @Override
        public Kernel makeCopy() {
            return new IndependenceProposal();
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            numCalls++;

            BridgeState current = (BridgeState) theState.getSubState(subStateName);
            ForwardStepBridge b = current.theBridge;
            double t0 = ((PositiveParameter)theState.getSubState("dispersion")).getValue();
            double clp;
            try {
                clp = b.computeIndependenceLogDensity(t0);
            }
            catch (AlgorithmException failed) {
                // Reject -- current state \emph{could not} have arisen from an indept proposal!
                if (DEBUG) System.out.println("Reject independence proposal due to current state. ");
                likelihoodIncrement = 0.0;
                numFails++;
                return Double.NEGATIVE_INFINITY;
            }
            double plp=0.0;
            double cll = b.getTotalLogLike();            
            try {
                plp = b.independenceProposal(t0);
                likelihoodIncrement = b.getTotalLogLike()-cll;
            }
            catch (AlgorithmException failed) {
                /* Now ensure proposal is rejected */
                numFails++;
                likelihoodIncrement = 0.0;
                return Double.NEGATIVE_INFINITY; 
            }

            
//System.out.println(clp+" "+plp+" "+b.getTotalLogLike()+" "+cll);
            
            return clp - plp;

        }
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(theBackup.getLogLikelihood()+likelihoodIncrement);
            return ratio;
        }


        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BridgeState curr = (BridgeState) theState.getSubState(subStateName);
            BridgeState back = (BridgeState) theBackup.getSubState(subStateName);
            if(whichWay) {
                curr.theBridge.copyTo(back.theBridge);
            }
            else {
                back.theBridge.copyTo(curr.theBridge);
            }
        }
        
        public double getPropFails() {
            return ((double) numFails)/((double) numCalls);
        }

    }
    
    
    /** Partial proposal */
    static protected class PartialBridgeProposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        /* Random number generators */
        private DoubleUniform unifDist;
        private CategoricalDistribution catDist;
        
        private int a, b;
        private double likelihoodIncrement;
        private int fixedLen=0;
        private int fixedStart=-1, fixedEnd=-1;
        
        private int numCalls = 0, numFails = 0;

        public PartialBridgeProposal() {
            name = "Partial bridge proposal";
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            catDist = null;
        }

        public PartialBridgeProposal(int l) {
            name = "Partial bridge proposal";
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            catDist = null;
            fixedLen = l;
            fixedStart = -1;
        }

        public PartialBridgeProposal(int a, int b) {
            name = "Partial bridge proposal";
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            catDist = null;
            fixedLen = -1;
            fixedStart = a;
            fixedEnd = b;
        }

        public PartialBridgeProposal(CategoricalDistribution cd) {
            name = "Partial bridge proposal";
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            catDist = cd;
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BridgeState)) throw new AlgorithmException("Bridge hierarchy proposal not compatible with given state.");
        }

        @Override
        public void resetRandomEngineSeed() {
            unifDist = new DoubleUniform(simulation.Random.getEngine());
        }
        @Override
        public Kernel makeCopy() {
            return new PartialBridgeProposal();
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            numCalls++;

            BridgeState current = (BridgeState) theState.getSubState(subStateName);
            double t0 = ((PositiveParameter)theState.getSubState("dispersion")).getValue();
           
            /* First sample where to make the replacement */
            int numSteps = current.theBridge.getNumSteps();
            int len=1;
            
            if (fixedStart>0) {
                a = fixedStart; 
                b = fixedEnd;
            } 
            else {
                if (catDist==null) {
                    if (fixedLen==0)    
                        len = unifDist.nextIntFromTo(1, numSteps-1);
                    else 
                        len = fixedLen;
                }
                else {
                    len = 1+catDist.sample();
                }

                if (len==numSteps-1) {
                    a = 1; b = numSteps;
                }
                else {
                    a = unifDist.nextIntFromTo(1, numSteps-len);
                    b = a+len;
                }
            }

            double clp;
            try {
                clp = current.theBridge.computePartialBridgeLogDensity(t0, a, b);
            }
            catch (AlgorithmException failed) {
                // Reject -- current state \emph{could not} have arisen from this partial proposal!
                //replacing the below so I can see an example of the partial prop being rejected
                //if (DEBUG) System.out.println("Reject partial proposal due to current state. ");
                if (DEBUG){
                    System.out.println("Reject partial proposal due to current state. ");
                    System.out.println(current.theBridge.getTreePathList().toString());
                    System.out.println(a);
                    System.out.print(" ");
                    System.out.print(b);
                }
                numFails++;
                likelihoodIncrement = 0.0;
                return Double.NEGATIVE_INFINITY;
            }
                        
            double cll = current.theBridge.getTotalLogLike();
            double plp=0.0;
                                    
            try {
                plp = current.theBridge.partialBridgeProposal(t0, a, b);
                likelihoodIncrement = current.theBridge.getTotalLogLike()-cll;
            }
            catch (AlgorithmException failed) {
                /* Now ensure proposal is rejected */
                likelihoodIncrement = 0.0;
                numFails++;
                return Double.NEGATIVE_INFINITY; 
            }
            
            double ratio =  clp-plp;

            return ratio;

        }
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(theBackup.getLogLikelihood()+likelihoodIncrement);
            return ratio;
        }


        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BridgeState curr = (BridgeState) theState.getSubState(subStateName);
            BridgeState back = (BridgeState) theBackup.getSubState(subStateName);
            if(whichWay) {
                curr.theBridge.copyTo(back.theBridge);
            }
            else {
                back.theBridge.copyTo(curr.theBridge);
            }
        }
        
        public double getPropFails() {
            return ((double) numFails)/((double) numCalls);
        }
    }
    
    
}
