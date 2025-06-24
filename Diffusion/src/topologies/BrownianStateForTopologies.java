/*
   BrownianStateForTopologies
    Copyright (C) 2018  Tom M. W. Nye

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

package topologies;

/**
    MCMC state representing the initial source tree and dispersion parameter for 
    a Brownian motion in tree-space.
    
    State parameters are 
    * x_0 -- tree;
    * t_0 -- dispersion
  
    Used for fitting to data sets of topologies
    
    Note: t_0 should be kept fixed due to identifiability. 
    Only x_0 changes. Proposals etc which use t_0 are included only for backwards compatibility
    
 */

import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import MCMC.PositiveParameter;
import MCMC.State;
import diffbase.BrownianState;
import diffbase.TreeState;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import simulation.ExponentialDistribution;
import simulation.GammaDistribution;
import simulation.NormalDistribution;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.TreeAsSplits;
import treedatasets.OperationsOnTreeAsSplits;


public class BrownianStateForTopologies extends BrownianState {

    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    final private static boolean DEBUG = true;

    
    /* CONSTRUCTORS --------------------------------------------------------- */

    public BrownianStateForTopologies(TreeAsSplits x0, double t0)  {
        super(x0, t0);
    }


    /* UTILITY METHODS ------------------------------------------------------ */
    
     
    /* PROPOSALS ------------------------------------------------------------ */
    
    /** Wrapper for any proposal on t0 alone.
     Default is log normal RW. 
     Do not use! */
    public class T0Proposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        private Kernel t0Proposal;
       
        /** Constructor  */
        public T0Proposal(Kernel t0Prop) {
            name = "Brownian state t0 proposal retaining bridges";
            t0Proposal = t0Prop;
        }

        public T0Proposal(double sd) {
            name = "Brownian state t0 proposal retaining bridges";
            t0Proposal = new PositiveParameter.LogNormalProposal(sd);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            t0Proposal.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new T0Proposal(t0Proposal);
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            BrownianStateForTopologies currentBrownianState = (BrownianStateForTopologies) theState.getSubState(subStateName);
            
            double t0PropRatio = t0Proposal.sample(theState, "dispersion");
                                     
            return t0PropRatio;    
        }
            
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            return ratio;
        }


        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            BrownianStateForTopologies back = (BrownianStateForTopologies) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store t0  */
                back.sett0(curr.gett0());
            }
            else {
                /* Restore t0  */
                curr.sett0(back.gett0());
            }
        }
        
    }
    
    
    /** Proposal on x0.
        Acts as a wrapper for a TreeState proposal. */
    public static class X0Proposal extends Kernel {

        private Kernel x0Proposal;

        /** Constructor  */
        public X0Proposal(Kernel x0prop) {
            x0Proposal = x0prop;
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            x0Proposal.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new X0Proposal(x0Proposal);
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            
            double x0PropRatio = x0Proposal.sample(theState, "source");
            return x0PropRatio;       
        }
            

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            BrownianStateForTopologies back = (BrownianStateForTopologies) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store x0 */
                back.setx0(curr.getx0().efficientClone());
            }
            else {
                /* Restore x0 */
                curr.setx0(back.getx0().efficientClone());
            }
        }
     

    }


    /** Proposal on x0.
        Get next edge. Use log-RW on edge length. */
    public static class LogNormalEdgeProposal extends Kernel {

        private NormalDistribution normDist;
        private int currentSplitInd;
        private double oldLength;
        private Split theSplit;

        /** Constructor  */
        public LogNormalEdgeProposal(double sd) {
            normDist = new NormalDistribution(0, sd);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            normDist.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new LogNormalEdgeProposal(normDist.getSigma());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            TreeAsSplits t = curr.getx0();

            ArrayList<Split> orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");
            
            theSplit = orderedSplits.get(currentSplitInd);
            oldLength = t.getSplitLength(theSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above

            // Now sample a new length for this edge
            double newLength = oldLength*Math.exp(normDist.sample());
            t.setSplitLength(theSplit, newLength);
            
            // Update splitInd
            currentSplitInd++;
            if (currentSplitInd>=orderedSplits.size()) {
                currentSplitInd=0;
            }

            return Math.log(newLength)-Math.log(oldLength); // Log proposal ratio
        }
            

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            if(whichWay) {
                /* Nothing to do */
            }
            else {
                /* Restore x0 */
                try {
                    curr.getx0().setSplitLength(theSplit, oldLength);
                }
                catch (AlgorithmError anErr) {
                    System.out.println("Algorithm error restoring split length after reject in BrownianStateForTopologies.LogNormalEdgeProposal. This should not be possible.");
                }
            }
        }
     

    }

     public static class GammaEdgeProposal extends Kernel {

        private GammaDistribution gammaDist;
        private int currentSplitInd;
        private double oldLength;
        private Split theSplit;

        /** Constructor  */
        public GammaEdgeProposal(double shape, double scale) {
            gammaDist = new GammaDistribution(shape, scale);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            gammaDist.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new GammaEdgeProposal(gammaDist.getShape(),gammaDist.getScale());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            TreeAsSplits t = curr.getx0();

            ArrayList<Split> orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");
            
            theSplit = orderedSplits.get(currentSplitInd);
            oldLength = t.getSplitLength(theSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above

            // Now sample a new length for this edge
            double newLength = oldLength*gammaDist.sample();
            t.setSplitLength(theSplit, newLength);
            
            // Update splitInd
            currentSplitInd++;
            if (currentSplitInd>=orderedSplits.size()) {
                currentSplitInd=0;
            }

            return (2.0*gammaDist.getShape()-2.0)*(Math.log(oldLength) - Math.log(newLength)) - 1/gammaDist.getScale()*(oldLength/newLength - newLength/oldLength); // Log proposal ratio
        }
            

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            if(whichWay) {
                /* Nothing to do */
            }
            else {
                /* Restore x0 */
                try {
                    curr.getx0().setSplitLength(theSplit, oldLength);
                }
                catch (AlgorithmError anErr) {
                    System.out.println("Algorithm error restoring split length after reject in BrownianStateForTopologies.LogNormalEdgeProposal. This should not be possible.");
                }
            }
        }
     

    }

    
        /** Proposal on x0.
        Get next edge. Use log-RW on edge length. */
    public static class ExponentialIndepProposal extends Kernel {

        private ExponentialDistribution expDist;
        private int currentSplitInd;
        private double oldLength;
        private Split theSplit;

        /** Constructor  */
        public ExponentialIndepProposal(double rate) {
            expDist = new ExponentialDistribution( rate);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            expDist.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new ExponentialIndepProposal(expDist.getRate());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            TreeAsSplits t = curr.getx0();

            ArrayList<Split> orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");
            
            theSplit = orderedSplits.get(currentSplitInd);
            oldLength = t.getSplitLength(theSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above

            // Now sample a new length for this edge
            double newLength = expDist.sample();
            t.setSplitLength(theSplit, newLength);
            
            // Update splitInd
            currentSplitInd++;
            if (currentSplitInd>=orderedSplits.size()) {
                currentSplitInd=0;
            }

            return -expDist.getRate()*(oldLength-newLength); // Log proposal ratio
        }
            

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            if(whichWay) {
                /* Nothing to do */
            }
            else {
                /* Restore x0 */
                try {
                    curr.getx0().setSplitLength(theSplit, oldLength);
                }
                catch (AlgorithmError anErr) {
                    System.out.println("Algorithm error restoring split length after reject in BrownianStateForTopologies.LogNormalEdgeProposal. This should not be possible.");
                }
            }
        }
     

    }
    
    /** Proposal on x0.
        Get next edge. Use RW on edge length with NNI. */
    public static class RWEdgeProposal extends Kernel {

        private NormalDistribution normDist;
        private int currentSplitInd;
        private double oldLength;
        private Split oldSplit, newSplit;

        /** Constructor  */
        public RWEdgeProposal(double sd) {
            normDist = new NormalDistribution(0, sd);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Proposal not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            normDist.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new RWEdgeProposal(normDist.getSigma());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            TreeAsSplits t = curr.getx0();

            ArrayList<Split> orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");
            
            oldSplit = orderedSplits.get(currentSplitInd);
            oldLength = t.getSplitLength(oldSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above

            // Now sample a new length for this edge
            double z = normDist.sample()+oldLength;
            try {
                if (z>=0.0) {
                    t.setSplitLength(oldSplit, z);
                    newSplit = oldSplit;
                }
                else {
                    Split[] nniSplits = OperationsOnTreeAsSplits.getNNISplits(t, oldSplit);
                    if (normDist.sample()>0.0) newSplit = nniSplits[0];
                    else newSplit = nniSplits[1];
                    t.remove(oldSplit);
                    t.add(newSplit, -z);
                    orderedSplits.set(currentSplitInd, newSplit);
                }
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error: unable to perform NNI in BrownianStateForTopologies diffusion proposal.");
            }
            
            currentSplitInd++;
            if (currentSplitInd==orderedSplits.size()) currentSplitInd = 0;

            return 0.0; // Symmetric
        }
            

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForTopologies curr = (BrownianStateForTopologies) theState.getSubState(subStateName);
            if(whichWay) {
                /* Nothing to do */
            }
            else {
                /* Restore x0 */
//                try {
                    curr.getx0().remove(newSplit);
                    curr.getx0().add(oldSplit, oldLength);
//                }
//                catch (AlgorithmError anErr) {
//                    System.out.println("Algorithm error restoring split length after reject in BrownianStateForTopologies.RWEdgeProposal. This should not be possible.");
//                }
            }
        }
     

    }


    /* PRIORS  -------------------------------------------------------------- */

    /** Prior based on combining indept priors on tree and variance param */
    static public class X0Prior implements DensityCalculator {

        /* Instance variables: a prior for the TreeState */
        private DensityCalculator treePrior;

        /* Constructors */

        /** Construct from a prior for TreeState  */
        public X0Prior(DensityCalculator t) {
            treePrior = t;
        }


        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            State subst = x.getSubState(subStateName);
            double l = treePrior.logDensity(x, subst.getSubState(TreeState.class).getName());
            return l*priorTemperature;
        }

        @Override
        public DensityCalculator makeCopy() {
            return new X0Prior(treePrior);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof BrownianState)) {
                throw new AlgorithmException("SimpleBrownianStatePrior not compatible with given state.");
            }
            treePrior.checkCompatibility(x,x.getSubState(subStateName).getSubState(TreeState.class).getName());
        }

    }

    
    /** Prior based on independent exponential on each edge length, scaled
     according to the (fixed) value of t_0 used. */
    static public class ExponentialEdgePrior implements DensityCalculator {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        protected double rate;
        protected boolean includePendants;

        public ExponentialEdgePrior(double mean, double sqrtt0, boolean ip) {
            rate = 1.0/(mean*sqrtt0);
            includePendants = ip;
        }

        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            State subst = x.getSubState(subStateName);
            TreeState xx = (TreeState) subst.getSubState(TreeState.class);
            double l = 0.0;
            Iterator<Map.Entry<Split,Double>> it = xx.getTree().getEntryIterator();
            while(it.hasNext()) {
                Map.Entry<Split, Double> theEntry = it.next();
                double length = theEntry.getValue();
                Split p = theEntry.getKey();
                if ((includePendants)||(p.isTerminal()==null)) {
                    l += ExponentialDistribution.logpdf(length, rate);
                }
            }
            return l*priorTemperature;
        }

        @Override
        public DensityCalculator makeCopy() {
            return new ExponentialEdgePrior(1.0/rate, 1.0, includePendants);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
                throw new AlgorithmException("Prior for branch lengths not compatible with given state.");
            }
        }

    } 

}