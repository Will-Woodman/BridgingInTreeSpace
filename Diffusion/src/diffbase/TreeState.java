/*
   TreeState
    Copyright (C) 2014  Tom M. W. Nye

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

package diffbase;

/**
 * MCMC state representing the source tree for a Brownian motion.
 */

import MCMC.State;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import simulation.NormalDistribution;
import simulation.LogNormalDistribution;
import simulation.GammaDistribution;
import cern.jet.random.tdouble.ChiSquare;
import treebase.Graph;
import treebase.Tree;
import treebase.TreeWithTopologicalOperations;
import treebase.AlgorithmException;
import java.util.Iterator;
import java.util.ArrayList;
import cern.jet.random.tdouble.DoubleUniform;
import java.util.Collections;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.Split;
import treebase.TreeAsSplits;
import treedatasets.OperationsOnTreeAsSplits;


public class TreeState extends State {

    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /** Instance variables */
    private TreeAsSplits theTree;
    


    /* CONSTRUCTORS --------------------------------------------------------- */

    public TreeState(TreeAsSplits t, String stateName) {
        super(stateName);
        theTree = t;
    }


    /* UTILITY METHODS ------------------------------------------------------ */

    public TreeAsSplits getTree() {
        return theTree;
    }
    
    public void setTree(TreeAsSplits x) {
        theTree = x;
    }

    @Override
    /** Print Newick String */
    public String getValueString() {
        return theTree.toString();
    }

    @Override
    public void setValueFromString(String str) throws AlgorithmException {
        theTree = new TreeAsSplits(new Tree(str));
    }
    @Override
    public int getLengthStringRepr() {
        return 1; // Single block of text
    }

    @Override
    public void updateRunningMean(State runningMean, int runningCount) throws AlgorithmException {
        throw new UnsupportedOperationException("Not supported.");
    }


    /* PRIORS  -------------------------------------------------------------- */

    /** Typical phylogenetic prior: uniform over topologies, independent gamma on edge length. */
    static public class StandardTreePrior implements DensityCalculator {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        protected double shape;
        protected double scale;
        protected boolean includePendants;

        public StandardTreePrior(double shape, double scale, boolean ip) {
            this.shape = shape;
            this.scale = scale;
            includePendants = ip;
        }

        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            TreeState xx = (TreeState) x.getSubState(subStateName);
            double l = 0.0;
            Iterator<Map.Entry<Split,Double>> it = xx.getTree().getEntryIterator();
            while(it.hasNext()) {
                Map.Entry<Split, Double> theEntry = it.next();
                double length = theEntry.getValue();
                Split p = theEntry.getKey();
                if ((includePendants)||(p.isTerminal()==null)) {
                    l += GammaDistribution.logpdf(length, shape, scale);
                }
            }
            return l*priorTemperature;
        }

        @Override
        public DensityCalculator makeCopy() {
            return new StandardTreePrior(shape, scale, includePendants);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof TreeState)) {
                throw new AlgorithmException("Prior for branch lengths not compatible with given state.");
            }
        }

    } // End standard tree prior


    /** Prior based on tree geometry. */
    static public class SphericalTreePrior implements DensityCalculator {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        protected double shape;
        protected double scale;
        protected boolean includePendants;
        protected int nTaxa;

        public SphericalTreePrior(double shape, double scale, boolean p, int numTaxa) {
            this.shape = shape;
            this.scale = scale;
            this.includePendants = p;
            this.nTaxa = numTaxa;
        }

        /** Default constructor: Gamma(1/2, rate=3.3175*4/N) where N is number of species; ensures 99% quantile is D_N^2 */
        public SphericalTreePrior(int numTaxa) {
            this.shape = 0.5;
            this.scale = numTaxa/13.27;
            this.includePendants = false;
            this.nTaxa = numTaxa;
        }
        public SphericalTreePrior(double gammaScale,int numTaxa) {
            this.shape = 0.5;
            this.scale = gammaScale;
            this.includePendants = false;
            this.nTaxa = numTaxa;
        }
        
        public String toString() {
            return "gamma prior on distance, shape="+String.format("%7.7f",shape)+" scale="+String.format("%7.7f",scale);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            TreeState xx = (TreeState) x.getSubState(subStateName);
            // Get distances from origin
            double distSqu=0.0, l;
            Iterator<Map.Entry<Split,Double>> it = xx.getTree().getEntryIterator();
            while(it.hasNext()) {
                Map.Entry<Split, Double> theEntry = it.next();
                l = theEntry.getValue();
                Split p = theEntry.getKey();
                if ((includePendants)||(p.isTerminal()==null)) {
                    distSqu += l*l;
                }
            }

            return priorTemperature*(GammaDistribution.logpdf(distSqu, shape, scale)-Math.log(distSqu)*(nTaxa-4)/2);
        }

        @Override
        public DensityCalculator makeCopy() {
            return new SphericalTreePrior(shape, scale, includePendants,nTaxa);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof TreeState)) {
                throw new AlgorithmException("Prior for branch lengths not compatible with given state.");
            }
        }

    } // End spherical tree prior

/** Prior based 95% of edge lengths less than 1 and finite density at the origin. */
    static public class NormalTreePrior implements DensityCalculator {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        protected double sd;
        protected boolean includePendants;

        public NormalTreePrior(double sd, boolean p) {
            this.sd = sd;
            this.includePendants = p;
        }

        /** Default constructor: sd=0.51 ensures 95% of edge lengths less than 1 */
        public NormalTreePrior() {
            this.sd = 0.51;
            this.includePendants = false;
        }
        public NormalTreePrior(double standDev) {
            this.sd = standDev;
            this.includePendants = false;
        }
        
        public String toString() {
            return "Reflected normal prior on edge lengths, sd="+String.format("%7.7f",sd);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            TreeState xx = (TreeState) x.getSubState(subStateName);
            // Get distances from origin
            double theLogpdf=0.0, l;
            Iterator<Map.Entry<Split,Double>> it = xx.getTree().getEntryIterator();
            while(it.hasNext()) {
                Map.Entry<Split, Double> theEntry = it.next();
                l = theEntry.getValue();
                Split p = theEntry.getKey();
                if ((includePendants)||(p.isTerminal()==null)) {
                    theLogpdf += Math.log(2)+priorTemperature*NormalDistribution.logpdf(l, 0, sd);
                }
            }
            //line to check it gives the correct density:
            //System.out.println(theLogpdf);
            return theLogpdf;
           
        }

        @Override
        public DensityCalculator makeCopy() {
            return new NormalTreePrior(sd, includePendants);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof TreeState)) {
                throw new AlgorithmException("Prior for branch lengths not compatible with given state.");
            }
        }

    } // End normal tree prior   
    
    static public class ChiSquTreePrior implements DensityCalculator {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        //sd is the standard deviation of the normal distribution that we are summing over, so normally 0.51 as 
        //in normal prior
        protected double sd;
        protected boolean includePendants;
        

        public ChiSquTreePrior(double sd, boolean p) {
            this.sd = sd;
            this.includePendants = p;
        }

        /** Default constructor: sd=0.51 ensures 95% of edge lengths less than 1 */
        public ChiSquTreePrior() {
            this.sd = 0.51;
            this.includePendants = false;
        }
        public ChiSquTreePrior(double standDev) {
            this.sd = standDev;
            this.includePendants = false;
        }
        
        public String toString() {
            return "Reflected normal prior on edge lengths, sd="+String.format("%7.7f",sd);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            TreeState xx = (TreeState) x.getSubState(subStateName);
            ChiSquare theChiSqu = new ChiSquare(xx.getTree().getNumTaxa()-3,Random.getEngine());
            // Get distances from origin
            double distSqu=0.0, l;
            Iterator<Map.Entry<Split,Double>> it = xx.getTree().getEntryIterator();
            while(it.hasNext()) {
                Map.Entry<Split, Double> theEntry = it.next();
                l = theEntry.getValue();
                Split p = theEntry.getKey();
                if ((includePendants)||(p.isTerminal()==null)) {
                    distSqu += l*l;
                }
            }
            //line to check it gives the correct density:
            //System.out.println(theLogpdf);
            return Math.log(theChiSqu.pdf(distSqu/(sd*sd)));
           
        }

        @Override
        public DensityCalculator makeCopy() {
            return new NormalTreePrior(sd, includePendants);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof TreeState)) {
                throw new AlgorithmException("Prior for branch lengths not compatible with given state.");
            }
        }

    } // End n
    
    

    /* PROPOSALS  ----------------------------------------------------------- */


    /** Proposal to change one edge length at a time,
     with NNI if the proposed new length is negative. */
    public static class SingleStepRWProposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        private NormalDistribution norm;
        private int currentSplitInd;
        private ArrayList<Split> orderedSplits;  
        private Split oldSplit, proposedSplit;
        private double oldLength;

        public SingleStepRWProposal(double normSD) {
            name = "Single step Gaussian RW proposal";
            norm = new NormalDistribution(0.0, normSD);
            currentSplitInd = 0;
        }
        protected SingleStepRWProposal(NormalDistribution n) {
            norm = n;
            currentSplitInd = 0;
        }

        public double getSD() {
           return norm.getSigma();
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof TreeState)) throw new AlgorithmException("NNI proposal not compatible with given state.");
        }

        @Override
        public void resetRandomEngineSeed() {
            norm.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new SingleStepRWProposal(norm.getSigma());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            TreeState curr = (TreeState) theState.getSubState(subStateName);

            TreeAsSplits t = curr.getTree();
            orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");

            Split theSplit = orderedSplits.get(currentSplitInd);
            
            oldSplit = theSplit;
            oldLength = t.getSplitLength(theSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above

            // Now sample a new length for this edge
            double z = norm.sample()+oldLength;
            try {
                if (z>=0.0) {
                    t.setSplitLength(theSplit, z);
                    proposedSplit = oldSplit;
                }
                else {
                    Split[] nniSplits = OperationsOnTreeAsSplits.getNNISplits(t, theSplit);
                    if (norm.sample()>0.0) proposedSplit = nniSplits[0];
                    else proposedSplit = nniSplits[1];
                    t.remove(theSplit);
                    t.add(proposedSplit, -z);
                    orderedSplits.set(currentSplitInd, proposedSplit);
                }
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error: unable to perform NNI in TreeState diffusion proposal.");
            }
            
            currentSplitInd++;
            if (currentSplitInd==orderedSplits.size()) currentSplitInd = 0;

            return 0.0; // Symmetric
        }

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            TreeState xt = (TreeState) theState.getSubState(subStateName);
            TreeAsSplits t = xt.getTree();
            if(whichWay) {
                // Nothing to do
            }
            else {
                // Reverse the NNI if nec
                if (oldSplit!=proposedSplit) {
                    t.remove(proposedSplit);
                    int ind = orderedSplits.indexOf(proposedSplit);
                    orderedSplits.set(ind, oldSplit);
                }
                t.add(oldSplit, oldLength);
            }
        }
    }


    /** Proposal to change one edge length at a time via log normal,
     with no changes in topology. 
     Every time proposal is called, move on to next internal edge -- call proposal
     repeatedly in a loop / sweep. */
    public static class EdgeLengthProposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        private LogNormalDistribution sampler;
        private int currentSplitInd;
        private double oldLength;
        private double sigma;
        private Split theSplit;

        public EdgeLengthProposal(double sd) {
            sampler = new LogNormalDistribution(0.0, 1.0);
            name = "Edge length proposal";
            currentSplitInd = 0;
            sigma = sd;
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof TreeState)) throw new AlgorithmException("Edge length proposal not compatible with given state.");
        }

        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new EdgeLengthProposal(sigma);
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            TreeState curr = (TreeState) theState.getSubState(subStateName);
            TreeAsSplits t = curr.getTree();

            ArrayList<Split> orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");
            
            theSplit = orderedSplits.get(currentSplitInd);
            oldLength = t.getSplitLength(theSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above

            // Now sample a new length for this edge
            double newLength = sampler.sample(Math.log(oldLength), sigma);
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
            TreeState xt = (TreeState) theState.getSubState(subStateName);
            if(whichWay) {
                // Nothing to do
            }
            else {
                try {
                    xt.getTree().setSplitLength(theSplit, oldLength);
                } catch (AlgorithmError ex) {
                    System.out.println("Error restoring split length following reject in edge length proposal. This should not be possible.");
                }
            }

        }

    }


    /** Proposal to change one edge at a time via NNI move and Gamma proposal for new length.
     Call repeatedly within a sweep. */
    public static class NNIWithGammaProposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        private GammaDistribution sampler;
        private DoubleUniform unif;
        private int currentSplitInd;
        private ArrayList<Split> orderedSplits;  
        private Split oldSplit, proposedSplit;
        private double oldLength;


        public NNIWithGammaProposal(double shape, double scale) {
            name = "NNI with Gamma branch length proposal";
            sampler = new GammaDistribution(shape, scale);
            currentSplitInd = 0;
            unif = new DoubleUniform(Random.getEngine());
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof TreeState)) throw new AlgorithmException("NNI proposal not compatible with given state.");
        }

        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new NNIWithGammaProposal(sampler.getShape(), sampler.getScale());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            
            TreeState curr = (TreeState) theState.getSubState(subStateName);
            TreeAsSplits t = curr.getTree();
            orderedSplits = new ArrayList();
            orderedSplits.addAll(t.getNonTrivialSplits());
            Collections.sort(orderedSplits);
            if (!t.fullyResolved()) System.out.println("Warning: initial tree for a diffusion is not fully resolved.");

            Split theSplit = orderedSplits.get(currentSplitInd);
            
            oldSplit = theSplit;
            oldLength = t.getSplitLength(theSplit);
            // NB: this will completely hang if you pass in a star tree!!! -- see warning above
            
            // Now sample a new length and replacement for this edge
            Split[] nniSplits = OperationsOnTreeAsSplits.getNNISplits(t, theSplit);
            double newLength = sampler.sample();
            if (unif.nextBoolean()) proposedSplit = nniSplits[0];
            else proposedSplit = nniSplits[1];

            t.remove(theSplit);
            t.add(proposedSplit, newLength);
            orderedSplits.set(currentSplitInd, proposedSplit);
            
            currentSplitInd++;
            if (currentSplitInd==orderedSplits.size()) currentSplitInd = 0;

            return (sampler.getShape()-1.0)*(Math.log(oldLength)-Math.log(newLength))-(oldLength-newLength)/sampler.getScale(); // Ratio for gamma independence sampler
        }

        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            TreeState xt = (TreeState) theState.getSubState(subStateName);
            TreeAsSplits t = xt.getTree();
            if(whichWay) {
                // Nothing to do
            }
            else {
                // Reverse the NNI
                t.remove(proposedSplit);
                t.add(oldSplit, oldLength);
                int ind = orderedSplits.indexOf(proposedSplit);
                orderedSplits.set(ind, oldSplit);
            }

        }
        
        public static void main(String[] args) throws AlgorithmException {
            //test normal tree prior
            String testTreeString = "(((C:1.0000000,B:1.0000000):0.5,A:1.0000000):0.5,E:1.0000000,D:1.0000000);";
            Tree testTree= new Tree(testTreeString);
            //but now I need to create a global state etc right? Could I just print it out in the MCMC I am already doing?
        }

    }

}
