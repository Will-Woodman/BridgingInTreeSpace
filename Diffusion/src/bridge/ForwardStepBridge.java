/*
   ForwardStepBridge
    Copyright (C) 2013  Tom M. W. Nye

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
import cern.jet.random.tdouble.ChiSquare;
import geodesics.Geodesic;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.NormalDistribution;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
    Bridge by stepping forward between two trees. 
    Include a penalty term for every kink in the geodesic
 */

public abstract class ForwardStepBridge {
    
    /* The bridge itself -- trees and the number of steps */
    protected final int numSteps;
    protected TreeAsSplits[] theTrees;
    protected TreeAsSplits treeA, treeB;

    /* Storage of numerical values */
    protected double totalLogLike;
    protected double[] logLike;
    protected double[] logStepPropDensity;
    
    //Added two extra fields for efficient updating of the likelihoods when updating t0 value
    protected double[] stepLengths;
    protected double[] boundaryCrosses;

    protected NormalDistribution normDist;
    protected DoubleUniform unifDist;
    protected ChiSquare chiSquaredDist;
    
       
    
    /** Constructor */
    public ForwardStepBridge(TreeAsSplits tA, TreeAsSplits tB, int n) {
        treeA = (tA!=null) ? tA.efficientClone() : null;
        treeB = (tB!=null) ? tB.efficientClone() : null;
        numSteps = n;
        theTrees = new TreeAsSplits[numSteps+1];
        theTrees[0] = treeA;
        theTrees[numSteps] = treeB;

        /* Set up log like and log prop densities */
        totalLogLike = Double.NaN;
        logLike = new double[numSteps+1];
        logStepPropDensity = new double[numSteps+1]; // NB: final step to treeB might have a log prop density
        //Auxilliary information for fast update of likelihoods for new t0
        stepLengths=new double[numSteps+1];
        boundaryCrosses= new double[numSteps+1];
        
        /* Random number generators */
        normDist = new NormalDistribution(0.0, 1.0);
        unifDist = new DoubleUniform(Random.getEngine());
        
    }
    
    /** Constructor from an array list of trees as splits - using in the MCMC for the quasistatic stepping stone sampling*/
    public ForwardStepBridge(TreeAsSplits[] startTrees) {
        //add all the trees
        numSteps = startTrees.length-1;
        treeA = startTrees[0].efficientClone();
        treeB = startTrees[numSteps].efficientClone();
        theTrees = new TreeAsSplits[numSteps+1];
        theTrees[0]=treeA;
        theTrees[numSteps]=treeB;
        for(int i=1;i<numSteps;i++){
        theTrees[i] = startTrees[i].efficientClone();
        }

        /* Set up log like and log prop densities */
        totalLogLike = Double.NaN;
        logLike = new double[numSteps+1];
        logStepPropDensity = new double[numSteps+1]; // NB: final step to treeB might have a log prop density
        //Auxilliary information for fast update of likelihoods for new t0
        stepLengths=new double[numSteps+1];
        boundaryCrosses= new double[numSteps+1];
        
        /* Random number generators */
        normDist = new NormalDistribution(0.0, 1.0);
        unifDist = new DoubleUniform(Random.getEngine());
        
    }
    
    /* Abstract methods ----------------------------------------------------- */
    
    /** Compute the log likelihood and log prop density associated with a particular step. */
    protected abstract double[] calcStepLogLike(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException;

    
    /** Abstract method: copy the current bridge (end points and number of steps)
     for use in MCMC proposals. */
    public abstract ForwardStepBridge makeEmptyCopy();
    public abstract ForwardStepBridge makeEmptyCopy(TreeAsSplits newTree, boolean whichEnd);
    public abstract ForwardStepBridge makeFreshBridge(TreeAsSplits tA, TreeAsSplits tB);

    
    
    /* Util methods --------------------------------------------------------- */
    
    /** Compute likelihood */
    public double sumLogLike() {
        double l = 0.0;
        for (int i=0; i<numSteps+1; i++) {
            l += logLike[i];
        }
        return l;
    }

    public double sumStepLogPropDensity() {
        double p = 0.0;
        for (int i=0; i<numSteps+1; i++) {
            p += logStepPropDensity[i];
        }
        return p;
    }
    
    public double sumLogLike(int a, int b) {
        double l = 0.0;
        for (int i=a; i<b; i++) {
            l += logLike[i];
        }
        return l;
    }

    public double sumStepLogPropDensity(int a, int b) {
        double p = 0.0;
        for (int i=a; i<b; i++) {
            p += logStepPropDensity[i];
        }
        return p;
    }
    
        
    public double getTotalLogLike() {
        
        
        if (Double.isNaN(totalLogLike)) {
            totalLogLike = sumLogLike();
        }
                
        return totalLogLike;
       
    }
    
    /*the below is for stepping stone marginal likelihood when beta=0 -- since then we need
    to compute the GGF log likelihood from scratch to calculate the estimate
    It can also be used in the CHib and tunnel estimators when the likelihood for the bridge is the independence proposal
    likelihood
    */
    public double getTotalGGFLogLike(double t0) throws AlgorithmError, AlgorithmException{
        double GGFLogLike=0;
        
        for(int i=0;i<numSteps;i++){
            Geodesic g = new Geodesic(theTrees[i], theTrees[i+1]);
            if(!g.isSimple()) System.out.println("Non simple geodesic somewhere");
            double[] GGFStepLogLike=calcStepLogLikeGGF(theTrees[i], theTrees[i+1], g, t0);
            GGFLogLike+=GGFStepLogLike[0];
        }
        return(GGFLogLike);
    }
    
        protected double[] calcStepLogLikeGGF(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException {
        int nprime = tA.getNumTaxa()-3;
        double likeStepSD = Math.sqrt(t0/numSteps);
        if (!g.isSimple()) {
            throw new AlgorithmException("Non-simple segment in BridgeWithApproxMVNLike");
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
                    //GGF Likelihood
                    ll += NormalDistribution.logpdf(la+lb, 0.0, likeStepSD)-Math.log(2.0);
                }
                else {
                    // Shared split
                    double lb = tB.getSplitLength(p);
                    //GGF Likelihood
                     ll += NormalDistribution.logpdf(la-lb, 0.0, likeStepSD);
                } 
            }
        }

        return new double[] {ll, 0.0};
        
    }
            
    
    /** Given Two trees compute number of boundaries traversed by geodesic. */
    public static int countBoundaries(TreeAsSplits a, TreeAsSplits b) {
        /* Assume both trees fully resolved -- and geodesic a,b is simple! */
        if ((a.getNumSplits()!=2*a.getNumTaxa()-3)||(b.getNumSplits()!=2*b.getNumTaxa()-3)) {
            System.out.println("Warning: unresolved tree when counting boundaries.");
        }

        Iterator<Split> it = a.getSplitIterator();
        int count = 0;
        while (it.hasNext()) {
            Split p = it.next();
            if (!b.contains(p)) {
                // Increase counter only if the split in a is away from a boundary
                if (a.getSplitLength(p)>0.0) count++;
            }
        }
        return count;
    }
    
    
    public int getNumSteps() {
        return numSteps;
    }
    
    public double[] getBoundaryCrosses(){
        return(boundaryCrosses);
    }
    
    
    /** Return array list of trees along the bridge. */
    public ArrayList<TreeAsSplits> getTreePathList() {
        ArrayList<TreeAsSplits> res = new ArrayList();
        res.addAll(Arrays.asList(theTrees));
        return res;
    }
    public TreeAsSplits[] getTreePathArray() {
        TreeAsSplits[] res = Arrays.copyOf(theTrees, theTrees.length);
        return res;
    }

    
    /** copyTo and setFrom used in storing and restoring during MCMC */
    public void copyTo(ForwardStepBridge b) {
        b.setFrom(theTrees, logLike, logStepPropDensity,stepLengths,boundaryCrosses);
    } 
     private void setFrom(TreeAsSplits[] t, double[] l, double[] p,double[] sl, double[] bc) {
        theTrees = Arrays.copyOf(t, numSteps+1);
        treeA = t[0]; treeB = t[numSteps];
        logLike = Arrays.copyOf(l, numSteps+1);
        logStepPropDensity = Arrays.copyOf(p, numSteps+1);
        totalLogLike = Double.NaN;
        stepLengths = Arrays.copyOf(sl, numSteps+1);
        boundaryCrosses = Arrays.copyOf(bc, numSteps+1);
    }
    
    
    /** Sum squ distances DEBUG! */
    public double sumSquaredDistances() {
        double dd = 0.0;
        for (int k=1; k<=numSteps; k++) {
            try {
                Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
                double d = g.getInternalLength();
                dd += d*d;
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error summing squared lengths for a bridge. \nError report --"+anErr.getMessage());
            }
             
        }
        return dd;
    }
    
    //adding the below for t0 joint proposal -- may not be needed.
    public double sumDistances() {
        double dd = 0.0;
        for (int k=1; k<=numSteps; k++) {
            try {
                Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
                double d = g.getInternalLength();
                dd += d;
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error summing squared lengths for a bridge. \nError report --"+anErr.getMessage());
            }
             
        }
        return dd;
    }
    
    /* Classes and methods for a single step forward -------------------------*/
    
    /** private class for returning subdivision results */
    private static class StepResult {
        public TreeAsSplits tree;
        public double logDensity;
        
        public StepResult(TreeAsSplits x, double l) {
            tree =  x;
            logDensity = l;
        }     
    }
    
    
    /** Utility method for subdividing any given interval. 
     The returned log density includes perturbation, traversing boundaries, and resolving. */
    static protected StepResult stepForwardViaMVN(TreeAsSplits tA, TreeAsSplits tB, double sd, int numRemainingSteps, NormalDistribution normDist, DoubleUniform unifDist) {
        
        GeodesicForBridging g = null;
        try {
            g = new GeodesicForBridging(tA, tB);
        } catch (AlgorithmError ex) {
            System.out.println("Error subdividing in Bridge");
        }
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        
        /* If penalty != 0 you might like to "boost" the sd of the perturbation below. */
//        if (penalty[0]>0) {
//            sd = sd*2;
//        }
        
        double resLogDens = 0.0; // Log density
        if (!nextTree.fullyResolved()) {
            try {
                resLogDens = TreeResolver.resolveTree(nextTree, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        /* Perturb the tree */
        StepResult info = null;
        try {        
            info = sampleDirectMVN(nextTree, sd, normDist, unifDist);
        } catch (AlgorithmException ex) {
            System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
        }
        info.logDensity += resLogDens; // Add on the log density for resolving        

        return info;
        
    }

    
     
    
    /** Compute step density. */
    public static double computeMVNStepDensity(TreeAsSplits treeA, TreeAsSplits treeB, TreeAsSplits realization, double sigma, int numRemainingSteps, DoubleUniform unifDist) throws AlgorithmException {
        
        double logDensity = 0.0;
        int nprime = treeA.getNumTaxa()-3;
        
        GeodesicForBridging g = new GeodesicForBridging(treeA, treeB);
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        nextTree.removeZeroLengthSplits();
        
        if (!nextTree.fullyResolved()) {
            try {
                logDensity += TreeResolver.resolveTree(nextTree, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        Geodesic h = new Geodesic(nextTree, realization);
        if (!h.isSimple()) {
            throw new AlgorithmException("Non-simple geodesic in density calculation for sampled tree.");
        }
        
        /* Measure distance */
        double len = h.getInternalLength();
        logDensity += -nprime*Math.log(sigma)-0.5*len*len/(sigma*sigma); //  multivariate normal density
        
        /* Take account of boundaries */
        int nb = countBoundaries(nextTree, realization);
        logDensity += -0.6931472*nb;

        return logDensity;
    }
    
    
    
    
    /** Log density for the central gaussian  */
    private static double logCentralDensity(TreeAsSplits x, double sd) {
        int nprime = x.getNumTaxa()-3;
        double l = 0.0;
        for (int k=1; k<=nprime-2; k++) {
            l -= Math.log(2*k-1);
        }
        l += nprime*Math.log(2.0);
        l -= nprime*Math.log(sd);
        l -= 0.5*x.sumSquaredLengths(true)/sd/sd;
        return l;
    }
    
    /** Sample a central gaussian */
    public static TreeAsSplits sampleCentralGaussian(double sd, NormalDistribution normDist, DoubleUniform unifDist) {
        StarTree s = StarTree.getInstance();
        TreeAsSplits theTree = s.getTreeAsSplits();
        try {
            double l = TreeResolver.resolveTree(theTree, s.getTree(), unifDist);
        } catch (AlgorithmException ex) {
            System.out.println("Error resolving star tree.");
        }
        try {
            Iterator<Split> it = theTree.getSplits().iterator();
            while (it.hasNext()) {
                Split p = it.next();
                if (p.isTerminal()==null) {
                    theTree.setSplitLength(p, Math.abs(normDist.sample(0.0, sd)));
                }
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error setting edge length when sampling central Gaussian.");
        }
        
        return theTree;
    }
    
    
    /** Return weight of central gaussian and its SD given current tree being perturbed */
     static protected double[] calcMixtureWeight(TreeAsSplits x, double sd) {
        
        /* weight = 1 - chisqu of sum of squ edge lengths with a max value of 0.999*/
        double dd = x.sumSquaredLengths(true);
        
        int nprime = x.getNumTaxa()-3;
        double cdf = NormalDistribution.cdf(dd, sd*sd*nprime, Math.sqrt(2.0*nprime)*sd*sd); // This is a normal approximation to chi squared on nprime df
        if(cdf>=0.999) cdf=0.999;
        return new double[] {(1.0-cdf), sd};

    }
    

    
    /* The following is in case a different mixture is to be used for the independence proposal distribution in marginal likelihood
     calcs (currently it uses standard mixture)
     */
    static protected StepResult stepForwardViaFixedMixture(TreeAsSplits tA, TreeAsSplits tB, double sd, double sdRW, int numRemainingSteps, NormalDistribution normDist, DoubleUniform unifDist) {
        
        
        GeodesicForBridging g = null;
        try {
            g = new GeodesicForBridging(tA, tB);
        } catch (AlgorithmError ex) {
            System.out.println("Error subdividing in Bridge");
        }
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        nextTree.removeZeroLengthSplits();
        double resLogDens=0.0;
        if (!nextTree.fullyResolved()) {
            try {
                double nS = nextTree.getNumInternalSplits();//number of splits in the next tree
                resLogDens = TreeResolver.resolveTree(nextTree, unifDist);
                resLogDens+=(tA.getNumTaxa()-3-nS)*Math.log(2);

            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        /* We have the tree to perturb.
           Now compute mixture weights and SD of central gaussian */
        double[] weightAndSD = calcMixtureWeight(nextTree, sd);
        //double[] weightAndSD = {0.01,sd};
              
        //centralDens here refers to the part coming from the GGF(y_{j-1},sd^2)
        double centralDens=0.0, mvnDens=0.0;
        StepResult info = null;
        
        double w1, w2;
        if (unifDist.nextDouble()>weightAndSD[0]) {
            /* Sample from the point on the geodesic */
            try {        
                info = sampleDirectMVN(nextTree, sd, normDist, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
            }
            mvnDens = info.logDensity+resLogDens; // Add on the log density for resolving

            /* Get central density --code basically stolen from the else part */
            Geodesic h=null;
            try {
                h = new Geodesic(tA, info.tree);
            } catch (AlgorithmError ex) {
                System.out.println("Error making geodesic when sampling mixture.");
            }
            if (h.isSimple()) {
                w2 = 1.0-weightAndSD[0];
                /* Measure distance */
                double len = h.getInternalLength();
                centralDens = -(tA.getNumTaxa()-3)*Math.log(sdRW)-0.5*len*len/(sdRW*sdRW); //  multivariate normal density
                //centralDens+= resLogDens;

                /* Take account of boundaries */
                int nb = countBoundaries(tA, info.tree);
                centralDens += -0.6931472*nb;
            }
            /*
            if h is not simple then we exit the algorithm for building the bridge or it has zero likelihood
            so nothing to take care of
            */
            w1 = weightAndSD[0];
            w2 = 1.0-w1;
            
        }
        else {
            /* Sample MVN from y_{j-1} -- which is tA */
            try {        
                //sampleDirectMVN  modifies the input tree so need to make a copy of tA to input
                TreeAsSplits treeA=tA.clone();
                info = sampleDirectMVN(treeA, sdRW, normDist, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
            }
            centralDens = info.logDensity; // Don't need to resolve - y_{j-1} should be fully resolved
            TreeAsSplits centralTree=info.tree;
            w1 = weightAndSD[0];
            w2 = 0.0; // Point cannot have come from MVN, unless geod is found to be simple
            
            /* Compute density as if central tree came from MVN from nextTree */
            Geodesic h=null;
            try {
                h = new Geodesic(nextTree, centralTree);
            } catch (AlgorithmError ex) {
                System.out.println("Error making geodesic when sampling mixture.");
            }
            if (h.isSimple()) {
                w2 = 1.0-weightAndSD[0];
                /* Measure distance */
                double len = h.getInternalLength();
                mvnDens = -(centralTree.getNumTaxa()-3)*Math.log(sd)-0.5*len*len/(sd*sd); //  multivariate normal density
                mvnDens+= resLogDens;

                /* Take account of boundaries */
                int nb = countBoundaries(nextTree, centralTree);
                mvnDens += -0.6931472*nb;
                info.tree=centralTree;
            }
            else {
                info.logDensity = Math.log(w1)+centralDens;
                info.tree=centralTree;
                return info;
            }
 
        }
        
        /* Take care of numerical stability */
        double u1 = Math.log(w1);
        double u2 = Math.log(w2);
        if ((u1+centralDens)>(u2+mvnDens)) {
            info.logDensity = u1+centralDens;
            info.logDensity += Math.log(1+Math.exp(u2-u1+mvnDens-centralDens));
        }
        else {
            info.logDensity = u2+mvnDens;
            info.logDensity += Math.log(1+Math.exp(u1-u2-mvnDens+centralDens));
        }
        
        return info;  
    }
    
        // make a mixture step that uses fixed weights and generates the mixture part from y_{j-1}
     // more efficient version of stepForwardViaFixedMixture for computing prop simple in the marginal likelihood methods
    // here we only need to know if the bridge has simple steps, not its proposal density or log likelihood
    static protected StepResult stepForwardViaFixedMixturePropSimple(TreeAsSplits tA, TreeAsSplits tB, double sd, double sdRW, int numRemainingSteps, NormalDistribution normDist, DoubleUniform unifDist) {
        
        
        GeodesicForBridging g = null;
        try {
            g = new GeodesicForBridging(tA, tB);
        } catch (AlgorithmError ex) {
            System.out.println("Error subdividing in Bridge");
        }
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        nextTree.removeZeroLengthSplits();
        double resLogDens=0.0;
        if (!nextTree.fullyResolved()) {
            try {
                resLogDens = TreeResolver.resolveTree(nextTree, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        /* We have the tree to perturb.
           Now compute mixture weights and SD of central gaussian */
        //double[] weightAndSD = {0.01,sd};
        double[] weightAndSD = calcMixtureWeight(nextTree, sd);
              
        //centralDens here refers to the part coming from the GGF(y_{j-1},sd^2)
        double centralDens=0.0, mvnDens=0.0;
        StepResult info = null;
        
        double w1, w2;
        if (unifDist.nextDouble()>weightAndSD[0]) {
            /* Sample MVN */
            try {        
                info = sampleDirectMVN(nextTree, sd, normDist, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
            }
        }
        else {
            /* Sample MVN from y_{j-1} -- which is tA */
            try {        
                //sampleDirectMVN  modifies the input tree so need to make a copy of tA to input
                TreeAsSplits treeA=tA.clone();
                info = sampleDirectMVN(treeA, sd, normDist, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
            }
            
 
        }
        return info;  
    }
    
    /** Compute step density when using fixed weights and generating the mixture part from y_{j-1}. */
    public static double computeFixedMixtureStepDensity(TreeAsSplits treeA, TreeAsSplits treeB, TreeAsSplits realization, double sigma, double sdRW, int numRemainingSteps, DoubleUniform unifDist)  {
        GeodesicForBridging g = null;
        try {
            g = new GeodesicForBridging(treeA, treeB);
        } catch (AlgorithmError ex) {
            System.out.println("Error subdividing in Bridge");
        }
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        nextTree.removeZeroLengthSplits();
        double resLogDens=0.0;
        if (!nextTree.fullyResolved()) {
            try {
                double nS = nextTree.getNumInternalSplits();//number of splits in the next tree
                resLogDens = TreeResolver.resolveTree(nextTree, unifDist);
                resLogDens+=(treeA.getNumTaxa()-3-nS)*Math.log(2);
            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        //double[] weightAndSD = {0.01,sigma};
        double[] weightAndSD = calcMixtureWeight(nextTree, sigma);
        double centralDens=0.0, mvnDens=0.0;
        
        //geodesic between next tree and realisation -- main part
        Geodesic h=null;
        try {
            h = new Geodesic(nextTree, realization);
        } catch (AlgorithmError ex) {
            System.out.println("Error making geodesic when computing mixture density.");
        }
        //geodesic between y_{j-1} and realisation -- 'central' part
        Geodesic mh=null;
        try {
            mh = new Geodesic(treeA, realization);
        } catch (AlgorithmError ex) {
            System.out.println("Error making geodesic when computing mixture density.");
        }
        
        
        
        //calculate the density of the part coming from GGF(y_{j-1})
        
        /* Compute mixture log density */
            double mlen = mh.getInternalLength();
            centralDens = -(treeA.getNumTaxa()-3)*Math.log(sdRW)-0.5*mlen*mlen/(sdRW*sdRW); //  multivariate normal density
        
            /* Take account of boundaries */
            int mnb = countBoundaries(treeA, realization);
            centralDens += -0.6931472*mnb;
            
            /* Take account of resolving --shouldn't need to do this, again check? */
            //centralDens += resLogDens;
        
        
        if (h.isSimple()) {

            /* Compute mixture log density */
            double len = h.getInternalLength();
            mvnDens = -(nextTree.getNumTaxa()-3)*Math.log(sigma)-0.5*len*len/(sigma*sigma); //  multivariate normal density
        
            /* Take account of boundaries */

            int nb = countBoundaries(nextTree, realization);
            mvnDens += -0.6931472*nb;
            
            /* Take account of resolving */
            mvnDens += resLogDens;
            
            /* Take care of numerical stability */
            double logDensity;
            double u1 = Math.log(weightAndSD[0]);
            double u2 = Math.log(1.0-weightAndSD[0]);
            if ((u1+centralDens)>(u2+mvnDens)) {
                logDensity = u1+centralDens;
                logDensity += Math.log(1+Math.exp(u2-u1+mvnDens-centralDens));
            }
            else {
                logDensity = u2+mvnDens;
                logDensity += Math.log(1+Math.exp(u1-u2-mvnDens+centralDens));
            }
            return logDensity;
        }
        else {
            return Math.log(weightAndSD[0])+centralDens;
        }
    }
    

    // make a mixture step that uses variable weights and generates the mixture part from y_{j-1}
    static protected StepResult stepForwardViaMixture(TreeAsSplits tA, TreeAsSplits tB, double sd, double sdRW, int numRemainingSteps, NormalDistribution normDist, DoubleUniform unifDist) {
        
        
        GeodesicForBridging g = null;
        try {
            g = new GeodesicForBridging(tA, tB);
        } catch (AlgorithmError ex) {
            System.out.println("Error subdividing in Bridge");
        }
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        
        
        nextTree.removeZeroLengthSplits();
        //System.out.println(nextTree.toString());
        double resLogDens=0.0;
        if (!nextTree.fullyResolved()) {
            try {
                double nS = nextTree.getNumInternalSplits();//number of splits in the next tree
                resLogDens = TreeResolver.resolveTree(nextTree, unifDist);
                resLogDens+=(tA.getNumTaxa()-3-nS)*Math.log(2);
                //System.out.println(nextTree.getTree().toTopologyString());
            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        /* We have the tree to perturb.
           Now compute mixture weights and SD of central gaussian */
        double[] weightAndSD = calcMixtureWeight(nextTree, sd);
              
        //centralDens here refers to the part coming from the GGF(y_{j-1},sd^2)
        double centralDens=0.0, mvnDens=0.0;
        StepResult info = null;
        
        double w1, w2;
        if (unifDist.nextDouble()>weightAndSD[0]) {
            /* Sample MVN */
            try { 
                info = sampleDirectMVN(nextTree, sd, normDist, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
            }
            mvnDens = info.logDensity+resLogDens; // Add on the log density for resolving

            /* Get central density -- modified from the else part */
            Geodesic h=null;
            try {
                h = new Geodesic(tA, info.tree);
            } catch (AlgorithmError ex) {
                System.out.println("Error making geodesic when sampling mixture.");
            }
            if (h.isSimple()) {
                w2 = 1.0-weightAndSD[0];
                /* Measure distance */
                double len = h.getInternalLength();
                centralDens = -(tA.getNumTaxa()-3)*Math.log(sdRW)-0.5*len*len/(sdRW*sdRW); //  multivariate normal density

                /* Take account of boundaries */
                int nb = countBoundaries(tA, info.tree);
                centralDens += -0.6931472*nb;
            }                //test:
            
            /* Set weights */
            w1 = weightAndSD[0];
            w2 = 1.0-w1;
            
        }
        else {
            /* Sample MVN from y_{j-1} -- which is tA */
            try {        
                //sampleDirectMVN  modifies the input tree so need to make a copy of tA to input
                TreeAsSplits treeA=tA.clone();
                info = sampleDirectMVN(treeA, sdRW, normDist, unifDist);
            } catch (AlgorithmException ex) {
                System.out.println("AlgorithmError sampling direct MVN. "+ex.getMessage());
            }
            centralDens = info.logDensity; // Don't need to resolve - y_{j-1} should be fully resolved
            TreeAsSplits centralTree=info.tree;

            w1 = weightAndSD[0];
            w2 = 0.0; // Point cannot have come from MVN, unless geod is found to be simple
            
            /* Compute density as if central tree came from MVN from nextTree */
            Geodesic h=null;
            try {
                h = new Geodesic(nextTree, centralTree);
            } catch (AlgorithmError ex) {
                System.out.println("Error making geodesic when sampling mixture.");
            }
            if (h.isSimple()) {
                w2 = 1.0-weightAndSD[0];
                /* Measure distance */
                double len = h.getInternalLength();
                mvnDens = -(centralTree.getNumTaxa()-3)*Math.log(sd)-0.5*len*len/(sd*sd); //  multivariate normal density
                mvnDens+= resLogDens;

                /* Take account of boundaries */
                int nb = countBoundaries(nextTree, centralTree);
                mvnDens += -0.6931472*nb;
                info.tree=centralTree;
            }
            else {
                info.logDensity = Math.log(w1)+centralDens;
                info.tree=centralTree;
                return info;
            }
 
        }
        
        /* Take care of numerical stability */
        double u1 = Math.log(w1);
        double u2 = Math.log(w2);
        if ((u1+centralDens)>(u2+mvnDens)) {
            info.logDensity = u1+centralDens;
            info.logDensity += Math.log(1+Math.exp(u2-u1+mvnDens-centralDens));
        }
        else {
            info.logDensity = u2+mvnDens;
            info.logDensity += Math.log(1+Math.exp(u1-u2-mvnDens+centralDens));
        }
        
        return info;  
    }
    
    /** Compute step density when using variable weights and generating the mixture part from y_{j-1}. */
    public static double computeMixtureStepDensity(TreeAsSplits treeA, TreeAsSplits treeB, TreeAsSplits realization, double sigma, double sdRW, int numRemainingSteps, DoubleUniform unifDist)  {
        GeodesicForBridging g = null;
        try {
            g = new GeodesicForBridging(treeA, treeB);
        } catch (AlgorithmError ex) {
            System.out.println("Error subdividing in Bridge");
        }
        
        int[] penalty = new int[1];
        TreeAsSplits nextTree = getStepLocation(g, numRemainingSteps, penalty);
        nextTree.removeZeroLengthSplits();
        
        double resLogDens=0.0;
        if (!nextTree.fullyResolved()) {
            try {
                double nS = nextTree.getNumInternalSplits();//number of splits in the next tree
                resLogDens = TreeResolver.resolveTree(nextTree, unifDist);
                resLogDens+=(treeA.getNumTaxa()-3-nS)*Math.log(2);
            } catch (AlgorithmException ex) {
                System.out.println(ex.getMessage()+"\nError resolving tree during bridge construction.");
            }
        }
        
        double[] weightAndSD = calcMixtureWeight(nextTree, sigma);

        double centralDens=0.0, mvnDens=0.0;
        
        //geodesic between next tree and realisation -- main part
        Geodesic h=null;
        try {
            h = new Geodesic(nextTree, realization);
        } catch (AlgorithmError ex) {
            System.out.println("Error making geodesic when computing mixture density.");
        }
        //geodesic between y_{j-1} and realisation -- 'central' part
        Geodesic mh=null;
        try {
            mh = new Geodesic(treeA, realization);
        } catch (AlgorithmError ex) {
            System.out.println("Error making geodesic when computing mixture density.");
        }
        
        /* Compute mixture log density */
            double mlen = mh.getInternalLength();
            centralDens = -(treeA.getNumTaxa()-3)*Math.log(sdRW)-0.5*mlen*mlen/(sdRW*sdRW); //  multivariate normal density
        
            /* Take account of boundaries */
            int mnb = countBoundaries(treeA, realization);
            centralDens += -0.6931472*mnb;
        
        
        if (h.isSimple()) {

            /* Compute mixture log density */
            double len = h.getInternalLength();
            mvnDens = -(nextTree.getNumTaxa()-3)*Math.log(sigma)-0.5*len*len/(sigma*sigma); //  multivariate normal density
        
            /* Take account of boundaries */
            int nb = countBoundaries(nextTree, realization);
            mvnDens += -0.6931472*nb;
            
            /* Take account of resolving */
            mvnDens += resLogDens;
            
            /* Take care of numerical stability */
            
            double logDensity;
            double u1 = Math.log(weightAndSD[0]);
            double u2 = Math.log(1.0-weightAndSD[0]);
            if ((u1+centralDens)>(u2+mvnDens)) {
                logDensity = u1+centralDens;
                logDensity += Math.log(1+Math.exp(u2-u1+mvnDens-centralDens));
            }
            else {
                logDensity = u2+mvnDens;
                logDensity += Math.log(1+Math.exp(u1-u2-mvnDens+centralDens));
            }
            //return logDensity;
            if(logDensity==Double.NEGATIVE_INFINITY) System.out.println("something "+u1+" "+u2+" "+mvnDens+" "+centralDens);
            
            return logDensity;
            
            
        }
        else {
            return Math.log(weightAndSD[0])+centralDens;
        }
    }
    
    
 
    
    /** Get the [0,1] coord of the next point to perturb, based on penalizing geodesic kinks.  */
    public static TreeAsSplits getStepLocation(GeodesicForBridging g, int numSteps, int[] penalty) {
        
        double[] lambda = g.getCriticalLambda();
        int[][] sizes = g.getPartitionSizes();
        
        if (lambda.length==0) {
            penalty[0] = 0;
            return g.getTree(1.0/((double)numSteps));
        }
        
        int totalPenalty = 0;
        for (int i=0; i<lambda.length; i++) {
            totalPenalty += calcPenalty(sizes[0][i], sizes[1][i]);
        }
        
        boolean critical = false;
        double location = 0.0;
        int lambdaInd=-1;
        if (numSteps<=totalPenalty+2) {
            /* The location corresponds to first critical lambda value
            changing the below - changed it so you still move along the goedesic but don't automatically skip to 
            the next boundary point
            */
             //location = lambda[0];
            location=1.0/(double)(numSteps);
            // Find any critical lamba values
            for (int i=0; i<lambda.length; i++) {
                if ((lambda[i]<=location)&&(calcPenalty(sizes[0][i], sizes[1][i])>0)) {
                    location = lambda[i];
                    critical = true;
                    penalty[0] = calcPenalty(sizes[0][i], sizes[1][i]);
                    lambdaInd = i;
                    break;
                }
                if (lambda[i]>location) {
                    break;
                }
            }  
        }
        else {
            /* The location might be a stepping point */
            location = 1.0/((double)(numSteps-totalPenalty));
            // Find any critical lamba values
            for (int i=0; i<lambda.length; i++) {
                if ((lambda[i]<=location)&&(calcPenalty(sizes[0][i], sizes[1][i])>0)) {
                    location = lambda[i];
                    critical = true;
                    penalty[0] = calcPenalty(sizes[0][i], sizes[1][i]);
                    lambdaInd = i;
                    break;
                }
                if (lambda[i]>location) {
                    break;
                }
            }
        }
        
        /* Get the tree */
        if (critical) {
            /* location corresponds to a kink */
            return g.getLimitTree(lambdaInd, true);
        }
        else {
            penalty[0] = 0;
            return g.getTree(location);
        }
    }
    
    /** Compute the penalty (measured in number of steps) for a particular "kink" in a geodesic */
    protected static int calcPenalty(int k1, int k2) {
        int k = (k1>k2) ? k1 : k2;
        int r = (k>1) ? k : 0;
        return 1*r;
        //return 5;
    }
    
    
    /* Sampling methods ----------------------------------------------------- */
    
    /** Class for sorting edges -- used by MVN sampler */
    static public class IntDoublePair implements Comparable {
        public int theInt;
        public double theDouble;
        
        public IntDoublePair(int i, double d) {
            theInt = i;
            theDouble = d;
        }

        public int compareTo(Object t) {
            IntDoublePair arg = (IntDoublePair) t;
            if (theDouble<arg.theDouble) return -1;
            if (theDouble>arg.theDouble) return 1;
            return 0;
        }        
    }
    
    
    /** Perturb tree via "direct" MVN Step.
        Argument sigma is SD of isotropic multivariate normal distro used to sample direction and distance.
    */
    static public StepResult sampleDirectMVN(TreeAsSplits theTree, double sigma, NormalDistribution normDist, DoubleUniform unifDist) throws AlgorithmException {

        normDist.setParams(0.0, sigma);
        int nprime = theTree.getNumTaxa()-3;
        double logDens = 0.0;
        
        /* Make an array of splits */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(theTree.getNonTrivialSplits()); 
        /* SORT array list if necessary */
//        Collections.sort(splits); // Sort to ensure repeatability with same random seed       
        
        /* Make a uniform direction */
        double[] u = new double[nprime];
        double lenSqu = 0.0;
        for (int i=0; i<nprime; i++) {
            u[i] = normDist.sample();
            lenSqu += u[i]*u[i];
        }
        /* Normalize and compute edge lengths */
        double[] l = new double[nprime];
        double len = Math.sqrt(lenSqu);
        logDens = -nprime*Math.log(sigma)-0.5*lenSqu/(sigma*sigma); //  multivariate normal density

        ArrayList<ForwardStepBridge.IntDoublePair> boundaries = new ArrayList();
        for (int i=0; i<nprime; i++) {
            Split p = splits.get(i);
            double y = theTree.getSplitLength(p);
            l[i] = y+u[i];
            try {
                // Set edge lengths
                theTree.setSplitLength(p, Math.abs(l[i]));
            } catch (AlgorithmError ex) {
                System.out.println("Error in sample direct MVN in Bridge algorithm: missing split.");
            }
            if ((l[i]<0.0)&&(u[i]<0.0)) {
                boundaries.add(new ForwardStepBridge.IntDoublePair(i,-(y*len)/u[i])); // -y*len/u[i] is the "time" when you hit the boundary. Order by time!
                // Modify log density
                if (y>0.0) { 
                    /* y will be zero if the tree was unresolved and splits have been added
                    If y>0 then we have crossed a genuine boundary: add -log(2) */
                    logDens -= 0.6931472;
                }
            } 
        }
                
        if (boundaries.size()!=0) {
            /* We crossed at least one boundary */
            Collections.sort(boundaries);
            for (int i=0; i<boundaries.size(); i++) {
                /* Do the nni's in turn  */
                int ind = boundaries.get(i).theInt;
                Split p = splits.get(ind);
                // Sanity check: does theTree contain p?
                if (!theTree.contains(p)) {
                    throw new AlgorithmError("Missing p in direct MVN sampler.");
                }
                Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(theTree, p);
                Split q = nniSplit[unifDist.nextIntFromTo(0, 1)];
                theTree.remove(p);
                theTree.add(q, Math.abs(l[ind]));    // UNCOMMENT AND REPLACE LINE BELOW AFTER DEBUG
                //theTree.addWithCheck(q, Math.abs(l[ind]));
            }
        }
                
        StepResult info = new StepResult(theTree, logDens);
     
        return info;
    }

    
    /* Methods for bridge construction -------------------------------------- */
    
    
    /** Independence proposal */
    public double independenceProposal(double t0) throws AlgorithmException {
        
        totalLogLike = Double.NaN;
        theTrees[0] = treeA;
        theTrees[numSteps] = treeB;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        double totalLogPropDensity = 0.0;
        
        for (int k=1; k<numSteps; k++) {
            
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            StepResult res = stepForwardViaMixture(theTrees[k-1], theTrees[numSteps], sd, RWsd, numSteps+1-k, normDist, unifDist);
            theTrees[k] = res.tree;
            Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
            if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing bridge in step "+Integer.toString(k));
            totalLogPropDensity += res.logDensity;
            /* Compute and store log like */
            double[] step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0);
            logLike[k] = step[0];
            logStepPropDensity[k] = step[1];
            totalLogPropDensity += step[1]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
            stepLengths[k]=g.getInternalLength();
            boundaryCrosses[k]=g.numOrthants()-1;
        }
        
        /* Final step: build link between last two trees */
        Geodesic g = new Geodesic(theTrees[numSteps-1], theTrees[numSteps]);
        if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing bridge in final step");
        /* Compute final log like and totalLogPropDensity */
        double[] step = calcStepLogLike(theTrees[numSteps-1], theTrees[numSteps], g, t0);
        logLike[numSteps] = step[0];
        logStepPropDensity[numSteps] = step[1];
        totalLogPropDensity += step[1];
        stepLengths[numSteps]=g.getInternalLength();
        boundaryCrosses[numSteps]=g.numOrthants()-1;
        
        return totalLogPropDensity;
    }
    
    
    /** Compute independence proposal log density */
    public double computeIndependenceLogDensity(double t0) throws AlgorithmException {
        
        double totalLogPropDensity = 0.0;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        for (int k=1; k<numSteps; k++) {         
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            totalLogPropDensity += computeMixtureStepDensity(theTrees[k-1], theTrees[numSteps], theTrees[k], sd, RWsd, numSteps+1-k, unifDist);
            totalLogPropDensity += logStepPropDensity[k]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
        }
        
        totalLogPropDensity += logStepPropDensity[numSteps];
        
        return totalLogPropDensity;
    }
    
    

    /** Sample partial bridge.
        * a is between 1 and numSteps-1,
        * b is between a+1 and numSteps, so b is where you are heading
        * re-sample trees a up to b-1*/
    public double partialBridgeProposal(double t0, int a, int b) throws AlgorithmException {
        
        totalLogLike = Double.NaN;
        double var = t0/numSteps;    
        double totalLogPropDensity = 0.0;
        double RWsd = Math.sqrt(var);
        
        for (int k=a; k<b; k++) {
            double sd = Math.sqrt(var*((double)(b-k))/((double)(b-k+1)));
            StepResult res = stepForwardViaMixture(theTrees[k-1], theTrees[b], sd, RWsd, b+1-k, normDist, unifDist);
                        
            theTrees[k] = res.tree;
            Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
            if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing partial bridge.");
            totalLogPropDensity += res.logDensity;
            
            /* Compute and store log like */
            double[] step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0);
            logLike[k] = step[0];
            logStepPropDensity[k] = step[1];
            totalLogPropDensity += step[1]; // This term corresponds to any proposed quantities for spanning a simple geodesic.     
            stepLengths[k]=g.getInternalLength();
            boundaryCrosses[k]=g.numOrthants()-1;
        }
        
        /* Final step: build link between last two trees */
        Geodesic g = new Geodesic(theTrees[b-1], theTrees[b]);
        if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing partial bridge.");
        /* Compute final log like and totalLogPropDensity */
        double[] step = calcStepLogLike(theTrees[b-1], theTrees[b], g, t0);
        logLike[b] = step[0];
        logStepPropDensity[b] = step[1];
        totalLogPropDensity += step[1];
        stepLengths[b]=g.getInternalLength();
        boundaryCrosses[b]=g.numOrthants()-1;
       
        return totalLogPropDensity;
        
    }
    
    /** Compute independence proposal log density */
    public double computePartialBridgeLogDensity(double t0, int a, int b) throws AlgorithmException {
        
        double totalLogPropDensity = 0.0;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        for (int k=a; k<b; k++) {         
            double sd = Math.sqrt(var*((double)(b-k))/((double)(b-k+1)));
            totalLogPropDensity += computeMixtureStepDensity(theTrees[k-1], theTrees[b], theTrees[k], sd, RWsd, b+1-k, unifDist);
            totalLogPropDensity += logStepPropDensity[k]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
        }
        
        totalLogPropDensity += logStepPropDensity[b];
        
        return totalLogPropDensity;
    }
    
    
    /** Recompute log like and re-propose steps for new t0 */
    
    public void updateForNewt0(double t0) {
        double sigma = Math.sqrt(t0/numSteps);
        for (int k=1; k<=numSteps; k++) {
            try {
                Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
                double[] step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0);
                logLike[k] = step[0];
                logStepPropDensity[k] = step[1];
               
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error updating a bridge for new t0 value. This should not happen: the step should be simple. \nError report --"+anErr.getMessage());
            }
             
        }
        totalLogLike = Double.NaN;
    }
    
    
    
    /** Recompute log like and re-propose steps for new t0 */
    
    public void updateForNewt0new(double t0) throws AlgorithmError, AlgorithmException {
        double sigma = Math.sqrt(t0/numSteps);
        double NPrime=(treeA.getNumTaxa()-3);
        for (int k=1; k<=numSteps; k++) {
            logLike[k] = boundaryCrosses[k]*Math.log(0.5)-NPrime/2*Math.log(2*Math.PI)-NPrime*Math.log(sigma)-0.5*stepLengths[k]*stepLengths[k]/(sigma*sigma);
             Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
             double step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0)[0];
        }
        
        totalLogLike = Double.NaN;
    }
    
   
    
    /** Propose a new bridge given new x0, by shifting trees 1 .. bridgeLength */
    
    public double updateForNewx0(TreeAsSplits x0, double t0, int bridgeLength) throws AlgorithmException {
        
        treeA = x0.efficientClone();
        theTrees[0] = treeA;
        double logPropDensity = partialBridgeProposal(t0, 1, bridgeLength+1);
        totalLogLike = Double.NaN;
        return logPropDensity;
    }
    
    
    //Different class for when simulating directly from the independence proposal in marginal likelihood
    // estimation as fixed mixture weights are used here instead of the variable ones
     public double MargLikeIndependenceProposal(double t0) throws AlgorithmException {
        
        totalLogLike = Double.NaN;
        theTrees[0] = treeA;
        theTrees[numSteps] = treeB;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        double totalLogPropDensity = 0.0;
        
        for (int k=1; k<numSteps; k++) {
            
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            StepResult res = stepForwardViaFixedMixture(theTrees[k-1], theTrees[numSteps], sd, RWsd, numSteps+1-k, normDist, unifDist);
            theTrees[k] = res.tree;
            Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
            if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing bridge in step "+Integer.toString(k));
            totalLogPropDensity += res.logDensity;
            //double logPropDensityCheck =computeFixedMixtureStepDensity(theTrees[k-1], theTrees[numSteps], theTrees[k], sd, numSteps+1-k, unifDist);
            //System.out.println(logPropDensityCheck+" "+res.logDensity);
            /* Compute and store log like */
            double[] step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0);
            logLike[k] = step[0];
            logStepPropDensity[k] = step[1];
            totalLogPropDensity += step[1]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
            stepLengths[k]=g.getInternalLength();
            boundaryCrosses[k]=g.numOrthants()-1;
        }
        
        /* Final step: build link between last two trees */
        Geodesic g = new Geodesic(theTrees[numSteps-1], theTrees[numSteps]);
        if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing bridge in final step");
        /* Compute final log like and totalLogPropDensity */
        double[] step = calcStepLogLike(theTrees[numSteps-1], theTrees[numSteps], g, t0);
        logLike[numSteps] = step[0];
        logStepPropDensity[numSteps] = step[1];
        totalLogPropDensity += step[1];
        stepLengths[numSteps]=g.getInternalLength();
        boundaryCrosses[numSteps]=g.numOrthants()-1;
        
        return totalLogPropDensity;
    }
     
      public double MargLikeIndependenceProposalNoStop(double t0) throws AlgorithmException {
        
        totalLogLike = Double.NaN;
        theTrees[0] = treeA;
        theTrees[numSteps] = treeB;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        double totalLogPropDensity = 0.0;
        
        for (int k=1; k<numSteps; k++) {
            
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            StepResult res = stepForwardViaFixedMixture(theTrees[k-1], theTrees[numSteps], sd, RWsd, numSteps+1-k, normDist, unifDist);
            theTrees[k] = res.tree;
            Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
            totalLogPropDensity += res.logDensity;
            //double logPropDensityCheck =computeFixedMixtureStepDensity(theTrees[k-1], theTrees[numSteps], theTrees[k], sd, numSteps+1-k, unifDist);
            //System.out.println(logPropDensityCheck+" "+res.logDensity);
            /* Compute and store log like */
            double[] step = calcStepLogLike(theTrees[k-1], theTrees[k], g, t0);
            logLike[k] = step[0];
            logStepPropDensity[k] = step[1];
            totalLogPropDensity += step[1]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
            stepLengths[k]=g.getInternalLength();
            boundaryCrosses[k]=g.numOrthants()-1;
        }
        
        /* Final step: build link between last two trees */
        Geodesic g = new Geodesic(theTrees[numSteps-1], theTrees[numSteps]);
        /* Compute final log like and totalLogPropDensity */
        double[] step = calcStepLogLike(theTrees[numSteps-1], theTrees[numSteps], g, t0);
        logLike[numSteps] = step[0];
        logStepPropDensity[numSteps] = step[1];
        totalLogPropDensity += step[1];
        stepLengths[numSteps]=g.getInternalLength();
        boundaryCrosses[numSteps]=g.numOrthants()-1;
        
        return totalLogPropDensity;
    }
    
    
    
    //Different class for when calculating the independence proposal density in marginal likelihood
    // estimation as fixed mixture weights are used here instead of the variable ones
    public double computeMargLikeIndependenceLogDensity(double t0) throws AlgorithmException {
        
        double totalLogPropDensity = 0.0;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        for (int k=1; k<numSteps; k++) {         
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            totalLogPropDensity += computeFixedMixtureStepDensity(theTrees[k-1], theTrees[numSteps], theTrees[k], sd, RWsd, numSteps+1-k, unifDist);
            totalLogPropDensity += logStepPropDensity[k]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
        }
        
        totalLogPropDensity += logStepPropDensity[numSteps];
        
        return totalLogPropDensity;
    }
    
    public double computeMargLikeIndependenceLogDensityNoStop(double t0) throws AlgorithmException {
        
        double totalLogPropDensity = 0.0;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        for (int k=1; k<numSteps; k++) {         
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            totalLogPropDensity += computeFixedMixtureStepDensity(theTrees[k-1], theTrees[numSteps], theTrees[k], sd, RWsd, numSteps+1-k, unifDist);
            totalLogPropDensity += logStepPropDensity[k]; // This term corresponds to any proposed quantities for spanning a simple geodesic. 
        }
        
        totalLogPropDensity += logStepPropDensity[numSteps];
        
        return totalLogPropDensity;
    }
    
    
    
    // more efficient version of MargLikeIndependenceProposal for computing prop simple in the marginal likelihood methods
    // here we only need to know if the bridge has simple steps, not its proposal density of log likelihood
    public double MargLikeIndependenceProposalPropSimple(double t0) throws AlgorithmException {
        
        totalLogLike = Double.NaN;
        theTrees[0] = treeA;
        theTrees[numSteps] = treeB;
        double var = t0/numSteps;
        double RWsd = Math.sqrt(var);
        
        double totalLogPropDensity = 0.0;
        
        for (int k=1; k<numSteps; k++) {
            
            double sd = Math.sqrt(var*((double)(numSteps-k))/((double)(numSteps-k+1)));
            StepResult res = stepForwardViaFixedMixturePropSimple(theTrees[k-1], theTrees[numSteps], sd, RWsd, numSteps+1-k, normDist, unifDist);
            theTrees[k] = res.tree;
            Geodesic g = new Geodesic(theTrees[k-1], theTrees[k]);
            if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing bridge in step "+Integer.toString(k));
        }   
        
        /* Final step: build link between last two trees */
        Geodesic g = new Geodesic(theTrees[numSteps-1], theTrees[numSteps]);
        if (!g.isSimple()) throw new AlgorithmException("Non-simple link constructing bridge in final step");
        /* Compute final log like and totalLogPropDensity */
        
        return totalLogPropDensity;
    }
    
    
    
    /* Test area ------------------------------------------------------------ */

    public static void main(String[] args) {
        
        
//        ChiSquare c = new ChiSquare(20, Random.getEngine());
//        for (int i=0; i<=400; i++) {
//           System.out.println(incompleteGamma(0.5*df,0.5*x);
//            System.out.println(c.cdf(0.1*i));
//        }

        // Data set
        String dataFilename = "/home/ntmwn/research/diffusion/bridge/analysis/yeast.nex"; 
       
        // Initial tree x0 and t0
//        String initialTreeFilename = "/home/ntmwn/research/diffusion/bridge/sims/randomtree_48_gamma_2_20.txt"; 
        String initialTreeFilename = "/home/ntmwn/research/diffusion/bridge/analysis/mean_tree.txt"; 
        double sd = 0.5;
        final double t0 = sd*sd;
        final int numSteps = 50;
        
        Random.setEngine(78129);

        /* End inputs ------------------------------------------------------- */

        
        /* Read in the data */
        File inFile = new File(dataFilename);
        TreeAsSplitsDataSet theData = null;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }

        /* Read in an initial tree */
        Tree initialTree = null;
        try {
            initialTree = new Tree(new File(initialTreeFilename));
            initialTree.removeDegreeTwoVertices();
        }
        catch (AlgorithmException anError) {
            System.out.println("Bad initial tree. "+anError.getMessage());
            System.exit(1);
        }
        TreeAsSplits x0 = new TreeAsSplits(initialTree);
        try {
            StarTree.getInstance().setTree(x0);
        } catch (AlgorithmError ex) {
            System.out.println("Error making star tree.");
            System.exit(1);
        }
        
        DoubleUniform unifDist = new DoubleUniform(Random.getEngine());
        
        for (int k=0; k<theData.numTrees; k++) {
            boolean building = true;
            int count = 0;

            while (building) {

                count++;
                BridgeWithApproxMVNLike fsb = new BridgeWithApproxMVNLike(x0, theData.getTree(k), numSteps);
                try {
                    fsb.independenceProposal(t0);
                    building = false;
                } catch (AlgorithmException ex) {
                    //System.out.println(ex.toString());
                }    
            }
            System.out.println("Tree "+k+" "+count);
        }
    }
}
