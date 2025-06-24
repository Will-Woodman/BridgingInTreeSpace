/*
    InferBrownianParamsMCMC
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

/**
 *  
 */

import MCMC.Chain;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import MCMC.KernelWrapper;
import MCMC.MetropolisHastingsKernelWrapper;
import MCMC.PositiveParameter;
import static MCMC.PosteriorAnalysis.countTopologies;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.TreeState;
import geodesics.Geodesic;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import simulation.CategoricalDistribution;
import simulation.Random;
import static topologies.BacakAlgorithm.unweightedFM;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;


public class InferBrownianParamsMCMC {
    
    /*Version of InferBrownianParamsMCMC that implements the methods of choosing initial values from the paper
    -Initial t0 is the Frechet variance
    -Initial x0 is the closest data point to the Frechet mean
    */
    
        public static void displaySyntaxMessage() {
        System.out.println("Syntax: InferBrownianParams data_file initial_tree_file initial_sqrt_t0 numSteps seed\n"
                + "    [-n num_its] [-b burn_its] [-t thin] [-o outfile] \n"
                + "    [-x0sph shape scale] [-x0std shape scale] [-t0gam shape scale]\n"
                + "    [-pbg geom_fac] [-pbl len] [-pbf start end] \n"
                + "    [-ptr sigma] \n"
                + "    [-pxf sigma length] [-pxg sigma geom]\n");
                
    }

    public static void parseCommandLine(String[] args) throws AlgorithmException {

        int seed = 1;
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        int nits = 100, burnits = 0, thin = 1, numSteps = 50;
        double t0 = 0.01;
        File outFile = null;
        
        /* Parse arguments */
        if ((args.length < 4)||(args.length > 32)) {
            displaySyntaxMessage();
            System.exit(1);
        }

        int n = args.length;
        String dataFileName = args[0];
        String initialTreeFileName = args[1];
        
        /* Read in data */
        File inFile = new File(dataFileName);
        TreeAsSplitsDataSet theData = null;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        
        /* Read in an initial tree */
        TreeAsSplits FrechMean =null;
        Tree initialTree = null;
        if("closest".equals(initialTreeFileName)){
            FrechMean = unweightedFM(theData.theTrees,1000);
            initialTree = SourceClosestToFM(theData.theTrees,FrechMean).getTree();
            System.out.println("Initial x0 = "+initialTree.toString());
        }
        else if("random".equals(initialTreeFileName)){
           
            initialTree =RandomSourceFromModalTop(theData.theTrees).getTree();
            System.out.println("Initial x0 = "+initialTree.toString());
        }
        else{
        try {
            initialTree = new Tree(new File(initialTreeFileName));
            initialTree.removeDegreeTwoVertices();
            System.out.println(initialTree.numVertices());
        }
        catch (AlgorithmException anError) {
            System.out.println("Bad initial tree. "+anError.getMessage());
            System.exit(1);
        }
        }
        TreeAsSplits x0 = new TreeAsSplits(initialTree);
        try {
            double sd = Double.parseDouble(args[2]);
            if(sd==0){//then we want to use the Frechet variance
                //if(FrechMean==null) FrechMean = unweightedFM(theData.theTrees,1000);
                //output the Frechet variance
                t0 = getFrechetVariance(x0, theData.theTrees);
                
            }
            else t0 = sd*sd;
            System.out.println("Initial t0 = "+String.format("%7.7f", t0));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read initial t0 from command line.");
            System.exit(1);
        }
        try {
            numSteps = Integer.parseInt(args[3]);
            System.out.print("Number of random walk steps in bridge = "+String.format("%d%n", numSteps));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of RW steps from command line.");
            System.exit(1);
        }
        
        try {
            seed = Integer.parseInt(args[4]);
            System.out.print("Seed = "+String.format("%d%n", seed));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read seed from command line.");
            System.exit(1);
        }
        Random.setEngine(seed);
        TreeResolver.resolveTree(x0, new DoubleUniform(Random.getEngine()));
        
        /* Build the prior */
        
        // Prior on sigma squ
        double rateForT0Prior = 18.44*((double)(initialTree.numTaxa()-3))/((double)(initialTree.numTaxa()));
        DensityCalculator t0Prior = new PositiveParameter.ExponentialPrior(rateForT0Prior);
        System.out.println("Exponential prior on t0 with rate = "+String.format("%7.7f", rateForT0Prior));
        // Prior on x0
        DensityCalculator treePrior = new TreeState.SphericalTreePrior(initialTree.numTaxa());
        //DensityCalculator treePrior = new TreeState.NormalTreePrior();
        //DensityCalculator treePrior = new TreeState.ChiSquTreePrior();
        System.out.println("Normal prior on x0: "+treePrior.toString());
        // Combined prior
        DensityCalculator thePrior = new BrownianStateForBridging.SimpleBrownianStatePrior(treePrior, t0Prior); 
        /* Build initial bridges */
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForBridging theBrownianState = null; 
        ForwardStepBridge template = new BridgeWithApproxMVNLike(x0, null, numSteps);
        try {
            for(int i=0; i<2; i++) {
                BrownianStateForBridging state = new BrownianStateForBridging(x0, t0, theData.theTrees, (i==0), template);
                if (i==0) theBrownianState = state;
                java.util.HashSet<State> set = new java.util.HashSet<State>();
                set.add(state);
                initialState[i] = new GlobalState(set);
                initialState[i].setLogLikelihood(state.getLogLike());
             }
        }
        catch (AlgorithmException ex) {
            System.out.println("Unable to initialize bridges: quitting.");
            System.exit(1);
        }
        System.out.println("Made the initial states...");
       
       
        // Loop thru arguments
        int i=5;
        while (i<n) {
            String a = args[i];
            
            if(a.equals("-n")) {
                try {
                    int temp = Integer.parseInt(args[i+1]);
                    nits =  (temp>1) ? temp : 1;
                    System.out.println("Number of iterations = "+nits);
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Unable to read number of iterations from command line.");
                    System.exit(1);
                }
                i=i+2;
            }
            
            else if(a.equals("-b")) {
                try {
                    int temp = Integer.parseInt(args[i+1]);
                    burnits =  (temp>1) ? temp : 1;
                    System.out.println("Burn iterations = "+burnits);
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Unable to read number of burn iterations from command line.");
                    System.exit(1);
                }
                i=i+2;
            }
            
            else if(a.equals("-t")) {
                try {
                    int temp = Integer.parseInt(args[i+1]);
                    thin =  (temp>1) ? temp : 1;
                    System.out.println("Thin iterations = "+thin);
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Unable to read number of thin iterations from command line.");
                    System.exit(1);
                }
                i=i+2;
            }
            
            else if(a.equals("-o")) {
                outFile = new File(args[i+1]);
                i=i+2;
            }
            
            else if(a.equals("-pbg")) {
                try {
                    double geomFac = Double.parseDouble(args[i+1]);
                    if ((geomFac<=0)||(geomFac>=1)) throw new NumberFormatException();
                    kernelWrapperList.add(theBrownianState.makePartialBridgeSweepProposal(makeTruncatedGeometricDistributionKTrials(geomFac, numSteps)));
                    System.out.println("Bridge proposal with geometric length added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with geometric probability for bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+2;
            }

            else if(a.equals("-ibp")) {
                try {
                    kernelWrapperList.add(theBrownianState.makeIndependenceBridgeSweepProposal());
                    System.out.println("Bridge proposal with fixed length added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with fixed length for bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+2;
            }
            
            else if(a.equals("-pbl")) {
                try {
                    int fixedLen = Integer.parseInt(args[i+1]);
                    if ((fixedLen<=0)||(fixedLen>numSteps)) throw new NumberFormatException();
                    kernelWrapperList.add(theBrownianState.makePartialBridgeSweepProposal(fixedLen));
                    System.out.println("Bridge proposal with fixed length added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with fixed length for bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+2;
            }

            else if(a.equals("-pbf")) {
                try {
                    int fixedStart = Integer.parseInt(args[i+1]);
                    int fixedEnd = Integer.parseInt(args[i+2]);
                    if ((fixedStart<=0)||(fixedStart>=fixedEnd)||(fixedEnd>numSteps)) throw new NumberFormatException();
                    kernelWrapperList.add(theBrownianState.makePartialBridgeSweepProposal(fixedStart, fixedEnd));
                    System.out.println("Bridge proposal with fixed end points added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with fixed end points for bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+3;
            }

            else if(a.equals("-ptr")) {
                try {
                    double sigma = Double.parseDouble(args[i+1]);
                    if (sigma<=0.0) throw new NumberFormatException();
                    Kernel t0RetainProp = theBrownianState.new T0RandomWalkProposal(sigma);
                    KernelWrapper t0PartialWrapper = new MetropolisHastingsKernelWrapper(t0RetainProp, thePrior, new BrownianStateForBridging.BrownianStateLikelihoodCalculator(), theBrownianState.getName());
                    kernelWrapperList.add(t0PartialWrapper);
                    System.out.println("t0 random walk proposal added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for t0 random walk proposal on command line.");
                    System.exit(1);
                }
                i=i+2;
            }

            else if(a.equals("-pxf")) {
                try {
                    double sigma = Double.parseDouble(args[i+1]);
                    if (sigma<=0.0) throw new NumberFormatException();
                    int fixedLen = Integer.parseInt(args[i+2]);
                    if ((fixedLen<0)||(fixedLen>numSteps)) throw new NumberFormatException();
                    Kernel x0BridgeProp = theBrownianState.new X0DirectMVNProposalWithHierarchySubsample(sigma, fixedLen);
                    KernelWrapper x0Wrapper = new MetropolisHastingsKernelWrapper(x0BridgeProp, thePrior, new BrownianStateForBridging.BrownianStateLikelihoodCalculator(), theBrownianState.getName());
                    kernelWrapperList.add(x0Wrapper);
                    System.out.println("x0 proposal with fixed length bridge proposal added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 proposal with fixed length bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+3;
            }

            else if(a.equals("-pxg")) {
                try {
                    double sigma = Double.parseDouble(args[i+1]);
                    if (sigma<=0.0) throw new NumberFormatException();
                    double geomFac = Double.parseDouble(args[i+2]);
                    if ((geomFac<=0)||(geomFac>=1)) throw new NumberFormatException();
                    Kernel x0BridgeProp = theBrownianState.new X0DirectMVNProposalWithHierarchySubsample(sigma, makeTruncatedGeometricDistribution(geomFac, numSteps));
                    KernelWrapper x0Wrapper = new MetropolisHastingsKernelWrapper(x0BridgeProp, thePrior, new BrownianStateForBridging.BrownianStateLikelihoodCalculator(), theBrownianState.getName());
                    kernelWrapperList.add(x0Wrapper);
                    System.out.println("x0 proposal with geometric length bridge proposal added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 proposal with geometric length bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+3;
            }
            
            else if(a.equals("-pxgi")) {
                try {
                    double sigma = Double.parseDouble(args[i+1]);
                    if (sigma<=0.0) throw new NumberFormatException();
                    double geomFac = Double.parseDouble(args[i+2]);
                    if ((geomFac<=0)||(geomFac>=1)) throw new NumberFormatException();
                    //Kernel x0BridgeProp = theBrownianState.new X0CGDProposalWithHierarchySubsample(sigma, makeTruncatedGeometricDistribution(geomFac, numSteps));
                    Kernel x0BridgeProp = theBrownianState.new X0CGDProposalWithHierarchySubsample(sigma, numSteps-1);
                    KernelWrapper x0Wrapper = new MetropolisHastingsKernelWrapper(x0BridgeProp, thePrior, new BrownianStateForBridging.BrownianStateLikelihoodCalculator(), theBrownianState.getName());
                    kernelWrapperList.add(x0Wrapper);
                    System.out.println("x0 proposal with geometric length bridge proposal added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 proposal with geometric length bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+3;
            }
            
            else if(a.equals("-pxgii")) {
                try {
                    double sigma = Double.parseDouble(args[i+1]);
                    if (sigma<=0.0) throw new NumberFormatException();
                    
                    //Kernel x0BridgeProp = theBrownianState.new X0CGDProposalWithHierarchySubsample(sigma, makeTruncatedGeometricDistribution(geomFac, numSteps));
                    Kernel x0BridgeProp = theBrownianState.new X0CGDProposalWithNewBridges(sigma);
                    KernelWrapper x0Wrapper = new MetropolisHastingsKernelWrapper(x0BridgeProp, thePrior, new BrownianStateForBridging.BrownianStateLikelihoodCalculator(), theBrownianState.getName());
                    kernelWrapperList.add(x0Wrapper);
                    System.out.println("x0 proposal with geometric length bridge proposal added.");
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 proposal with geometric length bridge proposal on command line.");
                    System.exit(1);
                }
                i=i+3;
            }
            

            else if(a.equals("-x0sph")) {
                try {
                    double shape = Double.parseDouble(args[i+1]);
                    if (shape<=0.0) throw new NumberFormatException();
                    double scale = Double.parseDouble(args[i+2]);
                    if (scale<0.0) throw new NumberFormatException();
                    if (kernelWrapperList.size()>0) {
                        System.out.println("Error with input: priors must be specified before proposals.");
                        System.exit(1);
                    }
                    treePrior = new TreeState.SphericalTreePrior(shape, scale, false,initialTree.numTaxa()); // False means ignore pendants
                    thePrior = new BrownianStateForBridging.SimpleBrownianStatePrior(treePrior, t0Prior); 
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 spherical prior on command line.");
                    System.exit(1);
                }
                i=i+3;
            }
            
            else if(a.equals("-x0std")) {
                try {
                    double shape = Double.parseDouble(args[i+1]);
                    if (shape<=0.0) throw new NumberFormatException();
                    double scale = Double.parseDouble(args[i+2]);
                    if (scale<0.0) throw new NumberFormatException();
                    if (kernelWrapperList.size()>0) {
                        System.out.println("Error with input: priors must be specified before proposals.");
                        System.exit(1);
                    }
                    treePrior = new TreeState.StandardTreePrior(shape, scale, false); 
                    thePrior = new BrownianStateForBridging.SimpleBrownianStatePrior(treePrior, t0Prior); 
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 standard phylogenetic prior on command line.");
                    System.exit(1);
                }
                i=i+3;
            }
            //option to implement alternative prior on t0:
            else if(a.equals("-t0gam")) {
                try {
                    double shape = Double.parseDouble(args[i+1]);
                    if (shape<=0.0) throw new NumberFormatException();
                    double scale = Double.parseDouble(args[i+2]);
                    if (scale<0.0) throw new NumberFormatException();
                    if (kernelWrapperList.size()>0) {
                        System.out.println("Error with input: priors must be specified before proposals.");
                        System.exit(1);
                    }
                    t0Prior = new PositiveParameter.GammaPrior(shape, scale);
                    thePrior = new BrownianStateForBridging.SimpleBrownianStatePrior(treePrior, t0Prior); 
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with parameters for x0 standard phylogenetic prior on command line.");
                    System.exit(1);
                }
                i=i+3;
            }

            else {
                /* Unrecognized option */
                System.out.println("Bad argument: "+a+"\n");
                displaySyntaxMessage();
                System.exit(1);
            }


        } // Finished parsing
        
        
        /* Build the sweep proposal */
        KernelWrapper myKernelWrapper = new SweepKernelWrapper("Overall proposal", kernelWrapperList);
       
 
        /* Run the chain */
        Chain myChain = new Chain(initialState[0], initialState[1], myKernelWrapper);
        myChain.run(nits, burnits, thin, outFile);
        
    }

    
    /** Make a truncated geom dist, supported on 0 .. max-1 */
    private static CategoricalDistribution makeTruncatedGeometricDistribution(double pr, int max) {
        double[] p = new double[max];
        double s = 0.0;
        double q = 1-pr; 
        for (int i=0; i<max-1; i++) {
            p[i] = Math.pow(q, i);
            s += p[i];
        }
        for (int i=0; i<max-1; i++) {
            p[i] = p[i]/s;
        }
        CategoricalDistribution catDist = new CategoricalDistribution(p);
        return catDist;
    } 
    
      /** Make a truncated geom dist based on number of trials, supported on 1 .. max-1 
        used for the bridge proposal where l=0 should not be possible */
    private static CategoricalDistribution makeTruncatedGeometricDistributionKTrials(double pr, int max) {
        double[] p = new double[max];
        double s = 0.0;
        double q = 1-pr; 
        for (int i=1; i<max-1; i++) {
            p[i] = Math.pow(q, i);
            s += p[i];
        }
        for (int i=0; i<max-1; i++) {
            p[i] = p[i]/s;
        }
        CategoricalDistribution catDist = new CategoricalDistribution(p);
        return catDist;
    } 
    
    private static class TreePair{
      TreeAsSplits treeA;
      TreeAsSplits treeB;
      double squDist;
      
      private TreePair(TreeAsSplits tA, TreeAsSplits tB) throws AlgorithmError{
          treeA=tA;
          treeB=tB;
          Geodesic g = new Geodesic(tA,tB);
          double dist = g.getInternalLength();
          squDist=dist*dist;

             
      }
     
  }
    
    //method for choosing the initial value of x0 to be a random data point from the modal topology in the data set
     private static TreeAsSplits RandomSourceFromModalTop(ArrayList<TreeAsSplits> theTrees) throws AlgorithmError{
     
          LinkedHashMap<String,Integer> myList = countTopologies(theTrees);
          //take account of the possibility of multiple modal topologies
          ArrayList<String> modalTops= new ArrayList();
          Iterator<String> myIterator= myList.keySet().iterator();
          String modalTop = myIterator.next();
          int numInModalTop = myList.get(modalTop);
          int numInNextTop = 0;
          modalTops.add(modalTop);
          boolean foundAll = false;
          while(foundAll==false & myIterator.hasNext()){
              System.out.println(foundAll+" "+myIterator.hasNext());
              String nextTop = myIterator.next();
              System.out.println(nextTop);
              numInNextTop = myList.get(nextTop);
              if(numInNextTop==numInModalTop) modalTops.add(nextTop);
              else foundAll=true;
          }
          
          //choose one of the possible modal topologies at random
          int numModalTops= modalTops.size();
          System.out.println(numModalTops);
          CategoricalDistribution topDist = new CategoricalDistribution(numModalTops);
          int ind1 = topDist.sample();
          modalTop= modalTops.get(ind1);
          
          
          //Get the trees from the chosen modal topology
          ArrayList<TreeAsSplits> TreesInModalTop = new ArrayList();
          for(int i =0; i< theTrees.size();i++){
              if(theTrees.get(i).getTree().toTopologyString().equals(modalTop)) {
                  TreesInModalTop.add(theTrees.get(i));
              }
          }
           
         //choose a tree from the modal topology at random
          CategoricalDistribution catDist = new CategoricalDistribution(numInModalTop);
          int ind = catDist.sample();
          return(TreesInModalTop.get(ind));
     }
     
     
     private static TreeAsSplits RandomSourceFromDataSet(ArrayList<TreeAsSplits> theTrees) throws AlgorithmError{
         int numTrees= theTrees.size();
          CategoricalDistribution catDist = new CategoricalDistribution(numTrees);
          int ind = catDist.sample();
          return(theTrees.get(ind));
     }
     
     //method for choosing the initial value for x0 as the closest data point to the Frechet mean
     private static TreeAsSplits SourceClosestToFM(ArrayList<TreeAsSplits> theTrees, TreeAsSplits mu) throws AlgorithmError{
         int numTrees= theTrees.size();
         double minDist =100000000;
         double len =0;
         int ind =0;
         for(int i=0 ; i< numTrees;i++){
             Geodesic g = new Geodesic(theTrees.get(i),mu);
             len = g.getInternalLength();
             if(len<minDist){
                 minDist=len;
                 ind =i;
             }
         }
         TreeAsSplits outTree = theTrees.get(ind);
         System.out.println(outTree.toString());
         return outTree;
     }
     
     private static Double getFrechetVariance (TreeAsSplits mu, ArrayList<TreeAsSplits> theTrees) throws AlgorithmError{
         int numTrees= theTrees.size();
         double totalSquDist=0;
         for(int i=0;i<numTrees;i++){
            TreePair theTreePair = new TreePair(theTrees.get(i),mu);
         totalSquDist+= theTreePair.squDist;
         System.out.println(theTreePair.squDist);
        }
         System.out.println((double)((mu.getNumTaxa()-3)*numTrees));
         return(totalSquDist/((double)((mu.getNumTaxa()-3)*numTrees)));
         
     }

    public static void main(String[] args) throws AlgorithmException {
        
       for(String arg:args){
           System.out.println(arg);
       }
        parseCommandLine(args);


    }

}