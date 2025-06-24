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

package MarginalLikelihoods;

/**
 *  
 */

import bridge.*;
import MCMC.Chain;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import MCMC.KernelWrapper;
import MCMC.MetropolisHastingsKernelWrapper;
import MCMC.PositiveParameter;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.TreeState;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import simulation.CategoricalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/*
Simulate the independence proposals and simulate from the bridge distribuions using MCMC in order to calculate
the Chib and tunnel estimators of marginal likelihood in the fixed t0 case (potentially with multiple data points)
*/
public class ChibJeliMultPointsYeastFM {
    
        public static void displaySyntaxMessage() {
        System.out.println("Syntax: InferBrownianParams data_file initial_tree_file sqrt_t0 numSteps seed\n"
                + "    [-n num_its] [-b burn_its] [-t thin] [-o outfile] \n"
                + "    [-x0sph shape scale] [-x0std shape scale] [-t0gam shape scale]\n"
                + "    [-pbg geom_fac] [-pbl len] [-pbf start end] \n"
                + "    [-numProps num]\n");
                
    }
        
        /*
        Simulate from the bridge distributions using MCMC
        */
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
        Tree initialTree = null;
        try {
            initialTree = new Tree(new File(initialTreeFileName));
            initialTree.removeDegreeTwoVertices();
        }
        catch (AlgorithmException anError) {
            System.out.println("Bad initial tree. "+anError.getMessage());
            System.exit(1);
        }
        TreeAsSplits x0 = new TreeAsSplits(initialTree);
        try {
            double sd = Double.parseDouble(args[2]);
            t0 = sd*sd;
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
        
        /* Build initial bridges */
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForBridgingTunnel theBrownianState = null; 
        ForwardStepBridge template = new BridgeWithApproxMVNLike(x0, null, numSteps);
        try {
            for(int i=0; i<2; i++) {
                BrownianStateForBridgingTunnel state = new BrownianStateForBridgingTunnel(x0, t0, theData.theTrees, (i==0), template);
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
            /* directly sample bridges from the independence proposal for the denominator
       store the likelihoods and independence proposal densities
            */
       else if(a.equals("-numProps")) {
                try {
                    int numProps = Integer.parseInt(args[i+1]);
                    String propsFileName = args[i+2];
                    File propsFile = new File(propsFileName);
                    makeProposalData(propsFile, numProps, x0, theData, numSteps,t0);                   
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with fixed end points for bridge proposal on command line.");
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
        }
        
       
        /* Build the sweep proposal */
        KernelWrapper myKernelWrapper = new SweepKernelWrapper("Overall proposal", kernelWrapperList);
       
 
        /* Run the chain targetting the bridge distributions*/
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
    
    private static void makeProposalData(File propFile, int numProp,TreeAsSplits startTree, TreeAsSplitsDataSet Data, int m,double t0) throws AlgorithmException{
        //get the number of data points
        int n = Data.numTrees;
        ArrayList<TreeAsSplits> theTrees = Data.theTrees;
        ArrayList<ForwardStepBridge> theBridges = new ArrayList();
        double[][] outputs = new double[numProp][2*n];
        int nPrime =startTree.getNumTaxa()-3;
        for(int i=0; i<n;i++){
            ForwardStepBridge theBridge = new BridgeWithApproxMVNLike(startTree, theTrees.get(i), m);//initialise bridges with usual likelihood
            theBridges.add(theBridge);
            }
        double logLike = 0;
        double propDens = 0;
        for(int j=0; j<numProp;j++){
            for(int i=0; i<n;i++){ 
            logLike = 0;
            propDens = 0;
            try {
            theBridges.get(i).MargLikeIndependenceProposal(t0);//simulate the proposals
            propDens = theBridges.get(i).computeMargLikeIndependenceLogDensity(t0)-Math.log(2*Math.PI)*nPrime*(m-1)/2;
            logLike = theBridges.get(i).getTotalLogLike();
            }
            catch(treebase.AlgorithmException anErr){  
                logLike = Double.NEGATIVE_INFINITY;
            }
           if(propDens==Double.NEGATIVE_INFINITY){
               System.out.println("Should not produce a bridge with zero proposal density");
           }
            outputs[j][i*2+1] = logLike;
            outputs[j][i*2] =propDens;
            }
            
        }
        
        //print the outputs to a file
        PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(propFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+propFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
            
}

//print a header
for(int i=0;i<n;i++){
   if(i==0) out.print("PropDens"+i+" LogLike"+i); 
   else out.print(" PropDens"+i+" LogLike"+i); 
}

out.println("");

for(int j=0 ; j<numProp ; j++){
    for(int i=0; i<n; i++){
      
    if(i==0) out.print(outputs[j][2*i]+" "+outputs[j][2*i+1]);
    else  out.print(" "+outputs[j][2*i]+" "+outputs[j][2*i+1]);
    
    }
    out.println("");
    } 
 out.close();  
        
        
        
    }

    public static void main(String[] args) throws AlgorithmException {
        
        //arguments the MCMC to get the posterior sample and for directly simulating independence proposals
        /*
        args = new String[18];
        args[0] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/FixedDispMultPoints/4Taxa/DataPoints.txt";
        args[1] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/FixedDispMultPoints/4Taxa/x0.txt";
        args[2] = "0.5"; // Initial squ root t_0
        args[3] = "20"; // Num steps
        args[4] = "1260"; // Seed
        args[5] = "-n";
        args[6] = "15000"; // Num interations - before thin
        args[7] = "-t";
        args[8] = "20"; // thin
        args[9] = "-b";
        args[10] = "1000";// burn-in
        args[11] = "-o";
        args[12] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/MultiplePointTest/Post.txt";
        args[13] = "-pbg";
        args[14] = "0.01";//partial bridge proposals
        args[15] = "-numProps";
        args[16] = "10000"; //number of independence proposals to simulate directly for each data point
        args[17] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/MultiplePointTest/Props.txt"; //file to put the proposal data in
        
        for(int i=0; i<5;i++){
          args[4] = Integer.toString(174*i^2);  
          args[12] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/FixedDispMultPoints/4Taxa/ChibPostOutTest"+i+".txt";
          args[17] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/FixedDispMultPoints/4Taxa/ChibPropDataTest"+i+".txt";
          parseCommandLine(args);
        }
        
*/
            args = new String[18];
        args[0] = "/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints.txt";
        args[1] = "/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_OP.txt";
        args[2] = "0.128"; // Initial squ root t_0
        args[3] = "50"; // Num steps
        args[4] = "1260"; // Seed
        args[5] = "-n";
        args[6] = "1000000"; // Num interations - before thin
        args[7] = "-t";
        args[8] = "100"; // thin
        args[9] = "-b";
        args[10] = "100000";// burn-in
        args[11] = "-o";
        args[12] = "/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_FM_OP_ChibJeliPostv2.txt";
        args[13] = "-pbg";
        args[14] = "0.05";//partial bridge proposals
        args[15] = "-numProps";
        args[16] = "50000"; //number of independence proposals to simulate directly for each data point
        args[17] = "/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_FM_OP_ChibJeliPropsv2.txt"; //file to put the proposal data in
        
          parseCommandLine(args);
       
    }

}