/*
SteppingStoneSampler
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

package MarginalLikelihoods;

/**
 *  
 */

import bridge.*;
import MCMC.Chain;
import MCMC.GlobalState;
import MCMC.KernelWrapper;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import cern.jet.random.tdouble.DoubleUniform;
import java.io.File;
import java.util.ArrayList;
import simulation.CategoricalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;
import bridge.BridgeState;
import java.util.Arrays;

/*
Run the MCMC needed to compute a stepping stone sampling estimate of the marginal likelihood when t0 is fixed.
Use the quasistatic method so burnin is only needed to be run once. First distribution (independence proposals) is
simulated from directly
*/

public class SteppingStoneSampler {
    
        public static void displaySyntaxMessage() {
        System.out.println("Syntax: SteppingStoneFixed data_file initial_tree_file initial_sqrt_t0 numSteps seed\n"
                + "    [-n num_its] [-b burn_its] [-t thin] [-o outfile] \n"
                + "    [-pbg geom_fac] [-pbl len] [-pbf start end] \n");
                
    }

    public static void parseCommandLine(String[] args) throws AlgorithmException {
        Double[] betas=new Double[1001];
        /*
        set the values for beta_k at which to sample from the unnormalised distribution
        here using equally spaced beta_k's.
        */
        for(int i=0;i<101;i++){
    betas[i]=(double)i/100;
}
        int seed = 1;
        
        int nits = 100, burnits = 0, thin = 1, numSteps = 50;
        double t0 = 0.01;
        File outFile = null;
        
        /* Parse arguments */
        if ((args.length < 4)||(args.length > 22)) {
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
        
        /* Read in the source tree */
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
        
         ArrayList<BridgeState> finalBridges = new ArrayList();
        
        //get the size of the data set
        int numDataPoints = theData.numTrees;        
        
        /*
        now loop through and sample for each value of beta_k
        start at j=1 since we simulate from the proposal distribution directly
        leave the creation of proposals and the reading in of number iterations etc
        within the loop so created fresh each time + allows for small tweaks so that we could specifiy different
        parameters for different beta_k etc
        */
        for(int j=1; j<100; j++){
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForBridging theBrownianState = null; 
        BridgeWithStepStoneLike template = new BridgeWithStepStoneLike(x0, null, numSteps,betas[j]);
        
        
        /* If this is the first MCMC run build initial bridges -- otherwise we pass in last point from
        previous run*/
        if(j<2){
        try {
            for(int i=0; i<2; i++) {
                BrownianStateForBridgingStepStone state = new BrownianStateForBridgingStepStone(x0, t0, theData.theTrees, (i==0), template,betas[j+1],betas[j]);
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
        }
        else{
                //create the BrownianState using the previousy saved bridges:
                BrownianStateForBridgingStepStone state = new BrownianStateForBridgingStepStone(x0, t0, numSteps,finalBridges,betas[j+1],betas[j]);
                theBrownianState = state;
                java.util.HashSet<State> set = new java.util.HashSet<State>();
                set.add(state);
                initialState[0] = new GlobalState(set);
                initialState[0].setLogLikelihood(state.getLogLike());
                
                BrownianStateForBridgingStepStone state2 = new BrownianStateForBridgingStepStone(x0, t0, theData.theTrees, false, template,betas[j+1],betas[j]);
                java.util.HashSet<State> set2 = new java.util.HashSet<State>();
                set2.add(state2);
                initialState[1] = new GlobalState(set2);
                initialState[1].setLogLikelihood(state2.getLogLike());
             
        }
       
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
                //only burn in for the first chain
                if(j==1){
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
                else{
                    burnits=0;
                    i=i+2;
                }
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
                outFile = new File(args[i+1]+j+".txt");
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
       else {
                /* Unrecognized option */
                System.out.println("Bad argument: "+a+"\n");
                displaySyntaxMessage();
                System.exit(1);
            }
        }
        
       
        /* Build the sweep proposal */
        KernelWrapper myKernelWrapper = new SweepKernelWrapper("Overall proposal", kernelWrapperList);
       
        System.out.println(initialState[0].getLogLikelihood()+" p");
        /* Run the chain */
        Chain myChain = new Chain(initialState[0], initialState[1], myKernelWrapper);
        
        myChain.run(nits, burnits, thin, outFile);
        //store the bridges at the end of the run to initialise the new chain - no need for burnin then 
        finalBridges.clear();
        for(int s = 0; s<numDataPoints;s++){
        //construct the substate name
        String subStateName = "bridge_"+String.format("%03d", (s+1));
        String finalBridgeAsString = myChain.getState().getSubState("brownian_motion_parameters").getSubState(subStateName).getValueString();
        System.out.println(myChain.getState().getLogLikelihood());
        
        String[] finalBridgeAsStrings = finalBridgeAsString.split(";");
        TreeAsSplits[] finalBridgeTrees = new TreeAsSplits[numSteps+1];
        
        for(int k=0;k<numSteps+1;k++){
            finalBridgeAsStrings[k]=finalBridgeAsStrings[k]+";";
            finalBridgeTrees[k]=new TreeAsSplits(new Tree(finalBridgeAsStrings[k]));
            System.out.println(finalBridgeTrees[k].toString());
        }
        
        
        BridgeWithStepStoneLike finalBridge = new BridgeWithStepStoneLike(finalBridgeTrees,betas[j+1]);
        finalBridge.setLogLikes(t0);
        
        
        finalBridges.add(new BridgeState(finalBridge,subStateName));
        
        }
        }
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
    

    public static void main(String[] args) throws AlgorithmException {
       /*
        String outFilename = args[0];
        String SourceTreeFilename = args[1];
        String EndTreeFilename = args[2];
        Integer StartPoint = Integer.parseInt(args[3]);
        Integer EndPoint = Integer.parseInt(args[4]);
        String t0String = args[5];
        
        
        args = new String[15];
        args[0] = EndTreeFilename;
        args[1] = SourceTreeFilename;
        args[2] = t0String; // Initial squ root t_0
        args[3] = "50"; // Num steps
        args[4] = "1802"; // Seed
        args[5] = "-n";
        args[6] = "12000"; // Num interations - before thin
        args[7] = "-t";
        args[8] = "20"; // thin
        args[9] = "-b";
        args[10] = "20000";//burn-in
        args[11] = "-o";
        args[12] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/MultiplePointTest/SSQSTest";
        args[13] = "-pbg";
        args[14] = "0.02";//geometric length bridge prop
        
        //args for directly sampling from the proposal distribution
        
        
        for(int i=StartPoint; i<EndPoint;i++){
          args[4] = Integer.toString(24*i*i); 
          args2[0] = outFilename+i+"r0.txt";
          args[12] = outFilename+i+"r";
        String[] args2 = new String[7];
        args2[0] = "";//File to output the information about the proposals
        args2[1] = "300000";//number of Proposals to run
        args2[2] = args[1]; //File containing the source tree
        args2[3] = args[0]; // File containing the data set
        args2[4] = args[3]; // number of steps
        args2[5] = args[2]; // Square root dispersion
        args2[6] = "0.01"; // the first non zero value of beta_k
        ProposalsForSteppingStone.parseArgs(args2);
        parseCommandLine(args);
        }
*/
       String[] args1 = Arrays.copyOfRange(args, 0, 15);
       String[] args2 = new String[7];
        args2[0] = args[12]+"0.txt";//File to output the information about the proposals
        args2[1] = args[15];//number of Proposals to run
        args2[2] = args[1]; //File containing the source tree
        args2[3] = args[0]; // File containing the data set
        args2[4] = args[3]; // number of steps
        args2[5] = args[2]; // Square root dispersion
        args2[6] = args[16]; // the first non zero value of beta_k
       
        
        ProposalsForSteppingStone.parseArgs(args2);
        parseCommandLine(args1);
        
       
       
    }

}