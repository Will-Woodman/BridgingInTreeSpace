/*
SteppingStoneSamplerDisp
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


package MarginalLikelihoodsWithDisp;

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
import java.io.File;
import java.util.ArrayList;
import simulation.CategoricalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;
import bridge.BridgeState;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import MarginalLikelihoods.AuxilliaryMethods;

/*
Run the MCMC needed to compute a stepping stone sampling estimate of the marginal likelihood when t0 is unknown.
Use the quasistatic method so burnin is only needed to be run once
*/

public class SteppingStoneSamplingDisp {
    
        public static void displaySyntaxMessage() {
        System.out.println("Syntax: SteppingStoneSampling data_file initial_tree_file initial_sqrt_t0 numSteps"
                + "    ref_dist_mu ref_dist_rho seed\n"
                + "    [-n num_its] [-b burn_its] [-t thin] [-o outfile] \n"
                + "    [-pbg geom_fac] [-pbl len] [-pbf start end] \n"
                + "    [-ptr sigma] \n");
                
    }

        
    public static void parseCommandLine(String[] args) throws AlgorithmException {
       

        Double[] betas=new Double[101];
        /*
        set the values for beta_k at which to sample from the unnormalised distribution
        here using equally spaced beta_k's.
        */
        for(int i=0;i<101;i++){
    betas[i]=(double)i/100;
}
        int seed = 1;
        double mu=0.0;
        double sigmaRho=0.0;
        
        int nits = 100, burnits = 0, thin = 1, numSteps = 50;
        double t0 = 0.01;
        double centralt0 =0.01;
        double finalt0=0.01;
        File outFile = null;
        
        /* Parse arguments */
        if ((args.length < 4)||(args.length > 23)) {
            displaySyntaxMessage();
            System.exit(1);
        }

        int n = args.length;
        String dataFileName = args[0];
        String initialTreeFileName = args[1];
        
        /* Read in data */
        File inFile = new File(dataFileName);
        TreeAsSplitsDataSet theData = null;
        int numBridges =0;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
            numBridges = theData.numTrees;
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
        if("FV".equals(args[2])){//usually input FV so that the initial value for t0 is the Frechet variance about x0
         centralt0=AuxilliaryMethods.getFrechetVariance(x0, theData.theTrees);
         System.out.println(centralt0+"FrechVar");
         t0=centralt0;
        }
        else{
        try {
            double sd = Double.parseDouble(args[2]);
            t0 = sd*sd;
            System.out.println("Initial t0 = "+String.format("%7.7f", t0));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read initial t0 from command line.");
            System.exit(1);
        }
    }
        try {
            numSteps = Integer.parseInt(args[3]);
            System.out.print("Number of random walk steps in bridge = "+String.format("%d%n", numSteps));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of RW steps from command line.");
            System.exit(1);
        }
        //read in the parameters of the reference distribution:
        try {
            mu = Double.parseDouble(args[4]);
            sigmaRho = Double.parseDouble(args[5]);
            System.out.print("Reference dist params = "+mu+" "+sigmaRho);
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read reference dist params from command line.");
            System.exit(1);
        }
        
        try {
            seed = Integer.parseInt(args[6]);
            System.out.print("Seed = "+String.format("%d%n", seed));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read seed from command line.");
            System.exit(1);
        }
        Random.setEngine(seed);
        TreeResolver.resolveTree(x0, new DoubleUniform(Random.getEngine()));
        
        
         ArrayList<BridgeState> finalBridges = new ArrayList();
         
       /*
       if the reference distribution has not been set, set from the Frechet variance of the data:
       recommended to always set it and not do the following:
       */
       if(mu==0){                 
          centralt0=AuxilliaryMethods.getFrechetVariance(x0, theData.theTrees);
          System.out.println(centralt0+"FrechVar");
          mu = Math.log(1.25 * centralt0);
          sigmaRho = Math.sqrt(mu - Math.log(centralt0));
          } 
       //initial value of t0 is the mean of the reference distribution -- delete the option to have another value
        t0=Math.exp(mu +sigmaRho*sigmaRho );
        //usual rate for the prior on t0
        double rateForT0Prior = 18.44*((double)(initialTree.numTaxa()-3))/((double)(initialTree.numTaxa())); 
        /*
        Now loop through the beta_k's: leave the creation of proposals and the reading in of number iterations etc
        within the loop so created fresh each time + allows for small tweaks so that we could specifiy different
        parameters for different beta_k etc
        */
        for(int j=0; j<100; j++){
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForBridging theBrownianState = null; 
        BridgeWithStepStoneLike template = new BridgeWithStepStoneLike(x0, null, numSteps,betas[j]);
        
        
        DensityCalculator t0Prior = new PositiveParameter.ExponentialPrior(rateForT0Prior);

       System.out.println(mu+" "+sigmaRho);
       DensityCalculator t0Dist = new PositiveParameter.LogNormalPrior(mu,sigmaRho);
       
        //define the stepping stone prior with the current value of beta
        BrownianStateForBridgingStepStone.SteppingStoneBrownianStatePrior thePrior = new BrownianStateForBridgingStepStone.SteppingStoneBrownianStatePrior(t0Prior,t0Dist,betas[j]); 
        
        if(j<1){//then create new bridges
        try {
            for(int i=0; i<2; i++) {
                BrownianStateForBridgingStepStone state = new BrownianStateForBridgingStepStone(x0, t0, theData.theTrees, (i==0), template,betas[j+1],betas[j],thePrior);
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
        else{//read in final bridges and t0 from the previous beta
            
                BrownianStateForBridgingStepStone state = new BrownianStateForBridgingStepStone(x0, finalt0, numSteps,finalBridges,betas[j+1],betas[j],thePrior);
                theBrownianState = state;
                java.util.HashSet<State> set = new java.util.HashSet<State>();
                set.add(state);
                initialState[0] = new GlobalState(set);
                initialState[0].setLogLikelihood(state.getLogLike());
                
                BrownianStateForBridgingStepStone state2 = new BrownianStateForBridgingStepStone(x0, finalt0, theData.theTrees, false, template,betas[j+1],betas[j]);
                java.util.HashSet<State> set2 = new java.util.HashSet<State>();
                set2.add(state2);
                initialState[1] = new GlobalState(set2);
                initialState[1].setLogLikelihood(state2.getLogLike());
             
        }
       
        // Loop thru arguments
        int i=7;
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
                if(j<1){
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
       else {
                /* Unrecognized option */
                System.out.println("Bad argument: "+a+"\n");
                displaySyntaxMessage();
                System.exit(1);
            }
        }
        
       
        /* Build the sweep proposal */
        KernelWrapper myKernelWrapper = new SweepKernelWrapper("Overall proposal", kernelWrapperList);
       
        System.out.println(initialState[0].getLogLikelihood()+" initial log likelihood at beta= "+betas[j]);
        /* Run the chain */
        Chain myChain = new Chain(initialState[0], initialState[1], myKernelWrapper);
        myChain.run(nits, burnits, thin, outFile);
        
        //store the bridges at the end of the run to intitialise the new chain - no need for burnin then 
        finalBridges.clear();
        for(int l=0;l<numBridges;l++){
        String bridgeName="bridge_"+String.format("%03d", (l+1));
        //System.out.println(bridgeName);
        String finalBridgeAsString = myChain.getState().getSubState("brownian_motion_parameters").getSubState(bridgeName).getValueString();
        //System.out.println(myChain.getState().getLogLikelihood());
        finalt0 = Double.parseDouble(myChain.getState().getSubState("brownian_motion_parameters").getSubState("dispersion").getValueString());
        
        String[] finalBridgeAsStrings = finalBridgeAsString.split(";");
        TreeAsSplits[] finalBridgeTrees = new TreeAsSplits[numSteps+1];
        
        for(int k=0;k<numSteps+1;k++){
            finalBridgeAsStrings[k]=finalBridgeAsStrings[k]+";";
            finalBridgeTrees[k]=new TreeAsSplits(new Tree(finalBridgeAsStrings[k]));
            //System.out.println(finalBridgeTrees[k].toString()); - uncomment if you want to print out the final bridges for each beta
        }
        
        
        BridgeWithStepStoneLike finalBridge = new BridgeWithStepStoneLike(finalBridgeTrees,betas[j+1]);

        finalBridge.setLogLikes(finalt0);
        
        
        finalBridges.add(new BridgeState(finalBridge,bridgeName));
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
    
     
     
    

    public static void main(String[] args) throws AlgorithmException, IOException {
        
        //run the MCMC for a short time to get the posterior sample for which to estimate the reference distribution parameters
        String[] args2 = new String[17];
        args2[0] = args[0];
        args2[1] = args[1];
        args2[2] = args[2]; // Initial squ root t_0 - set to zero to use the Frechet variance -- recommended to use FV
        args2[3] = args[3]; // Num steps
        args2[4] = args[4]; // Seed
        args2[5] = args[5];
        args2[6] = args[6]; // Num interations - before thin
        args2[7] = args[7];
        args2[8] = args[8]; // thin
        args2[9] = args[9];
        args2[10] = args[10]; //burn-in
        args2[11] = args[11];
        args2[12] = args[12];
        args2[13] = args[13];
        args2[14] = args[14]; //geometric length bridge prop
        args2[15] = args[15];
        args2[16] = args[16]; //t0 random walk proposal
        
        
         String[] args3 = new String[8];
    //arguments for calculating the proportion of simple proposals
        args3[0] = args2[0];
        args3[1] = args2[1];
        args3[2] = args[17]; // OutFileName
        args3[3] = "0.0"; //Placeholder for reference dist params
        args3[4] = "0.0"; //Placeholder for reference dist params
        args3[5] = args[18]; //numDisps
        args3[6] = args[19]; //numProps
        args3[7] = args2[3]; //num steps
        
        
        String [] args1 = new String[19];
            args1[0] = args2[0];
            args1[1] = args2[1];
            args1[2] = args[20]; // write FV to use the Frechet variance about x0 or otherwise specidy a number -- recommended to use FV
            args1[3] = args2[3]; // Num steps
            args1[4] ="0.0"; //placeholder for reference dist parameters
            args1[5] ="0.0"; //placeholder for reference dist parameters
            args1[6] = args[21]; // Seed
            args1[7] = args[22];
            args1[8] = args[23]; // Num interations - before thin
            args1[9] = args[24];
            args1[10] = args[25]; // thin
            args1[11] = args[26];
            args1[12] = args[27];//burn-in
            args1[13] = args[28];
            args1[14] = args[29];
            args1[15] = args[30];
            args1[16] = args[31];//geometric length bridge prop
            args1[17] = args[32];
            args1[18] = args[33];//t0 random walk proposal

          SteppingStonePosteriorExploration.simPosterior(args2);//run the short run from the posterior to parametrise t0 ref dist
          try {
            Double[] RefDistParams = AuxilliaryMethods.getRefDistParameters(args2[12]);//extract the ref dist params from short posterior run file
            args[4]=Double.toString(RefDistParams[0]);
            args[5]=Double.toString(RefDistParams[1]);
            args3[3]=Double.toString(RefDistParams[0]);
            args3[4]=Double.toString(RefDistParams[1]);
            System.out.println(args[4]+" mu "+args[5]+" sigmaRho");
        } catch (IOException ex) {
            System.out.println("Failed to get reference dist params");

            Logger.getLogger(SteppingStoneSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }
        
         ProportionOfSimpleProposalsDispForTunnel.propSimple(args3);
          parseCommandLine(args1);//run the stepping stone sampler
        }
        

}

/*
//Old main class for looping through many estimates
    public static void main(String[] args) throws AlgorithmException, IOException {
        
        int Startpt=0;
        int Endpt =100;
        
        try{
            Startpt=Integer.parseInt(args[0]);
            Endpt=Integer.parseInt(args[1]);
        }
        catch(ArrayIndexOutOfBoundsException e){
            System.out.println("Startpt and endpoint not supplied with args, resorting to default");
        }    
        
        //run the MCMC for a short time to get the posterior sample for which to estimate the reference distribution parameters
        String[] args2 = new String[17];
        args2[0] = "/home/c1032934/Documents/Netbeans/TopInf20241024/MarginalLikelihoods/4TaxaWith/DataPoints.txt";
        args2[1] = "/home/c1032934/Documents/Netbeans/TopInf20241024/MarginalLikelihoods/4TaxaWith/x1.txt";
        args2[2] = "0"; // Initial squ root t_0 - set to zero to use the Frechet variance -- recommended to use FV
        args2[3] = "20"; // Num steps
        args2[4] = "490"; // Seed
        args2[5] = "-n";
        args2[6] = "30000"; // Num interations - before thin
        args2[7] = "-t";
        args2[8] = "20"; // thin
        args2[9] = "-b";
        args2[10] = "5000";//burn-in
        args2[11] = "-o";
        args2[12] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/5Taxa/SSQSTest20231127.txt";
        args2[13] = "-pbg";
        args2[14] = "0.01";//geometric length bridge prop
        args2[15] = "-ptr";
        args2[16] = "0.5";//t0 random walk proposal
        
        
         String[] args3 = new String[8];
    //arguments for calculating the proportion of simple proposals
        args3[0] = args2[0];
        args3[1] = args2[1];
        args3[2] = ""; // OutFileName
        args3[3] = "0.0"; //Placeholder for reference dist params
        args3[4] = "0.0"; //Placeholder for reference dist params
        args3[5] = "200"; //numDisps
        args3[6] = "500"; //numProps
        args3[7] = args2[3]; //num steps
        
        
    args = new String[19];
        args[0] = args2[0];
        args[1] = args2[1];
        args[2] = "FV"; // write FV to use the Frechet variance about x0 or otherwise specidy a number -- recommended to use FV
        args[3] = args2[3]; // Num steps
        args[4] ="0.0"; //placeholder for reference dist parameters
        args[5] ="0.0"; //placeholder for reference dist parameters
        args[6] = "490"; // Seed
        args[7] = "-n";
        args[8] = "3000"; // Num interations - before thin
        args[9] = "-t";
        args[10] = "20"; // thin
        args[11] = "-b";
        args[12] = "5000";//burn-in
        args[13] = "-o";
        args[14] = "";
        args[15] = "-pbg";
        args[16] = "0.01";//geometric length bridge prop
        args[17] = "-ptr";
        args[18] = "0.5";//t0 random walk proposal
       
        for(int i=Startpt; i<Endpt;i++){
          args[6] = Integer.toString(47658*i+5678); 
          args2[4] = Integer.toString(746373*i+5678); 
          args2[12]="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/StepStone/x1_"+i+"Exploration.txt";
          args[14] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/StepStone/x1_"+i+"r";
          args3[2] ="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/StepStone/x1_"+i+"PropSimple.txt";
          SteppingStonePosteriorExploration.simPosterior(args2);//run the short run from the posterior to parametrise t0 ref dist
          try {
            Double[] RefDistParams = AuxilliaryMethods.getRefDistParameters(args2[12]);//extract the ref dist params from short posterior run file
            args[4]=Double.toString(RefDistParams[0]);
            args[5]=Double.toString(RefDistParams[1]);
            args3[3]=Double.toString(RefDistParams[0]);
            args3[4]=Double.toString(RefDistParams[1]);
            System.out.println(args[4]+" mu "+args[5]+" sigmaRho");
        } catch (IOException ex) {
            System.out.println("Failed to get reference dist params");

            Logger.getLogger(SteppingStoneSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }
            //no need to do the below on 4 taxa:
          System.out.println(ProportionOfSimpleProposalsDispForTunnel.propSimple(args3));
          parseCommandLine(args);//run the stepping stone sampler
        }
        


       
       
    }

}
*/