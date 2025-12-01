/*
TunnelSamplingDisp
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
import MCMC.OutputStringProcessing;
import MCMC.PositiveParameter;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import MarginalLikelihoods.AuxilliaryMethods;
import static MarginalLikelihoodsWithDisp.TunnelSamplingDisp.parseCommandLine;
import cern.jet.random.tdouble.DoubleUniform;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.CategoricalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;


public class TunnelSamplingDisp {
    
        public static void displaySyntaxMessage() {
        System.out.println("Syntax: InferBrownianParams data_file initial_tree_file initial_sqrt_t0 numSteps seed\n"
                + "    [-n num_its] [-b burn_its] [-t thin] [-o outfile] \n"
                + "    [-x0sph shape scale] [-x0std shape scale] [-t0gam shape scale]\n"
                + "    [-pbg geom_fac] [-pbl len] [-pbf start end] \n"
                + "    [-ptr sigma] \n"
                + "    [-pxf sigma length] [-pxg sigma geom]\n");
                
    }
      /*
     Run the MCMC targetting the posterior distribution on t0 and bridges when calculating marginal likelihoods
     over t0 and bridges using either Chib one block or tunnel estimators  
     Then BrownianStateForBridgingTunnel outputs the necessary information for computing the estimates
        */ 

    public static void parseCommandLine(String[] args) throws AlgorithmException {
        System.out.println("Tunnel 1");
        int seed = 1;
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        int nits = 100, burnits = 0, thin = 1, numSteps = 50;
        double t0 = 0.01;
        File outFile = null;
        Integer parametrise = 0;
        String RefDistFileName ="";
        String NewT0FileName ="";
        
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
        //define the prior on t_0
        double rateForT0Prior = 18.44*((double)(initialTree.numTaxa()-3))/((double)(initialTree.numTaxa()));
        DensityCalculator t0Prior = new PositiveParameter.ExponentialPrior(rateForT0Prior);
       
       //calculate the initial value for t0
       double centralt0=AuxilliaryMethods.getFrechetVariance(x0, theData.theTrees);
       t0=centralt0;
       
       System.out.println(centralt0+" Frechet variance of data about x0");
       //get the reference distribution on t0: this is overwritten when we parametrise the reference distribution after the run
       double mu = Math.log(1.25 * centralt0);
       double sigmaRho = Math.sqrt(mu - Math.log(centralt0));
       DensityCalculator t0Dist = new PositiveParameter.LogNormalPrior(mu,sigmaRho);
       
       DensityCalculator thePrior = new BrownianStateForBridgingTunnel.TunnelBrownianStatePrior(t0Prior,t0Dist,false);

       
       System.out.println("Prior: Exponential with rate = "+String.format("%7.7f", rateForT0Prior));
       System.out.println("Reference: Lognormal with mu = "+String.format("%7.7f", mu)+" and sigma = "+String.format("%7.7f", sigmaRho));
        
        
        /* Build initial bridges */
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForBridgingTunnel theBrownianState = null; 
        ForwardStepBridge template = new BridgeWithApproxMVNLike(x0, null, numSteps);//usual bridge likelihood
        try {
            for(int i=0; i<2; i++) {
                BrownianStateForBridgingTunnel state = new BrownianStateForBridgingTunnel(x0, t0, theData.theTrees, (i==0), template,t0Prior,t0Dist,false);
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
            else if(a.equals("-prd")) {
                
                try {
                    parametrise = Integer.parseInt(args[i+1]);
                    if (parametrise==1) {
                        RefDistFileName = args[i+2];
                        NewT0FileName = args[i+3];
                    }
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with t0 reference distribution information on command line.");
                    System.exit(1);
                }
                i=i+4;
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
       
 
        /* Run the chain */
        Chain myChain = new Chain(initialState[0], initialState[1], myKernelWrapper);
        myChain.run(nits, burnits, thin, outFile);
        
        /*parametrise the reference distribution and calculate new density values for the posterior t0 sample
        *now that the we know the actual reference distribution on t0
        */
        if(parametrise==1){
            try {
                        Double[] RefDistParams = AuxilliaryMethods.getRefDistParameters(outFile);
                        File RefDistOutputFile = new File(RefDistFileName);
                        File NewT0DensFile = new File(NewT0FileName);
        //save reference dist parameters to a file       
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(RefDistOutputFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+RefDistOutputFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
            
}
             out.println("mu sigmaRho");
             out.println(RefDistParams[0]+" "+RefDistParams[1]);
             out.close();
             //if reference distribution is calculated in this way, then need to recalculate the values of the reference distribution
             //in the posterior sample:
             try {
                 FileWriter out1= new FileWriter(NewT0DensFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+NewT0DensFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
            
}
             DensityCalculator NewT0Dist = new PositiveParameter.LogNormalPrior(RefDistParams[0],RefDistParams[1]);
             ArrayList<Double> T0Posterior = AuxilliaryMethods.getT0Posterior(outFile);
             for(double Disp:T0Posterior){
                 initialState[0].getSubState("brownian_motion_parameters").getSubState("dispersion").setValueFromString(Double.toString(Disp));
                 out.println(NewT0Dist.logDensity(initialState[0], "dispersion"));
             }
            out.close();
        }
            catch (IOException ex) {
            System.out.println("Failed to get reference dist params");
            Logger.getLogger(TunnelSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }   catch (OutputStringProcessing.InsufficientInformationException ex) {
                Logger.getLogger(TunnelSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
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
        
        
        String[] args1 = new String[21];
        //arguments for the posterior MCMC sims -- may want to use different parameters for the posterior and reference dists
        args1[0] = args[0]; //data filename
        args1[1] = args[1]; //source tree filename
        args1[2] = "0.0"; // // Initial squ root t_0 - delete this option eventually
        args1[3] = args[2]; // Num steps
        args1[4] = args[3]; // Seed
        args1[5] = args[4];
        args1[6] = args[5]; // Num interations - before thin
        args1[7] = args[6];
        args1[8] = args[7]; // thin
        args1[9] = args[8];
        args1[10] = args[9]; //burn-in
        args1[11] = args[10];
        args1[12] = args[11];
        args1[13] = args[12];
        args1[14] = args[13]; //geometric length partial bridge proposal parameter
        args1[15] = args[14];
        args1[16] = args[15]; // dispersion log random walk proposal parameter
        args1[17] = args[16]; //whether to parametrise the lognormal reference distribution for t0: 1 for yes 0 for no -- advisable to do so
        args1[18] = args[17];
        args1[19] =args[18]; //file to store the parameters of the t0 reference distribution in
        args1[20] =args[19]; //file to store the new t0 densities (for the reference dist) in 



        String[] args2 = new String[19];
        //arguments for the ref dist MCMC sims -- may want to use different parameters for the posterior and reference dists
        args2[0] = args1[0];
        args2[1] = args1[1];
        args2[2] = "0.0"; //Initial squ root t_0 - delete this option eventually
        args2[3] = args1[3]; // Num steps
        args2[4] = "0.0"; //Placeholder for reference dist params
        args2[5] = "0.0"; //Placeholder for reference dist params
        args2[6] = args[20]; // Seed
        args2[7] = args[21];
        args2[8] = args[22]; // Num interations - before thin
        args2[9] = args[23];
        args2[10] = args[24]; // thin
        args2[11] = args[25];
        args2[12] = args[26];//burnin
        args2[13] = args[27];
        args2[14] = args[28];//output filename for simple indep props via MCMC
        args2[15] = args[29];
        args2[16] = args[30]; //geometric length partial bridge proposal parameter
        args2[17] = args[31]; //dispersion log random walk proposal parameter
        args2[18] = args[32];

        String[] args3 = new String[8];
        //arguments for calculating the proportion of simple proposals
        args3[0] = args1[0];
        args3[1] = args1[1];
        args3[2] = args[33]; // OutFileName
        args3[3] = "0.0"; //Placeholder for reference dist params
        args3[4] = "0.0"; //Placeholder for reference dist params
        args3[5] = args[34]; //numDisps
        args3[6] = args[35]; //numProps
        args3[7] = args1[3]; //num steps
        
            System.out.println("Starting posterior sims");
            parseCommandLine(args1);//simulate from the posterior
            try {
            Double[] RefDistParams = AuxilliaryMethods.getRefDistParameters(args1[12]);
            args2[4]=Double.toString(RefDistParams[0]);
            args2[5]=Double.toString(RefDistParams[1]);
            args3[3]=Double.toString(RefDistParams[0]);
            args3[4]=Double.toString(RefDistParams[1]);
            //if reference distribution is calculated in this way, then need to recalculate the values of the reference distribution
            //in the posterior sample:

        } catch (IOException ex) {
            System.out.println("Failed to get reference dist params");
            Logger.getLogger(TunnelSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }
          
         
        try {
            System.out.println("Starting prop simple sims");
            System.out.println(ProportionOfSimpleProposalsDispForTunnel.propSimple(args3));
        } catch (IOException ex) {
            System.out.println("Error calculating proportion of simple proposals");
            Logger.getLogger(TunnelSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }

            System.out.println("Starting ref dist sims");
            TunnelSamplingProposal.simConditionalProposal(args2);//simulate from the reference distribution
         
       
    }

}


/* //old main class for looping through many proposals
   public static void main(String[] args) throws AlgorithmException {
        
        //run the MCMC to get the posterior sample
    int Startpt=0;
        int Endpt =100;
        
        try{
            Startpt=Integer.parseInt(args[0]);
            Endpt=Integer.parseInt(args[1]);
        }
        catch(ArrayIndexOutOfBoundsException e){
            System.out.println("Startpt and endpoint not supplied with args, resorting to default");
        }    
    args = new String[21];
    //arguments for the posterior MCMC sims -- may want to use different parameters for the posterior and reference dists
        args[0] = "/home/c1032934/Documents/Netbeans/TopInf20241024/MarginalLikelihoods/4TaxaWith/DataPoints.txt";
        args[1] = "/home/c1032934/Documents/Netbeans/TopInf20241024/MarginalLikelihoods/4TaxaWith/x1.txt";
        //args[0] = "/home/c1032934/Documents/Netbeans/TopInf20240115/MarginalLikelihoods/10Taxa20240222/10DataPoints.txt";
        //args[1] = "/home/c1032934/Documents/Netbeans/TopInf20240115/MarginalLikelihoods/10Taxa20240222/SourceTree1.txt";
        args[2] = "0.0"; // // Initial squ root t_0 - delete this option in final version
        args[3] = "20"; // Num steps
        args[4] = "344"; // Seed
        args[5] = "-n";
        args[6] = "250000"; // Num interations - before thin
        args[7] = "-t";
        args[8] = "20"; // thin
        args[9] = "-b";
        args[10] = "10000"; //burn-in
        args[11] = "-o";
        args[12] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/MultiplePointTest/Post.txt";
        args[13] = "-pbg";
        args[14] = "0.01";
        args[15] = "-ptr";
        args[16] = "0.5";
        args[17] = "-prd"; //whether to parametrise the lognormal reference distribution for t0: 1 for yes 0 for no -- advisable to do so
        args[18] = "1";
        args[19] =""; //file to store the parameters in
        args[20] =""; //file to store the new t0 densities in 

        
       
        String[] args2 = new String[19];
    //arguments for the ref dist MCMC sims -- may want to use different parameters for the posterior and reference dists
        args2[0] = args[0];
        args2[1] = args[1];
        args2[2] = "0.0"; //Initial squ root t_0 - delete this option in final version
        args2[3] = args[3]; // Num steps
        args2[4] = "0.0"; //Placeholder for reference dist params
        args2[5] = "0.0"; //Placeholder for reference dist params
        args2[6] = "897"; // Seed
        args2[7] = "-n";
        args2[8] = "250000"; // Num interations - before thin
        args2[9] = "-t";
        args2[10] = "20"; // thin
        args2[11] = "-b";
        args2[12] = "10000";//burnin
        args2[13] = "-o";
        args2[14] = "";//output filename for simple indep props via MCMC
        args2[15] = "-pbg";
        args2[16] = "0.01";
        args2[17] = "-ptr";
        args2[18] = "0.5";
        
        String[] args3 = new String[8];
    //arguments for calculating the proportion of simple proposals
        args3[0] = args[0];
        args3[1] = args[1];
        args3[2] = ""; // OutFileName
        args3[3] = "0.0"; //Placeholder for reference dist params
        args3[4] = "0.0"; //Placeholder for reference dist params
        args3[5] = "200"; //numDisps
        args3[6] = "500"; //numProps
        args3[7] = args[3]; //num steps
        
        
        for(int i=Startpt; i<Endpt;i++){
        
            args[4] = Integer.toString(i*i*42+9);
            args2[6] = Integer.toString(i*i*52+9);
            
            args[12] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_Post"+i+".txt";
            args[19] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_T0RefDistParams"+i+".txt";
            args3[2] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_PropSimple"+i+".txt";
            args2[14] ="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_Props"+i+".txt";
            args[20] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_NewT0RefDist"+i+".txt";

            parseCommandLine(args);//simulate from the posterior
            try {
            Double[] RefDistParams = AuxilliaryMethods.getRefDistParameters(args[12]);
            args2[4]=Double.toString(RefDistParams[0]);
            args2[5]=Double.toString(RefDistParams[1]);
            args3[3]=Double.toString(RefDistParams[0]);
            args3[4]=Double.toString(RefDistParams[1]);
            System.out.println(args2[4]+" mu "+args2[5]+" sigmaRho");
            //if reference distribution is calculated in this way, then need to recalculate the values of the reference distribution
            //in the posterior sample:

        } catch (IOException ex) {
            System.out.println("Failed to get reference dist params");
            Logger.getLogger(TunnelSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }
          
         
        try {
            System.out.println(ProportionOfSimpleProposalsDispForTunnel.propSimple(args3));
        } catch (IOException ex) {
            System.out.println("Error calculating proportion of simple proposals");
            Logger.getLogger(TunnelSamplingDisp.class.getName()).log(Level.SEVERE, null, ex);
        }


            TunnelSamplingProposal.simConditionalProposal(args2);//simulate from the reference distribution
        }
        

       
       
    }
*/