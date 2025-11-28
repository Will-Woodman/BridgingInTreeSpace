/*
ChibTwoBlockMultStar
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


import bridge.*;
import MCMC.Chain;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.KernelWrapper;
import MCMC.PositiveParameter;
import MCMC.RealParameter.RealParameterProposal;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import cern.jet.random.tdouble.DoubleUniform;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import simulation.CategoricalDistribution;
import simulation.LogNormalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;
import MarginalLikelihoods.AuxilliaryMethods;
/*
Sampler for the multiple t_0^* version of the Chib two block estimate of marginal likelihood
Simulates from the full posterior, the marginal posterior of bridges for each fixed t_0^* and the independence
proposals with t_0=t_0^*
*/

public class ChibTwoBlockMultStar {
    
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
        int numT0Rpts =1;
        double delta=0.1;//proposal parameter used in the estimation
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
        try {
            numT0Rpts = Integer.parseInt(args[5]);
            System.out.print("Number of times to repeat t0 proposal for the pi*rho dist = "+numT0Rpts);
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of times to repeat t0 proposal for the pi*rho dist from command line.");
            System.exit(1);
        }
        try {
            delta = Double.parseDouble(args[6]);
            System.out.println("t0 proposal param used in estimation = "+String.format("%7.7f", delta));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read t0 proposal param from command line.");
            System.exit(1);
        }
        Random.setEngine(seed);
        TreeResolver.resolveTree(x0, new DoubleUniform(Random.getEngine()));
        //usual prior on t0
        double rateForT0Prior = 18.44*((double)(initialTree.numTaxa()-3))/((double)(initialTree.numTaxa()));
       DensityCalculator t0Prior = new PositiveParameter.ExponentialPrior(rateForT0Prior);
       //proposal used in the calculation 
       RealParameterProposal t0Dist = new PositiveParameter.LogNormalProposal(delta);
       
       DensityCalculator thePrior = new BrownianStateForBridgingTunnel.ChibBrownianStatePrior(t0Prior,t0Dist);
        
       
       System.out.println("Prior: Exponential with rate = "+String.format("%7.7f", rateForT0Prior));
       System.out.println("Proposal for calculation: Lognormal with variance param: "+String.format("%7.7f", delta));
        /* Build initial bridges */
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForBridgingChibSecondBlockMultStar theBrownianState = null; 
        ForwardStepBridge template = new BridgeWithApproxMVNLike(x0, null, numSteps);
        try {
            for(int i=0; i<2; i++) {
                BrownianStateForBridgingChibSecondBlockMultStar state = new BrownianStateForBridgingChibSecondBlockMultStar(x0, t0, theData.theTrees, (i==0), template,t0Prior,t0Dist,t0,numT0Rpts);
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
                    //make a new method to do this
                    makeProposalData(propsFile, numProps, x0, theData, numSteps,t0);                   
                }
                catch (NumberFormatException anErr) {
                    System.out.println("Error with numProps points for bridge proposal on command line.");
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
       
 
        /* Run the chain */
        Chain myChain = new Chain(initialState[0], initialState[1], myKernelWrapper);
        myChain.run(nits, burnits, thin, outFile);

    
    }
    
        private static void makeProposalData(File propFile, int numProp,TreeAsSplits startTree, TreeAsSplitsDataSet Data, int m,double t0) throws AlgorithmException{
        //is it best to just read the results out to a file or just send the array back to the main method?
        ArrayList<ForwardStepBridge> theBridges = new ArrayList();
        int numBridges= Data.numTrees;
        double[][][] outputs = new double[numProp][numBridges][2];
        int nPrime =startTree.getNumTaxa()-3;
        for(int i=0; i<numBridges;i++){
            ForwardStepBridge theBridge = new BridgeWithApproxMVNLike(startTree, Data.theTrees.get(i), m);
            theBridges.add(theBridge);
            }
        for(int i=0; i<numBridges;i++){ 
        double logLike = 0;
        double propDens = 0;
        for(int j=0; j<numProp;j++){
            logLike = 0;
            propDens = 0;

            try {
            theBridges.get(i).MargLikeIndependenceProposal(t0);
            propDens += theBridges.get(i).computeMargLikeIndependenceLogDensity(t0)-Math.log(2*Math.PI)*nPrime*(m-1)/2;
            }
            catch(treebase.AlgorithmException anErr){  
                logLike = Double.NEGATIVE_INFINITY;
            }
            logLike += theBridges.get(i).getTotalLogLike();
            //propDens += theBridges.get(i).computeIndependenceLogDensity(t0)-Math.log(2*Math.PI)*nPrime*(m-1)/2;
           if(propDens==Double.NEGATIVE_INFINITY){
               System.out.println("Why?");
           }
            
            
            outputs[j][i][1] = logLike;
            outputs[j][i][0] =propDens;
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

//print the header
for(int i=0;i<numBridges;i++){    
String bridgeName ="bridge_"+String.format("%03d", (i+1));    
out.print(bridgeName+"PropDens "+bridgeName+"LogLike ");

}
out.println();

for(int j=0 ; j<numProp ; j++){
    for(int i=0;i<numBridges;i++){
out.print(outputs[j][i][0]+" "+outputs[j][i][1]+" ");
    }
out.println();
    } 
 out.close();  
        
        
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
                //arguments for the initial posterior run: short run to choose reasonable values for the fixed values of dispersion
        int numT0Stars = 4;
        Double[] t0Stars = new Double[numT0Stars];
        String[] args3 = new String[17];
        args3[0] = "/home/c1032934/Documents/Netbeans/TopInf20241024/MarginalLikelihoods/4TaxaWith/DataPoints.txt";
        args3[1] = "/home/c1032934/Documents/Netbeans/TopInf20241024/MarginalLikelihoods/4TaxaWith/x1.txt";
        args3[2] = "0.3"; // Initial squ root t_0
        args3[3] = "20"; // Num steps
        args3[4] = "344"; // Seed
        args3[5] = "-n";
        args3[6] = "15000"; // Num interations - before thin
        args3[7] = "-t";
        args3[8] = "20"; // thin
        args3[9] = "-b";
        args3[10] = "5000";
        args3[11] = "-o";
        args3[12] = "";//full posterior output filename (short run)
        args3[13] = "-pbg";
        args3[14] = "0.01";
        args3[15] = "-ptr";
        args3[16] = "0.3";

        
        
    
    args = new String[20];
    //arguments for the MCMC sims from the marginal posterior and direct sims from indep prop for fixed dispersion
        args[0] = args3[0];
        args[1] = args3[1];
        args[2] = "0.3"; // Initial squ root t_0
        args[3] = args3[3]; // Num steps
        args[4] = "344"; // Seed
        args[5] = "2"; //number of times to repeat the t0 proposal in the pi*q dist
        args[6] = "0.2";//parameter for the proposal used in the estimation
        args[7] = "-n";
        args[8] = "30000"; // Num interations - before thin
        args[9] = "-t";
        args[10] = "20"; // thin
        args[11] = "-b";
        args[12] = "5000";
        args[13] = "-o";
        args[14] = "";//conditional posterior output filename
        args[15] = "-pbg";
        args[16] = "0.01";
        args[17] = "-numProps";
        args[18] = "50000"; //number of proposals for the denominator
        args[19] = ""; //file to put the proposal data in

        
       
        String[] args2 = new String[18];
        //arguments for the main MCMC sims from the full posterior
        args2[0] = args3[0];
        args2[1] = args3[1];
        args2[2] = "0.3"; // Initial squ root t_0
        args2[3] = args[3]; // Num steps
        args2[4] = "897"; // Seed
        args2[5] = args[6];//parameter for the proposal used in the estimation -- does not have to be same as one used in the MCMC sims
        args2[6] = "-n";
        args2[7] = "60000"; // Num interations - before thin
        args2[8] = "-t";
        args2[9] = "20"; // thin
        args2[10] = "-b";
        args2[11] = "5000";
        args2[12] = "-o";
        args2[13] = "";//full posterior output filename
        args2[14] = "-pbg";
        args2[15] = "0.01";
        args2[16] = "-ptr";
        args2[17] = "0.3";
       
        
        
        //file to put the proposal data in
        
        
        
        for(int i=0; i<100;i++){
          args[4] = Integer.toString(i*47+9);  
          args2[4] = Integer.toString(i*49+9);  
          args3[4] = Integer.toString(i*422+9);  
          args2[13] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_PostOut"+i+".txt";
          args3[12] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_PostOutShortRun"+i+".txt";
          //args[17] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/4TaxaWIthDisp/FixedDispPropDataTwoBlockTestDispTwoPoints"+i+".txt";
          SteppingStonePosteriorExploration.simPosterior(args3);
          Double[] RefDistParams = AuxilliaryMethods.getRefDistParameters(args3[12]);
            //args[4]=Double.toString(RefDistParams[0]);
            //args[5]=Double.toString(RefDistParams[1]);
            //System.out.println(args[4]+" mu "+args[5]+" sigmaRho");
            LogNormalDistribution margPostEst = new LogNormalDistribution(RefDistParams[0],RefDistParams[1]);
            String L="";
            for(int j=1; j <=numT0Stars;j++){
            args[14] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_FixedDisp"+j+"PostOut"+i+".txt";
            args[19] = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_FixedDisp"+j+"PropData"+i+".txt";
            double quant = 0.2*j;
            System.out.println(quant);
            t0Stars[j-1]=Math.sqrt(margPostEst.quantile(quant));
            System.out.println(t0Stars[j-1]*t0Stars[j-1]+" t0Star");
            args[2]=Double.toString(t0Stars[j-1]);
            L+=args[2]+" ";
            parseCommandLine(args);
            }
          args2[2]=L;
          ChibFirstBlockMultStar.simFullPosterior(args2);
        }
        

       
       
    }

}