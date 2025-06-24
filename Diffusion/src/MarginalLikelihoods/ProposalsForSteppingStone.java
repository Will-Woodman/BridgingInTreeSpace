/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package MarginalLikelihoods;

import bridge.BridgeWithApproxMVNLike;
import bridge.ForwardStepBridge;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
 *Used to directly simulate the independence proposals (beta=0) for stepping stone sampling
 * Output individual values for each of the bridges so that marginal likelihoods can be estimated
 * individually for each data point
 */
public class ProposalsForSteppingStone {
public static void parseArgs(String[] args) throws AlgorithmException {
        int n = args.length;
        int numSteps=20;
        String dataFileName = args[3];
        String initialTreeFileName = args[2];
        String outFileName = args[0];
        Double beta=0.01;
        int numProps=10000;
        double t0=1.0;
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
        
        /* Create the output file */
        File outFile = new File(outFileName);
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
        /*Read in the fixed value of t0 */
        try {
            double sd = Double.parseDouble(args[5]);
            t0 = sd*sd;
            System.out.println("t0 for direct proposals = "+String.format("%7.7f", t0));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read t0 from command line for proposals.");
            System.exit(1);
        }
        try {
            numSteps = Integer.parseInt(args[4]);
            System.out.print("Number of random walk steps in bridge = "+String.format("%d%n", numSteps));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of RW steps from command line.");
            System.exit(1);
        }
        try {
            numProps = Integer.parseInt(args[1]);
            System.out.print("Number of indpendence proposals to run = "+String.format("%d%n", numProps));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of independence proposals from command line.");
            System.exit(1);
        }
        try {
            beta = Double.parseDouble(args[6]);
            System.out.println("First beta value = "+String.format("%7.7f", beta));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read first beta value from command line.");
            System.exit(1);
        }
        System.out.println(numProps+" "+numSteps+beta);
        makeProposalData(outFile,numProps,x0,theData,numSteps,t0,beta);
    
}    
    
    private static void makeProposalData(File propFile, int numProp,TreeAsSplits startTree, TreeAsSplitsDataSet Data, int m,double t0,double beta) throws AlgorithmException{
        ArrayList<ForwardStepBridge> theBridges = new ArrayList();
        int numBridges = Data.numTrees;
        double[][] outputs = new double[numProp][numBridges];
        int nPrime =startTree.getNumTaxa()-3;
        for(int i=0; i<numBridges;i++){
            ForwardStepBridge theBridge = new BridgeWithApproxMVNLike(startTree, Data.theTrees.get(i), m);
            theBridges.add(theBridge);
            }
        double logLike = 0;
        double propDens = 0;
        for(int j=0; j<numProp;j++){
        for(int i=0; i<numBridges;i++){ 
            logLike = 0;
            propDens = 0;

            try {
            theBridges.get(i).MargLikeIndependenceProposal(t0);
            }
            catch(treebase.AlgorithmException anErr){  
                /*If the proposla doesn't produce a full bridge because it exited due to a non-simple geodesic
            then likelihood is -infinity
                */
                logLike = Double.NEGATIVE_INFINITY;
                //arbitrary value for the proposal density
                propDens =1;
            }
           
            if(propDens!=1){
            propDens = theBridges.get(i).computeMargLikeIndependenceLogDensity(t0)-Math.log(Math.sqrt(2*Math.PI))*(m-1)*nPrime;
            logLike = theBridges.get(i).getTotalLogLike();
            }
         
           if(propDens==Double.NEGATIVE_INFINITY){
               //This should not happen -- and it does not seem to in practice
               System.out.println("Strange behaviour, a proposed bridge has zero proposal density");
           }
             //calculate the outputs needed for the estimate
             outputs[j][i] = beta*(logLike-propDens);
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

        String L="";
        for(int k=0 ;k<numBridges;k++){
        if(k==numBridges-1) L+="Ratio"+k;
        else L+="Ratio"+k+" ";
        }
        out.println(L);
        for(int j=0 ; j<numProp ; j++){
            for(int k=0 ;k<numBridges;k++){  
                if(k==numBridges-1){
                  if(Double.isFinite(outputs[j][k])) out.println(outputs[j][k]);
                  else out.println("NA"); //if the log likelihood is negative infinity, output some string
            }
               else{
                 if(Double.isFinite(outputs[j][k])) out.print(outputs[j][k]+" ");
                 else out.print("NA "); //if the log likelihood is negative infinity, output some string   
                }
            }

     
         } 
     
       
        
          
    out.close();  
        
    }

    

}
