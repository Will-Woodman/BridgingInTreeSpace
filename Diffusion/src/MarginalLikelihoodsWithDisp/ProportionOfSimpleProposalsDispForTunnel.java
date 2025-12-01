/*
ProportionOfSimpleProposalsDispForTunnel
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

import bridge.BridgeWithApproxMVNLike;
import bridge.ForwardStepBridge;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import simulation.LogNormalDistribution;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
 Estimate the normalising constant of the reference distribution for the tunnel and Chib one block estimators
 * that include t0(See paper)
 */
public class ProportionOfSimpleProposalsDispForTunnel {
    /*
    parse the arguments passed from the tunnel sampling file
    */
    public static double propSimple(String[] args) throws AlgorithmException, IOException {
        if(args.length!=8){
            System.out.println("Wrong number of arguments provided for estimating proportion of simple proposals");
            return 0;
        }
        
        int numDisps =0;
        int numProps=0;
        int m=20;//default number of steps
        try {
            numDisps = Integer.parseInt(args[5]);
            System.out.print("numDispsPropSimple= "+String.format("%d%n", numDisps));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of dispersions for which to estimate the proportion of simple proposals.");
            System.exit(1);
        }
        try {
            numProps = Integer.parseInt(args[6]);
            System.out.print("numPropsPropSimple = "+String.format("%d%n", numProps));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of proposals to use for the estimate of the proportion of simple proposals.");
            System.exit(1);
        }
        try {
            m = Integer.parseInt(args[7]);
            System.out.print("Number of steps in prop simple sims = "+String.format("%d%n", m));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of steps in prop simple sims.");
            System.exit(1);
        }
        
        
        double[] disp=new double[numDisps];
        double[] refDensities = new double[numDisps];
        
        String initialTreeFileName =args[1];
        String dataFileName=args[0];
        String outFileName=args[2];
        
        
        Double[] RefDistParams= new Double[2];
        try {
            RefDistParams[0] = Double.parseDouble(args[3]);
            RefDistParams[1] = Double.parseDouble(args[4]);
            //System.out.print("Ref dist params for prop simple = "+String.format("%7.7f", RefDistParams[0])+" "+String.format("%7.7f", RefDistParams[1]));
        }
        catch (NumberFormatException anErr) {
            System.out.println("Unable to read number of proposals to use for the estimate of the proportion of simple proposals.");
            System.exit(1);
        }
          

        // Read in data
        File inFile = new File(dataFileName);
        TreeAsSplitsDataSet theData = null;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        
        
        // Read in an initial tree
        Tree initialTree = null;
        TreeAsSplits initialTreeAsSplits = null;
        try {
            initialTree = new Tree(new File(initialTreeFileName));
            initialTree.removeDegreeTwoVertices();
            initialTreeAsSplits= new TreeAsSplits(initialTree);          
        }
        catch (AlgorithmException anError) {
            System.out.println("Bad initial tree. "+anError.getMessage());
            System.exit(1);
        }

       //System.out.println(RefDistParams[0]+" "+RefDistParams[1]);
       LogNormalDistribution t0Dist = new LogNormalDistribution(RefDistParams[0],RefDistParams[1]);
       
       //set up the fixed values of t_0 at which to estimate the proportion of simple proposals
       for(int j=1; j <=numDisps;j++){
           disp[j-1]= t0Dist.quantile((j)/(double)(numDisps+1));
           refDensities[j-1] = t0Dist.logpdf(disp[j-1]);
       }
       double theOutput = propSimple(disp,refDensities,numDisps,numProps,initialTreeAsSplits,theData,m,outFileName);
       return theOutput;
       } 
    
    //class to implement the quadrature over t0 and save to the output file
    private static double propSimple(double[] disp,double refDens[], int numDisps,int numProps, TreeAsSplits startTree, TreeAsSplitsDataSet Data, int m ,String Filename) throws AlgorithmException{
        double[] outProps = new double[numDisps];
        double[] outProp = new double[numProps];
        File outFile = new File(Filename);
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(outFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
            
}

        
        int n= Data.numTrees;
        //System.out.println(n);
        for(int i=0; i<numDisps;i++){
            outProp = propSimpleSingleDisp(numProps,startTree, Data,m,disp[i]);
            outProps[i]=0;
            for(int j=1;j<n;j++){
                outProps[i]+=outProp[j];//sum up the log proportions
            }
            //System.out.println(Math.exp(outProps[i]));
        }
        
        double output = 0;
        double indOutput = 0;
        for (int i=1 ;i<numDisps;i++){
            //implement the quadrature
            indOutput = (Math.exp(outProps[i]+refDens[i])+Math.exp(outProps[i-1]+refDens[i-1]))*(disp[i]-disp[i-1])/2;
        output += indOutput;
        out.println(indOutput);
    }
        out.close();
        return output;
    }
    
    //class to estimate the proportion of simple bridges at each value of t0    
     private static double[] propSimpleSingleDisp(int numProp,TreeAsSplits startTree, TreeAsSplitsDataSet Data, int m,double t0) throws AlgorithmException{
        ArrayList<ForwardStepBridge> theBridges = new ArrayList();
        int n = Data.numTrees;
        double[] outputs = new double[n];
        int nPrime =startTree.getNumTaxa()-3;
        for(int i=0; i<n;i++){
            ForwardStepBridge theBridge = new BridgeWithApproxMVNLike(startTree, Data.theTrees.get(i), m);
            theBridges.add(theBridge);
            }
        for(int j=0; j<numProp;j++){
            //simulate proposals for each bridge
            for(int i=0; i<n;i++){ 
            try {
            theBridges.get(i).MargLikeIndependenceProposalPropSimple(t0);//efficient version of the proposal does not calc log likes
            outputs[i]=outputs[i]+1;
            
            }
            catch(treebase.AlgorithmException anErr){  
                
            }
           
            
            }
            
        }
        for(int i=0;i<n;i++){
            outputs[i] =Math.log(outputs[i]/(double)numProp);
        }
        
        
        return(outputs);
        
        
        
    }
     


}
