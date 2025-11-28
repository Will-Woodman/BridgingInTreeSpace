/*
ProportionOfSimpleProposalsDisp
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
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import simulation.LogNormalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
 Calculate the normalising constant of the reference distribution for the marginal likelihood estimators
 * that include t0(See paper)
 */
public class ProportionOfSimpleProposalsDisp {
    
    //class to implement the quadrature over t0
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
        System.out.println(n);
        for(int i=0; i<numDisps;i++){
            
            outProp = propSimpleSingleDisp(numProps,startTree, Data,m,disp[i]);
            outProps[i]=0;
            for(int j=1;j<n;j++){
                outProps[i]+=outProp[j];
            }
            System.out.println(Math.exp(outProps[i]));
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
            for(int i=0; i<n;i++){ 
            try {
            theBridges.get(i).MargLikeIndependenceProposalPropSimple(t0);
            outputs[i]=outputs[i]+1;//if bridge built successfully, add one to the total
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
     
   
      
      //main class to test the method
      public static void main(String[] args) throws AlgorithmException, IOException {
        //number of times to repeat to see the variance in the estimate
        int numFiles =8;
        //number of proposals for each bridge and each t0
        int numProps =1000;

        //the number of values of t0 to estimate at: (could make this parallel to estimate different values of t0 on different cores
        int numDisps =9;
        double[] disp=new double[numDisps];
        double[] refDensities = new double[numDisps];
        
        
        String initialTreeFileName ="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/ProportionSimpleSims/10Taxa20240110/SourceTree1.txt";
        String dataFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/ProportionSimpleSims/10Taxa20240110/100DataPoints.txt";
        String outFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/ProportionSimpleSims/10Taxa20240110/100DataPoints_test1";
        String RefDistParamsFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/ProportionSimpleSims/10Taxa20240110/RefDistParams100DataPoint.txt";
        Double[] RefDistParams= new Double[2];
        
        BufferedReader br = new BufferedReader(new FileReader(RefDistParamsFileName));
        String s;
        //introduce t so we don't end up reading last line with the acceptance rates
        for(int i=0;i<2;i++){
        s = br.readLine();
            if(i ==1){
               String[] ss = s.split(" ");
            RefDistParams[0]=Double.parseDouble(ss[0]);
            RefDistParams[1]=Double.parseDouble(ss[1]);
            }
         
        }
        
        int seed=27771;
        Random.setEngine(seed);
          
        int m=20;
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
        
      
       //var about initial tree:
       System.out.println(RefDistParams[0]+" "+RefDistParams[1]);
       LogNormalDistribution t0Dist = new LogNormalDistribution(RefDistParams[0],RefDistParams[1]);
       
       for(int j=1; j <=numDisps;j++){
           double toPrint = (j)/(double) (numDisps+1);
           System.out.println((j)/(double)(numDisps+1));
           disp[j-1]= t0Dist.quantile((j)/(double)(numDisps+1));
           refDensities[j-1] = t0Dist.logpdf(disp[j-1]);
           System.out.println(refDensities[j-1]);
       }
        for(int j=0;j<numFiles;j++){
       System.out.println(propSimple(disp,refDensities,numDisps,numProps,initialTreeAsSplits,theData,m,outFileName+j+".txt")+" overall");  
       } 
       
       //String outFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/10Taxa/SteppingStone/IndepProps";
       
     
    }
    
     
  }
     



