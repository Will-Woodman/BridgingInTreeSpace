package InitialValueEstimation;
import diffbase.Simulator;
import treebase.Tree;
import simulation.NormalDistribution;
import cern.jet.random.tdouble.DoubleUniform;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;
import simulation.RandomTreeSampler;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

/**
 File to simulate a number of data sets at a time which we use to test potential methods of choosing initial values
 * for the MCMC.
 */
public class simulateDataSets {
    
    public static void main(String[] args) throws AlgorithmException, IOException{
   
    int numDataSets=100;
    
    int noOfTaxa =20;
    Random.setEngine(38493);
    String sourceTreeFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/initialValueEstimation/20Taxa/SourceTreept5_";
    
    double t_0 = 0.5;
    int m = 1700;
    int numToSim = 100;

    //sample a source tree -- here we use the same source tree for each data set but can move inside
    //the for loop if needed
    RandomTreeSampler theSampler = new RandomTreeSampler();
    Tree x_0 = theSampler.sampleCoalescent(noOfTaxa, 2.0, 0.05);
    x_0.removeDegreeTwoVertices();
    String outputFileName2 = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/initialValueEstimation/20Taxa/DataSetpt5_";
    
    //simulate the given number of different data sets
    for(int i=0; i<numDataSets;i++){
    
    //save the source tree
    File sourceTreeFile= new File (sourceTreeFileName+i+".txt");
    PrintWriter outtt;
             try {
                 FileWriter out1= new FileWriter(sourceTreeFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 outtt = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+sourceTreeFile.getName()+" failed: "
                         +"writing output to console instead.");
                 outtt = new PrintWriter(System.out);
            
}
             
 outtt.println(x_0);

 outtt.close();  
 
 DoubleUniform unif = new DoubleUniform(Random.getEngine());
 NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(t_0/m));
 TreeAsSplits TreeAsSplitsx_0 = new TreeAsSplits(x_0); 
 
 //simulate the random walks
 TreeAsSplits[] randomWalkTrees = Simulator.sampleRandomWalk(TreeAsSplitsx_0, t_0,m, numToSim,false,true);
 
  File outputFile2 = new File(outputFileName2+i+".txt");
    PrintWriter outt;
             try {
                 FileWriter outt1= new FileWriter(outputFile2, false);
                 BufferedWriter outt2 =new BufferedWriter(outt1);
                 outt = new PrintWriter(outt2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outputFile2.getName()+" failed: "
                         +"writing output to console instead.");
                 outt = new PrintWriter(System.out);
            
}

 //save the endpoints of the random walks            
 for(TreeAsSplits tree: randomWalkTrees){
     
     outt.println(tree.toString());
     

 }
 outt.close();
 
    }
    }
}

 


