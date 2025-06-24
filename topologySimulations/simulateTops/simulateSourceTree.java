package simulateTops;
import MCMC.PosteriorAnalysis;
import java.util.ArrayList;
import diffbase.Simulator;
import treebase.Tree;
import treebase.TreeWithTopologicalOperations;
import java.util.Iterator;
import simulation.NormalDistribution;
import cern.jet.random.tdouble.DoubleUniform;
import geodesics.Geodesic;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import randomwalks.RandomWalkSimulator;
import randomwalks.RandomWalkSimulatorGGF;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;
import topologies.LikelihoodForTopologies;
import simulation.RandomTreeSampler;
import treebase.Graph;
import treebase.Split;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

/**

 */
public class simulateSourceTree {
    
    public static void main(String[] args) throws AlgorithmException, IOException{
    /*main class to sample source tree from the sampleCoalescent- then need to check its somehow sensible.
        and then simulate some random walks started at the source tree
    */
    
    int noOfTaxa =5;
    Random.setEngine(38493);
    RandomTreeSampler theSampler = new RandomTreeSampler();
    Tree x_0 = theSampler.sampleCoalescent(noOfTaxa, 2.0, 0.05);
    x_0.removeDegreeTwoVertices();
    
    
    
    String sourceTreeFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/t0Dists/5Taxa/SourceTree1.txt";
    File sourceTreeFile= new File (sourceTreeFileName);
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
    
    //resample edge lengths from N(0,0.51)
    /*
    TreeAsSplits x_0_AsSplits = new TreeAsSplits(x_0);
    HashSet<Split> theSplits=x_0_AsSplits.getNonTrivialSplits();
    NormalDistribution norm = new NormalDistribution(0.0, 0.51);
    
    double len = 0;
    
    for(Split split : theSplits){
        x_0_AsSplits.setSplitLength(split, Math.abs(norm.sample()));
    }
    
    x_0=x_0_AsSplits.getTree();
    */
/*java.io.File theFile = new java.io.File(outputFileName);
        java.io.PrintWriter out = new java.io.PrintWriter(new java.io.BufferedWriter(new java.io.FileWriter(theFile)));    
*/
outtt.println(x_0);

 outtt.close();  
 
 double t_0 = 0.25;
 int m = 1700;
 int numToSim = 20;
 Random.setEngine(1245);
 DoubleUniform unif = new DoubleUniform(Random.getEngine());
 NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(t_0/m));
 TreeAsSplits TreeAsSplitsx_0 = new TreeAsSplits(x_0); 
 
 TreeAsSplits[] randomWalkTrees = Simulator.sampleRandomWalk(TreeAsSplitsx_0, t_0,m, numToSim,false,true);

 String outputFileName2 = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/t0Dists/5Taxa/SourceTree1_t0_0.25_20datapt.txt";
 
  File outputFile2 = new File(outputFileName2);
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

 for(TreeAsSplits tree: randomWalkTrees){
     
     outt.println(tree.toString());
     /*
     Geodesic g = new Geodesic(new TreeAsSplits(x_0),randomWalkTrees[i]);
     if((g.calcConePathDistance()-g.getLength())<Math.pow(10, -10)){
 System.out.println(i);
     }
     i=i+1;
*/
 }

 
 outt.close();
 
 Geodesic g = new Geodesic(new TreeAsSplits(x_0),randomWalkTrees[0]);
 System.out.println(g.calcConePathDistance()-g.getLength());
 
}
}

