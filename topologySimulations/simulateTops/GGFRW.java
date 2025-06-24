/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package simulateTops;

import cern.jet.random.tdouble.DoubleUniform;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import static randomwalks.RandomWalkSimulatorGGF.sampleRandomWalkGGF;
import simulation.NormalDistribution;
import simulation.Random;
import treebase.Tree;
import treebase.TreeAsSplits;
import MCMC.PosteriorAnalysis;
import geodesics.Geodesic;
import java.io.IOException;
import simulation.RandomTreeSampler;
import treebase.AlgorithmException;

/**
Main class to simulate a source tree and random walks started at the simulated source tree under the GGF model
 */
public class GGFRW {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, AlgorithmException {
    int noOfTaxa =10;
    Random.setEngine(7273);    
    double t_0 = 0.25;
    int m= 20;
    int numSamples=10;
    
    //simulate the source tree
    RandomTreeSampler theSampler = new RandomTreeSampler();
    Tree x_0 = theSampler.sampleCoalescent(noOfTaxa, 2.0, 0.05);
    x_0.removeDegreeTwoVertices();
    
    
    //save the source tree
    String sourceTreeFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/ProportionSimpleSims/10Taxa20240110/SourceTree1Check.txt";
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
             
             outtt.println(x_0.toString());
             outtt.close();
        
     Random.setEngine(5555);

     //sample the trees under the random walk model
     NormalDistribution norm = new NormalDistribution(0,Math.sqrt(t_0/(double) m));
     DoubleUniform unif = new DoubleUniform(Random.getEngine());
     String outputFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/ProportionSimpleSims/10Taxa20240110/10DataPointsCheck.txt";
     TreeAsSplits x0 = new TreeAsSplits(x_0);
     TreeAsSplits[] treeSample = sampleRandomWalkGGF(x0,m,numSamples,norm,unif);
     
     //save the simulated trees
     File outputFile = new File(outputFileName);
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(outputFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outputFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
            
}

//print some summaries of the simulated data set
for(int i=0 ; i<treeSample.length ; i++){
out.println(treeSample[i]);
Geodesic g = new Geodesic(x0,treeSample[i]);
System.out.println(g.getInternalLength());
    } 
 out.close();  
  PosteriorAnalysis.countTopologiesToScreen(new File(outputFileName));
 }
    }
    

