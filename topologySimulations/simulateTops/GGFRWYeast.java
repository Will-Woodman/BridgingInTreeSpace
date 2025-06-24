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
import geodesics.Geodesic;
import java.io.IOException;
import treebase.AlgorithmException;

/**
MMain class to simulate random walks started at a given source tree under the GGF model
 */
public class GGFRWYeast {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, AlgorithmException {
    String sourceTreeFileName = args[0];
    String outputFileName = args[1];
    int seed = Integer.parseInt(args[2]);  
    Double t_0 = Double.parseDouble(args[3]);
    Integer m= Integer.parseInt(args[4]);
    Integer numSamples=Integer.parseInt(args[5]);     
   
    Random.setEngine(seed);    

    
    //save the source tree
    File sourceTreeFile= new File (sourceTreeFileName);
    Tree x_0 = new Tree(sourceTreeFile);
    TreeAsSplits x0 = new TreeAsSplits( x_0);

     //sample the trees under the random walk model
     NormalDistribution norm = new NormalDistribution(0,Math.sqrt(t_0/(double) m));
     DoubleUniform unif = new DoubleUniform(Random.getEngine());
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

    //print the geodesic distances to an output file
    for(int i=0 ; i<treeSample.length ; i++){
        out.println(treeSample[i]);
        Geodesic g = new Geodesic(x0,treeSample[i]);
        System.out.println(g.getInternalLength());
        } 
     out.close();  
     }
    }
    

