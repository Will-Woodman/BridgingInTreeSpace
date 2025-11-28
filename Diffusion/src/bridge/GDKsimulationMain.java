/*
GDKSimulationMain
>>>>>>> 9326c42338de26f6a866256b3baa4f195ce299d6
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

package bridge;
import MCMC.PosteriorAnalysis;
import java.io.File;
import static bridge.TreeDistributions.sampleExponentialBumpViaMVNStep;
import cern.jet.random.tdouble.DoubleUniform;
import geodesics.Geodesic;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import treebase.TreeAsSplits;
import simulation.NormalDistribution;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Tree;

/**
 *
 * @author will
 */
public class GDKsimulationMain {

    /**
     Main file to simulate from the Gaussian kernel distribution by MCMC using the
     * sampleExponentialBumpViaMVNStep class
     */
    public static void main(String[] args) throws IOException, AlgorithmError, AlgorithmException {
       
        String initialTreeFilename = args[0];  //Tree x0 the 'centre' of the distribution
        String outputFilename = args[1]; //output file for the simulated trees
        String outputFilenameTops = args[2]; //output file for the topologies of the simulated trees
        String outputFilenameDist = args[3]; //output file for the geodesic distances of the simulated trees to the source        
        int seed = Integer.parseInt(args[4]); // see for the random engine
        double r = Double.parseDouble(args[5]); //value of t0 in the kernel
        double sd = Double.parseDouble(args[6]); // value for the proposal distribution which is a random walk proposal
        int nits = Integer.parseInt(args[7]); // number of iterations in the MCMC    
        int burnits = Integer.parseInt(args[8]); // number of burn-in iterations in the MCMC    
        int thin = Integer.parseInt(args[9]); // number of burn-in iterations in the MCMC    
        
        
        Random.setEngine(seed);
        //make the files
        File outFile = new File(outputFilename);
        File outFileTops = new File(outputFilenameTops);
        File outFileDist = new File(outputFilenameDist);
        
        
        /* Read in an initial tree */
        
        Tree initialTree = null;
        try {
            initialTree = new Tree(new File(initialTreeFilename));
            initialTree.removeDegreeTwoVertices();
        }
        catch (AlgorithmException anError) {
            System.out.println("Bad initial tree. "+anError.getMessage());
            System.exit(1);
        }
        TreeAsSplits x0 = new TreeAsSplits(initialTree);
        
        NormalDistribution normDist = new NormalDistribution(0.0, sd);
        DoubleUniform unifDist = new DoubleUniform(simulation.Random.getEngine());
        boolean[] CBound = new boolean[1];
       
       //run the sampler 
        ArrayList<TreeAsSplits> sample = sampleExponentialBumpViaMVNStep(x0, r, sd, burnits, thin, nits, normDist, unifDist,CBound);
       
        //print the tree sample
        PrintWriter out;
         if (outFile==null) {
             out = new PrintWriter(System.out);
         }
         else {
             try {
                 out = new PrintWriter(new BufferedWriter(new FileWriter(outFile, false)));
             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
             }
         }
         
         //print the sampled trees
         for (int i=1; i<sample.size(); i++) {
            TreeAsSplits nextTree = sample.get(i);
            out.println(nextTree.toString(true));
        }
         
        out.close();
        
       
        //print the distances of the sample points from x0
        PrintWriter out2;
         if (outFileDist==null) {
             out2 = new PrintWriter(System.out);
         }
         else {
             try {
                 out2 = new PrintWriter(new BufferedWriter(new FileWriter(outFileDist, false)));
             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out2 = new PrintWriter(System.out);
             }
         }
         
        //print out the distances between the source parameter x0 and the sampled points
        for (int i=1; i<sample.size(); i++) {
            Geodesic h = new Geodesic(sample.get(i), x0);
            out2.println(String.format("%7.7f", h.getInternalLength()));
        }
        
        out2.close();

        //print out the frequencies of the topologies in sample
       PosteriorAnalysis.countTopologiesToFile(sample, outFileTops);
      
        
    }
       
    
}
