/*
RandomWalkSimulatorGGF
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

package randomwalks;

/**
 * Set of static methods for performing simulations of diffusions -- WW modifying to GGF RW
 * 
 */

import simulation.NormalDistribution;
import simulation.GammaDistribution;
import cern.jet.random.tdouble.DoubleUniform;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Graph;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;
import treebase.TreeWithTopologicalOperations;
import treedatasets.OperationsOnTreeAsSplits;
//import Source Packages.MCMC.PosteriorAnalysis;


//import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class RandomWalkSimulatorGGF {

    /* Methods based on maps of splits -------------------------------------- */

    /** Diffuse a tree (represented as TreeAsSplits). Internal edges only.
      Do numSteps replacements FOR EACH EDGE.
      NB 1: tree must be fully resolved!
      NB 2: much faster to diffuse a true Tree object rather than TreeAsSplits
     */
    public static void randomWalkTreeGGF(TreeAsSplits theTree, int numSteps,NormalDistribution norm, DoubleUniform unif) throws AlgorithmException {
        /* Main loop */
        for (int i=0; i<numSteps; i++) {
           sampleDirectMVN(theTree,norm,unif);
        } // E
    }

    /** Generate many realizations from a single tree */
    public static TreeAsSplits[] sampleRandomWalkGGF(TreeAsSplits theTree, int numSteps, int numSamples, NormalDistribution norm, DoubleUniform unif) throws AlgorithmException {
        TreeAsSplits[] res = new TreeAsSplits[numSamples];
        for (int i=0; i<numSamples; i++) {
            TreeAsSplits t = theTree.clone();
            RandomWalkSimulatorGGF.randomWalkTreeGGF(t, numSteps, norm, unif);
            res[i] = t;
        }
        return res;
    }


//method below is modified from the one in ForwardStepBridge but modifies the treeAsSplits instead of outputting a stepResult
 static public void sampleDirectMVN(TreeAsSplits theTree, NormalDistribution normDist, DoubleUniform unifDist) throws AlgorithmException {
     
        int nprime = theTree.getNumTaxa()-3;

        /* Make an array of splits */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(theTree.getNonTrivialSplits()); 
        /* SORT array list if necessary */
//        Collections.sort(splits); // Sort to ensure repeatability with same random seed       
        
        /* Make a uniform direction */
        double[] u = new double[nprime];
        double lenSqu = 0.0;
        for (int i=0; i<nprime; i++) {
            u[i] = normDist.sample();
            lenSqu += u[i]*u[i];
        }
        /* Normalize and compute edge lengths */
        double[] l = new double[nprime];
        double len = Math.sqrt(lenSqu);
        

        ArrayList<IntDoublePair> boundaries = new ArrayList();
        for (int i=0; i<nprime; i++) {
            Split p = splits.get(i);
            double y = theTree.getSplitLength(p);
            l[i] = y+u[i];
            try {
                // Set edge lengths
                theTree.setSplitLength(p, Math.abs(l[i]));
            } catch (AlgorithmError ex) {
                System.out.println("Error in sample direct MVN in Bridge algorithm: missing split.");
            }
            if ((l[i]<0.0)&&(u[i]<0.0)) {
                boundaries.add(new IntDoublePair(i,-(y*len)/u[i])); // -y*len/u[i] is the "time" when you hit the boundary. Order by time!
            } 
        }
                
        if (boundaries.size()!=0) {
            /* We crossed at least one boundary */
            Collections.sort(boundaries);
            for (int i=0; i<boundaries.size(); i++) {
                /* Do the nni's in turn  */
                int ind = boundaries.get(i).theInt;
                Split p = splits.get(ind);
                // Sanity check: does theTree contain p?
                if (!theTree.contains(p)) {
                    throw new AlgorithmError("Missing p in direct MVN sampler.");
                }
                Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(theTree, p);
                Split q = nniSplit[unifDist.nextIntFromTo(0, 1)];
                theTree.remove(p);
                theTree.add(q, Math.abs(l[ind]));    // UNCOMMENT AND REPLACE LINE BELOW AFTER DEBUG
//            theTree.addWithCheck(q, Math.abs(l[ind]));
            }
        }
                
    }
 
 static public class IntDoublePair implements Comparable {
        public int theInt;
        public double theDouble;
        
        public IntDoublePair(int i, double d) {
            theInt = i;
            theDouble = d;
        }

        public int compareTo(Object t) {
            IntDoublePair arg = (IntDoublePair) t;
            if (theDouble<arg.theDouble) return -1;
            if (theDouble>arg.theDouble) return 1;
            return 0;
        }        
    }
 
 public static void main(String[] args) throws AlgorithmException{
     double Disp = 0.0189589;
     int m= 50;
     int numSamples=106;
     NormalDistribution norm = new NormalDistribution(0,Math.sqrt(Disp/(double) m));
     Random.setEngine(1000);
     DoubleUniform unif = new DoubleUniform(Random.getEngine());
     String outputFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/OldRokasData/InitialParamsForwardSim/m50_n106_seed_1000.txt";
     
     TreeAsSplits x0 = new TreeAsSplits(new Tree(new File("/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/OldRokasData/InitialEstTree.txt")));
     TreeAsSplits[] treeSample = sampleRandomWalkGGF(x0,m,numSamples,norm,unif);
     
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
/*java.io.File theFile = new java.io.File(outputFileName);
        java.io.PrintWriter out = new java.io.PrintWriter(new java.io.BufferedWriter(new java.io.FileWriter(theFile)));    
*/
for(int i=0 ; i<treeSample.length ; i++){
out.println(treeSample[i]);
        
    } 
 out.close();  
  //PosteriorAnalysis.countTopologiesToScreen(new File(outputFileName));
 }
}