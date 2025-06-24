/*
   Initializer
    Copyright (C) 2015  Tom M. W. Nye

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

    Contact the author at:  <tom.nye@ncl.ac.uk>
                            <http://www.mas.ncl.ac.uk/~ntmwn/>
 */

package bridge;

import diffbase.Simulator;
import geodesics.Geodesic;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
    Class used to find initial x_0 and t_0 (source and dispersion params) for an
    input data set of trees. 
 */

public class Initializer {
    
    /** Sturm / Bacak algorithm to find Frechet mean */
    static public TreeAsSplits computeApproxFrechetMean(ArrayList<TreeAsSplits> theTrees, int numIts) {
        cern.jet.random.tdouble.DoubleUniform unif = new cern.jet.random.tdouble.DoubleUniform(simulation.Random.getEngine());
        int n = theTrees.size();
        TreeAsSplits point = (TreeAsSplits) simulation.CategoricalDistribution.sampleFromListWithReplacement(theTrees, unif);
        try {
            for (int k=1; k<numIts; k++) {
                TreeAsSplits nextTree = (TreeAsSplits) simulation.CategoricalDistribution.sampleFromListWithReplacement(theTrees, unif);
                Geodesic h = new Geodesic(point, nextTree);
                TreeAsSplits newTree = h.getTree(1.0/((double)(k+1)));
                point = newTree;
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error with iterative method to find approx Frechet mean. "+anErr.getMessage());
        }
        return point;
    }
    
    
    /** Extract distances from x0 */
    static public double[] calcSquDistancesFromx0(TreeAsSplits x0, ArrayList<TreeAsSplits> theTrees) {
        double[] d = new double[theTrees.size()];
        for (int i=0; i<theTrees.size(); i++) {
            try {
                Geodesic g = new Geodesic(x0, theTrees.get(i));
                double tmp = g.getInternalLength();
                d[i] = tmp*tmp;
            } catch (AlgorithmError ex) {
                System.out.println("Error making geodesic.");
            }
        }
        return d;
    }
    
    /** Extract distances from x0 */
    static public double[] calcSquDistancesFromx0(TreeAsSplits x0, TreeAsSplits[] theTrees) {
        double[] d = new double[theTrees.length];
        for (int i=0; i<theTrees.length; i++) {
            try {
                Geodesic g = new Geodesic(x0, theTrees[i]);
                double tmp = g.getInternalLength();
                d[i] = tmp*tmp;
            } catch (AlgorithmError ex) {
                System.out.println("Error making geodesic.");
            }
        }
        return d;
    }
    
    
    private static double[] calcMeanAndVar(double[] x) {
        double s=0, ss=0;
        for (int i=0; i<x.length; i++) {
            s += x[i];
            ss += x[i]*x[i];
        }
        double mu = s/x.length;
        double var = ss/x.length-mu*mu;
        return new double[] {mu, var};
    }
    
    
    /** Simulate a bunch of points for given x0, t0 */
    static public double[] simulateDistances(TreeAsSplits x0, double t0, int numSteps, int numParticles) {
        
        /* Get RW trees */
        TreeAsSplits[] particles = Simulator.sampleRandomWalk(x0, t0, numSteps, numParticles, true, false);
        
        /* Compute distances */
        return calcSquDistancesFromx0(x0, particles);
        
    }
    
    
    /** Initialize! */
    public static void main(String[] args) {

        String dataFilename = "/home/ntmwn/research/diffusion/bridge/sims/sample5_0.1.txt"; 
//        String dataFilename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.1.txt"; 
        int numSteps = 50;
        int numItsFrechet = 100;
        int numParticles = 100;
        double testsd = 0.09;

        Random.setEngine((int)System.currentTimeMillis());
        
        /* End inputs ------------------------------------------------------- */
        
        /* Read in the data */
        File inFile = new File(dataFilename);
        TreeAsSplitsDataSet theData = null;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        
//        ArrayList<Split> splits = new ArrayList();
//        Split[][] topologies = new Split[15][2];
//        int topCount = 0;
//        splits.addAll(theData.allNonTrivialSplits);
//        Collections.sort(splits);
//        for (int i=0; i<splits.size(); i++) {
//            Split p = splits.get(i);
//            for (int j=i+1; j<splits.size(); j++) {
//                Split q = splits.get(j);
//                if (p.isCompatible(q)) {
//                    topologies[topCount][0] = p;
//                    topologies[topCount][1] = q;
//                    topCount++;
//                }
//            }
//        }
//        
//        for (int i=0; i<theData.numTrees; i++) {
//            TreeAsSplits theTree = theData.getTree(i);
//            ArrayList<Split> s = new ArrayList();
//            s.addAll(theTree.getNonTrivialSplits());
//            Collections.sort(s);
//            for (Split p : s) {
//                System.out.print(theTree.getSplitLength(p)+" ");
//            }
//            // Find the topology
//            for (int j=0; j<topologies.length; j++) {
//                if ((topologies[j][0].equals(s.get(0)))&&(topologies[j][1].equals(s.get(1)))) {
//                    System.out.print(j+1);
//                }
//            }
//            System.out.println("");
//        }
        
        
        TreeAsSplits x0 = computeApproxFrechetMean(theData.theTrees,numItsFrechet);
        System.out.println("Estimated x0 = ");
        System.out.println(x0.toString());
        
        
        double[] distSqu = calcSquDistancesFromx0(x0, theData.theTrees);
        
        System.out.println("\nSquared distances:-");
        for (int i=0; i<distSqu.length; i++) {
            System.out.println(String.format("%7.7f", distSqu[i]));
        }
        
        double[] stats = calcMeanAndVar(distSqu);
        System.out.println("\nMean = "+String.format("%7.7f", stats[0]));
        System.out.println("Var  = "+String.format("%7.7f", stats[1]));
        
        distSqu = simulateDistances(x0, testsd*testsd, numSteps, numParticles);
        
        System.out.println("\nSimulated squared distances:-");
        for (int i=0; i<distSqu.length; i++) {
            System.out.println(String.format("%7.7f", distSqu[i]));
        }

        
        stats = calcMeanAndVar(distSqu);
        System.out.println("\nSimulated mean = "+String.format("%7.7f", stats[0]));
        System.out.println("Simulated var  = "+String.format("%7.7f", stats[1]));
    }
        

}
