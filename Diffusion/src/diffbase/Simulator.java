/*
 * Simulator.java

    Copyright (C) 2013  Tom M. W. Nye

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

package diffbase;

/**
 * Set of static methods for performing simulations of diffusions
 */

import java.io.IOException;
import treedatasets.OperationsOnTreeAsSplits;
import geodesics.Geodesic;
import simulation.NormalDistribution;
import simulation.GammaDistribution;
import cern.jet.random.tdouble.DoubleUniform;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Graph;
import treebase.Graph.Edge;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;
import treebase.TreeWithTopologicalOperations;
//import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class Simulator {

    /* Methods based on maps of splits -------------------------------------- */
    
    private static Split singleStep(TreeAsSplits theTree, Split p, double x0, double sd, double proposed, DoubleUniform unif, boolean useExactDistrib) throws AlgorithmError {
        Split q = null;
        if (useExactDistrib) {
            
            // Use solution to heat equation on 3 spider
            if (proposed<0.0) {
                int topInd = unif.nextIntFromTo(1, 3);
                if (topInd==1) {
                    theTree.setSplitLength(p, -proposed);
                }
                else {
                    Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(theTree, p);
                    q = nniSplit[topInd-2]; // topInd-2 = 0 or 1
                    theTree.remove(p);
                    theTree.add(q, -proposed);
                }
            }
            else {
                // Positive edge length
                double u = unif.nextDouble();
                if (u*NormalDistribution.pdf(proposed, x0, sd)<NormalDistribution.pdf(-proposed, x0, sd)) {
                    // Cross boundary
                    int topInd = unif.nextIntFromTo(1, 3);
                    if (topInd==1) {
                        theTree.setSplitLength(p, proposed);
                    }
                    else {
                        Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(theTree, p);
                        q = nniSplit[topInd-2]; // topInd-2 = 0 or 1
                        theTree.remove(p);
                        theTree.add(q, proposed);
                    }
                }
                else {
                    theTree.setSplitLength(p, proposed);
                }
            }
            
        }
        else {
           
            // Use standard normal innovation
            if (proposed>0.0) {
                theTree.setSplitLength(p, proposed);
            }
            else {
                int topInd = unif.nextIntFromTo(1, 3);
                if (topInd==1) {
                    theTree.setSplitLength(p, -proposed);
                }
                else {
                    Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(theTree, p);
                    q = nniSplit[topInd-2]; // topInd-2 = 0 or 1
                    theTree.remove(p);
                    theTree.add(q, -proposed);
                }
            }
            
        }
        return q;
    } 

    /** Diffuse a tree (represented as TreeAsSplits). Internal edges only.
      Do numSteps replacements FOR EACH EDGE.
      NB 1: tree must be fully resolved!
      NB 2: much faster to diffuse a true Tree object rather than TreeAsSplits
     */
    public static void randomWalkTree(TreeAsSplits theTree, double var, int numSteps, boolean randomSweep, boolean useExactDistrib) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        Simulator.randomWalkTree(theTree, numSteps, randomSweep, useExactDistrib, norm, unif);
    }
    public static void randomWalkTree(TreeAsSplits theTree, int numSteps, boolean randomSweep, boolean useExactDistrib, NormalDistribution norm, DoubleUniform unif) {

        /** If you use this method make sure the normal distribution has variance = total var / numSteps */
        double sd = norm.getSigma();

        /* First get the non-trivial splits in a fixed order */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(theTree.getNonTrivialSplits());
        Collections.sort(splits);
        if (splits.size()<theTree.getNumTaxa()-3) {
            System.out.println("Warning: request made to random walk an unresolved tree. Trees must be randomly resolved first.");
        }

        /* Main loop */
        Split q;
        double l;
        int ind;
        for (int i=0; i<numSteps; i++) {
            for (int j=0; j<splits.size(); j++) {

                if (randomSweep)
                    ind = unif.nextIntFromTo(0, splits.size()-1);
                else
                    ind = j;

                Split p = splits.get(ind);
                double x0 = theTree.getSplitLength(p);
                l = x0+norm.sample();

                try {
                    q = singleStep(theTree, p, x0, sd, l, unif, useExactDistrib);
                    if (q!=null) {
                        splits.set(ind, q);
                    }
                    
                }
                catch (AlgorithmException err) {
                    System.out.println("Error setting edge length in a diffusion.");
                }
            } // End loop thru' splits
        } // End main loop

    }

    /** Generate one realization from a single tree.
     numSteps steps PER EDGE each with variance var/numSteps */
    public static TreeAsSplits sampleRandomWalk(TreeAsSplits theTree, double var, int numSteps, boolean randomSweep, boolean useExactDistrib) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        return Simulator.sampleRandomWalk(theTree, numSteps, randomSweep, useExactDistrib, norm, unif);
    }
    public static TreeAsSplits sampleRandomWalk(TreeAsSplits theTree, int numSteps, boolean randomSweep, boolean useExactDistrib, NormalDistribution norm, DoubleUniform unif) {
        /** If you use this method make sure the normal distribution has variance = total var / numSteps */
        TreeAsSplits t = theTree.clone();
        Simulator.randomWalkTree(t, numSteps, randomSweep, useExactDistrib, norm, unif);
        return t;
    }

    /** Generate many realizations from a single tree */
    public static TreeAsSplits[] sampleRandomWalk(TreeAsSplits theTree, double var, int numSteps, int numSamples, boolean randomSweep, boolean useExactDistrib) {
        TreeAsSplits[] res = new TreeAsSplits[numSamples];
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        for (int i=0; i<numSamples; i++) {
            TreeAsSplits t = theTree.clone();
            Simulator.randomWalkTree(t, numSteps, randomSweep, useExactDistrib, norm, unif);
            res[i] = t;
        }
        return res;
    }

    /* Methods based on trees ----------------------------------------------- */
    
    
    private static Edge singleStep(TreeWithTopologicalOperations theTree, Edge p, double x0, double sd, double proposed, DoubleUniform unif, boolean useExactDistrib) throws AlgorithmException {
    
        Edge q = null;
        if (useExactDistrib) {
            
            // Use solution to heat equation on 3 spider
            if (proposed<0.0) {
                int topInd = unif.nextIntFromTo(1, 3);
                if (topInd==1) {
                    p.setLength(-proposed);
                }
                else {
                    boolean choice = (topInd==2); // choice=true if topInd==2, choice=false if topInd==3
                    q = theTree.performNNI(p, choice, -proposed);
                }
            }
            else {
                // Positive edge length
                double u = unif.nextDouble();
                if (u*NormalDistribution.pdf(proposed, x0, sd)<NormalDistribution.pdf(-proposed, x0, sd)) {
                    // Cross boundary
                    int topInd = unif.nextIntFromTo(1, 3);
                    if (topInd==1) {
                        p.setLength(proposed);
                    }
                    else {
                        boolean choice = (topInd==2); // choice=true if topInd==2, choice=false if topInd==3
                        q = theTree.performNNI(p, choice, proposed);
                    }
                }
                else {
                    p.setLength(proposed);
                }
            }
            
        }
        else {
           
            // Use standard normal innovation
            if (proposed>0.0) {
                p.setLength(proposed);
            }
            else {
                int topInd = unif.nextIntFromTo(1, 3);
                if (topInd==1) {
                    p.setLength(-proposed);
                }
                else {
                    boolean choice = (topInd==2); // choice=true if topInd==2, choice=false if topInd==3
                    q = theTree.performNNI(p, choice, -proposed);
                }
            }
            
        }
        return q;
    } 


    /** Randomly walk a tree. Internal edges only.
      Do numSteps replacements ON EACH EDGE with variance var/numSteps.
      NB: tree must be fully resolved!  */
    public static void randomWalkTree(TreeWithTopologicalOperations theTree, double var, int numSteps, boolean randomSweep, boolean useExactDistrib) {
        randomWalkTree(theTree, numSteps, randomSweep, useExactDistrib, new NormalDistribution(0.0, Math.sqrt(var/numSteps)), new DoubleUniform(simulation.Random.getEngine()));
    }
    public static void randomWalkTree(TreeWithTopologicalOperations theTree, int numSteps, boolean randomSweep, boolean useExactDistrib, NormalDistribution norm, DoubleUniform unif) {

        /* First get the internal edges in a fixed order */
        ArrayList<Graph.Edge> edges = new ArrayList();
        Iterator<Graph.Edge> itE = theTree.getEdgeIterator();
        for (int i=0; i<theTree.numEdges(); i++) {
            Graph.Edge e = itE.next();
            if (!Tree.isEdgeTerminal(e)) edges.add(e);
        }

        if (edges.size()<theTree.numTaxa()-3) {
            System.out.println("Warning: request made to diffuse an unresolved tree. Trees must be randomly resolved first.");
        }

        /* Main loop */

        Edge f;
        double l, x0;
        int ind;
        double sd = norm.getSigma();
        
        for (int i=0; i<numSteps; i++) {
            for (int j=0; j<edges.size(); j++) {

                if (randomSweep)
                    ind = unif.nextIntFromTo(0, edges.size()-1);
                else
                    ind = j;

                Edge e = edges.get(ind);
                x0 = e.getLength();
                l = norm.sample()+x0;

                try {
                    f = singleStep(theTree, e, x0, sd, l, unif, useExactDistrib);
                    if (f!=null) {
                        edges.set(ind, f);
                    }
                }
                catch (AlgorithmException err) {
                    System.out.println("Error setting edge length in a diffusion."+err.getMessage());
                }
            } // End loop thru' splits
        } // End main loop

    }

    /** Generate one realization from a single tree */
    public static Tree sampleRandomWalk(Tree theTree, double var, int numSteps, boolean randomSweep, boolean useExactDistrib) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        return Simulator.sampleRandomWalk(theTree, numSteps, randomSweep, useExactDistrib, norm, unif);
    }
    public static Tree sampleRandomWalk(Tree theTree, int numSteps, boolean randomSweep, boolean useExactDistrib, NormalDistribution norm, DoubleUniform unif) {
        TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(theTree);
        randomWalkTree(t, numSteps, randomSweep, useExactDistrib, norm, unif);
        return t;
    }
    /** Generate many realizations from a single tree */
    public static Tree[] sampleRandomWalk(Tree theTree, double var, int numSteps, boolean randomSweep, boolean useExactDistrib, int numSamples) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        Tree[] res = new Tree[numSamples];
        for (int i=0; i<numSamples; i++) {
            TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(theTree);
            randomWalkTree(t, numSteps, randomSweep, useExactDistrib, norm, unif);
            res[i] = t;
        }
        return res;
    }
    
    
    

    /* Parallel method to sample large numbers of trees */
    public static HashSet<Tree> sampleRandomWalkParallel(Tree theTree, double var, int numSteps, boolean randomSweep, boolean useExactDistrib, int numSamples, int numProcessors)  {
        int chunk = Math.round(((float)numSamples)/((float)numProcessors));
        int start=0, end=chunk, seed;
        HashSet<Tree> res = new HashSet();

        ParallelRWChunk[] processes = new ParallelRWChunk[numProcessors];
        for (int i=0; i<numProcessors; i++) {
            // Create processes
            seed = simulation.Random.getEngine().nextInt();
            if (i==(numProcessors-1)) end=numSamples;
            ParallelRWChunk p = new ParallelRWChunk(theTree, numSamples, var, numSteps, randomSweep, useExactDistrib, seed);
            p.start();
            processes[i] = p;
            start = end;
            end = start+chunk;
        }
        // Wait
        try {
             for (int i=0; i<numProcessors; i++) {
                processes[i].join();
                res.addAll(processes[i].getTrees());
            }
        }
        catch (InterruptedException ex) {
           System.out.println("Error: something funny happened with a thread interrupt while projecting trees.");
        }
        return res;

    }

    public static class ParallelRWChunk extends Thread {

        private Tree theTree;
        private int numSamples, numSteps;
        private HashSet<Tree> res;
        private double sigma;
        private NormalDistribution norm;
        private DoubleUniform unif;
        private boolean randomSweep;
        private boolean useExactDistrib;


        public ParallelRWChunk(Tree t, int nSamples, double v, int nSteps, boolean rs, boolean ed, int seed) {
            theTree = t;
            numSamples = nSamples;
            sigma = Math.sqrt(v/nSteps);
            numSteps = nSteps;
            res = new HashSet();
            DoubleMersenneTwister engine = new DoubleMersenneTwister(seed);
            norm = new NormalDistribution(0.0, sigma, engine);
            unif = new DoubleUniform(engine);
            randomSweep = rs;
            useExactDistrib = ed;
        }

       public void run() {
            for (int i=0; i<numSamples; i++) {
                TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(theTree);
                randomWalkTree(t, numSteps, randomSweep, useExactDistrib, norm, unif);
                res.add(t);
            }
       }
       
       public HashSet<Tree> getTrees() {
           return res;
       }
    }
//
//
//    /** Use a gamma distribution to perturb pedant edges. */
//    public static void perturbPendantEdges(TreeAsSplits theTree, double var) {
//        GammaDistribution gamma = new GammaDistribution(1.0, 1.0); // Arguments are shape and scale
//        ArrayList<Split> splits = new ArrayList();
//        splits.addAll(theTree.getSplits());
//        splits.removeAll(theTree.getNonTrivialSplits());
//        Collections.sort(splits);
//
//        Iterator<Split> it = splits.iterator();
//        double mu, x;
//        Split p;
//        try {
//            while (it.hasNext()) {
//                p = it.next();
//                mu = theTree.getSplitLength(p);
//                if (Math.abs(mu)>1.0E-10) {
//                    x = gamma.sample(mu*mu/var, var/mu); // mean = shape*scale = mu; var = shape*scale^2 = sigma^2
//                    theTree.setSplitLength(p, x);
//                }
//            }
//        }
//        catch (AlgorithmError anErr) {
//            System.out.println("Error perturbing pendant edge lengths for a diffusion. This should not be possible. ");
//        }
//    }
//
//
//    /* "FAST" SIMULATION METHODS ---------------------------------------------*/
//    
//    /* Legacy code: not currently used! */
//
//    /** Simulate a sequence of X~N(0,sigma^2) such that sum X_i = s*/
//    public static double[] simulateGaussianBridge(NormalDistribution norm, int nSteps, double s) {
//        double sig = norm.getSigma();
//
//        double[] res = new double[nSteps];
//
//        double x = 0.0;
//        int k; // Num steps remaining
//        for (int i=0; i<(nSteps-1); i++) {
//            k = nSteps-i;
//            res[i] = norm.sample((s-x)/k, sig*Math.sqrt((k-1)/((double)k)));
//            x += res[i]; // running total
//        }
//
//        res[nSteps-1] = s-x;
//        return res;
//    }
//
//    /** Simulate single step of Gaussian bridge */
//    public static double simulateGaussianBridgeStep(NormalDistribution norm, double sig, int nRemainingSteps, double target, double current) {
//        double res = norm.sample((target-current)/nRemainingSteps, sig*Math.sqrt((nRemainingSteps-1)/((double)nRemainingSteps)));
//        return res;
//    }
//
    /* Debug area ------------------------------------------------------------*/

    private static void doIO(treebase.Tree[] trees, String filename) throws IOException {
        java.io.File theFile = new java.io.File(filename);
        java.io.PrintWriter out = new java.io.PrintWriter(new java.io.BufferedWriter(new java.io.FileWriter(theFile)));
        for (int i=0; i<trees.length; i++) {
            out.println(trees[i].toString());
        }
        out.close();

        java.io.File inFile = new java.io.File(filename);
        treedatasets.TreeAsSplitsDataSet theData = new treedatasets.TreeAsSplitsDataSet(inFile);
    }
    
    private static int calcRFDistance(treebase.Tree t, HashSet<Split> x) {
        HashSet<Split> y = new HashSet();
        Iterator<treebase.Graph.Edge> it = t.getEdgeIterator();
        while (it.hasNext()) try {
            y.add(t.getSplit(it.next()));
            } catch (AlgorithmError ex) {
                Logger.getLogger(Simulator.class.getName()).log(Level.SEVERE, null, ex);
            }
        
        int rf = 0;
        Iterator<Split> itx = x.iterator();
        while (itx.hasNext()) if (!y.contains(itx.next())) rf++;
        Iterator<Split> ity = y.iterator();
        while (ity.hasNext()) if (!x.contains(ity.next())) rf++;
        return rf;
    }
    
    private static void topologicalSummary(treebase.Tree[] trees, treebase.Tree x)  {
        System.out.println("\nRF distances:-");
        HashSet<Split> xsplits = new HashSet();
        Iterator<treebase.Graph.Edge> it = x.getEdgeIterator();
        while (it.hasNext()) try {
            xsplits.add(x.getSplit(it.next()));
            } catch (AlgorithmError ex) {
                Logger.getLogger(Simulator.class.getName()).log(Level.SEVERE, null, ex);
            }
        int s=0;
        for (int i=0; i<trees.length; i++) {
            int r = calcRFDistance(trees[i],xsplits);
            //System.out.println(r);
            s += r;
        }
        double m = ((double)s)/((double)trees.length*2.0);
        System.out.println("Mean = "+String.format("%7.7f",m ));
    }


    public static void main(String[] args) throws AlgorithmException, IOException {

        Random.setEngine(7);

        int numSamples=100;
        int numStepsPerEdge = 500;
        
        treebase.Tree x0 = new treebase.Tree("((A:1,B:1):0.1,(C:1,D:1):0.1);");
        x0.removeDegreeTwoVertices();
        String filename = "/home/ntmwn/research/diffusion/bridge/sims/sample4_0.2.txt";
        double sig = 0.2;
        treebase.Tree[] trees = sampleRandomWalk(x0, sig*sig, numStepsPerEdge, false, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, x0);

        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample4_0.4.txt";
        sig = 0.4;
        trees = sampleRandomWalk(x0, sig*sig, numStepsPerEdge, false, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, x0);

        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample4_0.6.txt";
        sig = 0.6;
        trees = sampleRandomWalk(x0, sig*sig, numStepsPerEdge, false, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, x0);


/*
        treebase.Tree x0 = new treebase.Tree("((A:1,B:1):0.1,C:1,(D:1,E:1):0.1);");
        String filename = "/home/ntmwn/research/diffusion/bridge/sims/sample5_0.05.txt";
        double sig = 0.05;
        treebase.Tree[] trees = sampleRandomWalk(x0, sig*sig, numStepsPerEdge, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, x0);

        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample5_0.1.txt";
        sig = 0.1;
        trees = sampleRandomWalk(x0, sig*sig, numStepsPerEdge, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, x0);

        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample5_0.2.txt";
        sig = 0.2;
        trees = sampleRandomWalk(x0, sig*sig, numStepsPerEdge, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, x0);
*/
        
        /*
        double sd = 0.5;
        treebase.Tree x0 = new treebase.Tree("((B:1.6586161,0:1.0000000):0.005,(D:1.8688598,(C:1.8704023,A:1.6019577):0.01):0.005);");
        x0.removeDegreeTwoVertices();
        treebase.Tree[] sample = sampleRandomWalk(x0, sd*sd, numStepsPerEdge, false, numSamples);
        doIO(sample, "/home/ntmwn/research/diffusion/bridge/analysis/gene-trees/sample_from_model.txt");
                
        treebase.Tree tree12exp = new treebase.Tree("(((10:0.05573371849,11:0.06642649375):0.2177203249,((1:0.1283550009,5:0.04093855931):0.05870581456,((2:0.03039954835,((7:0.1999242181,3:0.2111989623):0.07716522449,6:0.08479326898):0.161405007):0.01258620429,12:0.1022567336):0.01441222775):0.3385212578):0.04243862512,((9:0.03129900956,8:0.01618438188):0.2239585316,4:0.03890367225):0.0413265734);");
        treebase.Tree tree12gam = new treebase.Tree("((2:0.06940520928,((9:0.09499913099,6:0.1571981876):0.1954443317,10:0.1999821016):0.06664133088):0.046020603,((8:0.09139039396,5:0.1143910741):0.03315466986,((4:0.1319361885,(1:0.09569206748,11:0.05746595255):0.05218749854):0.1400176217,((7:0.2224401651,3:0.1437059455):0.05803236345,12:0.06408690304):0.1576082849):0.05489020873):0.04642169123);");
        treebase.Tree tree48gam = new treebase.Tree("((32:0.1927397277,42:0.08967935642):0.3471614028,((((((23:0.2345202743,25:0.2180703431):0.1479831202,((27:0.06739848559,8:0.0663348358):0.08766451388,(((45:0.0575773101,34:0.05165455834):0.04623604109,(2:0.007416652952,31:0.08326316312):0.06188288899):0.02936877451,12:0.08356665788):0.08925752262):0.06909573967):0.1018810664,18:0.09270694558):0.1464025609,(29:0.02936599367,(33:0.3721056951,(7:0.2225955202,6:0.04868868004):0.02759559462):0.1308542598):0.04608838515):0.08819969861,(19:0.08144320552,((((47:0.09632498525,30:0.1381094247):0.2227428254,48:0.09283881166):0.2371501551,28:0.09275759878):0.160036746,3:0.06583195405):0.02157228917):0.0144296224):0.1324568182,(((((46:0.3471614028,39:0.1927397277):0.08158717432,(37:0.09547001175,38:0.1324568182):0.08967935642):0.1321578983,(((40:0.1479831202,(22:0.2180703431,26:0.06909573967):0.2345202743):0.1018810664,((5:0.0663348358,44:0.08925752262):0.06739848559,35:0.02936877451):0.08766451388):0.1464025609,(((43:0.06188288899,21:0.007416652952):0.05165455834,(1:0.08356665788,10:0.09270694558):0.08326316312):0.0575773101,36:0.04608838515):0.04623604109):0.08819969861):0.1291252694,(((((4:0.04868868004,17:0.0144296224):0.2225955202,41:0.08144320552):0.02759559462,(9:0.160036746,24:0.2371501551):0.02157228917):0.3721056951,(20:0.09632498525,13:0.1381094247):0.2227428254):0.1308542598,(14:0.09275759878,15:0.06583195405):0.09283881166):0.02936599367):0.2295623152,(11:0.2295623152,16:0.1291252694):0.02218412876):0.02218412876):0.09547001175);");
        tree12exp.removeDegreeTwoVertices();
        tree12gam.removeDegreeTwoVertices();
        tree48gam.removeDegreeTwoVertices();
        
        String filename;
        double sig;
        treebase.Tree[] trees;
        
        Random.setEngine(24848);
        
        treebase.Tree yeastTree = new treebase.Tree("(((((Calb:0.1,Sklu:0.3639600):0.0813599,Scas:0.3693700):0.2361292,Sbay:0.0906500):0.0047396,Skud:0.0757200):0.0112125,(Smik:0.0766800,(Scer:0.0588900,Spar:0.0361500):0.0280271):0.0112125);");
        yeastTree.removeDegreeTwoVertices();
        filename = "/home/ntmwn/research/diffusion/bridge/analysis/yeast_sim.txt";
        sig = 0.12;
        trees = sampleRandomWalk(yeastTree, sig*sig, numStepsPerEdge, false, numSamples);
        doIO(trees, filename);
        topologicalSummary(trees, yeastTree);
                */
        
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.05.txt";
//        sig = 0.05;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//        topologicalSummary(trees, tree48gam);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.1.txt";
//        sig = 0.1;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//        topologicalSummary(trees, tree48gam);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.2.txt";
//        sig = 0.2;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//        topologicalSummary(trees, tree48gam);

//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.0001.txt";
//        sig = 0.0001;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_exp_0.01.txt";
//        sig = 0.01;
//        trees = sampleRandomWalk(tree12exp, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_exp_0.025.txt";
//        sig = 0.025;
//        trees = sampleRandomWalk(tree12exp, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_exp_0.05.txt";
//        sig = 0.05;
//        trees = sampleRandomWalk(tree12exp, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_exp_0.1.txt";
//        sig = 0.1;
//        trees = sampleRandomWalk(tree12exp, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_gam_0.01.txt";
//        sig = 0.01;
//        trees = sampleRandomWalk(tree12gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_gam_0.025.txt";
//        sig = 0.025;
//        trees = sampleRandomWalk(tree12gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_gam_0.05.txt";
//        sig = 0.05;
//        trees = sampleRandomWalk(tree12gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample12_gam_0.1.txt";
//        sig = 0.1;
//        trees = sampleRandomWalk(tree12gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.01.txt";
//        sig = 0.01;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.015.txt";
//        sig = 0.015;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.02.txt";
//        sig = 0.02;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.03.txt";
//        sig = 0.03;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
//
//        filename = "/home/ntmwn/research/diffusion/bridge/sims/sample48_gam_0.04.txt";
//        sig = 0.04;
//        trees = sampleRandomWalk(tree48gam, sig*sig, numStepsPerEdge, false, numSamples);
//        doIO(trees, filename);
        

    }
}

