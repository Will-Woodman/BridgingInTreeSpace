 /*
    TreeDistributions
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

import cern.jet.random.tdouble.DoubleUniform;
import geodesics.Geodesic;
import java.io.File;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import simulation.NormalDistribution;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;

/**
 * Static methods for simulating distributions on tree-space
 */

public class TreeDistributions {
    
    static private final double defaultRadiusFac = 0.33; // These params are used to control the default samplers
    static private final int defaultStepsPerEdge = 10;
    
    /** The following methods all sample from the uniform distribution on a ball
     centred at x0 and radius r. 
     The methods all rely on an MH-MCMC scheme where the innovation is the 
     standard random walk single step (symmetric Guassian RW).
     Since the distribution is uniform, the acc/rej step is particularly simple:
     reject if the proposed tree lies outside the ball. 
     An approx upper bound for the distance from x0 is maintained at all times. 
     This saves computing the actual distance at every iteration, */
    
    /** Obtain a single sample, using default values for sd and number of steps. */
    static public TreeAsSplits sampleUniformBall(TreeAsSplits x0, double r, NormalDistribution normDist, DoubleUniform unifDist, boolean[] crossedBoundary) {
        return sampleUniformBall(x0, r, defaultRadiusFac*r, defaultStepsPerEdge*(x0.getNumTaxa()-3), normDist, unifDist, crossedBoundary);
    }
    
    /** Obtain a single sample. */
    static public TreeAsSplits sampleUniformBall(TreeAsSplits x0, double r, double sd, int nits, NormalDistribution normDist, DoubleUniform unifDist, boolean[] crossedBoundary) {
        ArrayList<TreeAsSplits> res =  sampleUniformBall(x0, r, sd, nits, 0, -1, normDist, unifDist, crossedBoundary);
        return res.get(0);
    }

    /** Run MCMC, but no need to pass in boundary cross detection.
     This is for use in applications other than the bridge algorithms. */
    static public ArrayList<TreeAsSplits> sampleUniformBall(TreeAsSplits x0, double r, double sd, int burnits, int thin, int nits, NormalDistribution normDist, DoubleUniform unifDist) {
        boolean[] crossedBoundary = new boolean[1];
        return sampleUniformBall(x0, r, sd, burnits, thin, nits, normDist, unifDist, crossedBoundary);
    }
    
    /** Run MH-MCMC algorithm to sample from ball. 
     Maintain an estimate of distance from x0. 
     
     crossedBoundary[] is used to return boolean indicating whether any orthant boundaries were crossed. 
     This is important for the bridge algorithms. 
     
     Pass in thin = -1 if you only want the final tree
     */
    static public ArrayList<TreeAsSplits> sampleUniformBall(TreeAsSplits x0, double r, double sd, int burnits, int thin, int nits, NormalDistribution normDist, DoubleUniform unifDist, boolean[] crossedBoundary) {
        
        /* The set of trees to return */
        ArrayList<TreeAsSplits> res = new ArrayList();
        
        TreeAsSplits current = x0.efficientClone();
        double bound = 0.0, propBound;
        normDist.setParams(0.0, sd);
        
       /* Get the non-trivial splits in a fixed order */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(current.getNonTrivialSplits());
        // Sort the splits to obtain reproducible results -- but time consuming!
//        java.util.Collections.sort(splits);
        if (splits.size()<current.getNumTaxa()-3) {
            System.out.println("Warning: request made to randomly walk an unresolved tree in uniformball sampler. Trees must be randomly resolved first.");
        }
        crossedBoundary[0] = false;
        
        /* Main MCMC loop */
        int splitInd, topInd, acc = 0;
        Split p, q;
        double x, l, bl;
        Geodesic g = null;
        for (int i=0; i<(nits+burnits); i++) {
            
            /* Pick a split at random */
            splitInd = unifDist.nextIntFromTo(0, splits.size()-1);
            p = splits.get(splitInd);
            /* Perturb edge length */
            x = normDist.sample();
            bl = current.getSplitLength(p); // Store current edge length
            l = x+bl;
            propBound = bound+Math.abs(x); // Use triangle inequality to get max possible distance from x0
            q = null;

            try {
                if (l>0.0) { // No boundary crossed
                    current.setSplitLength(p, l);
                }
                else {
                    topInd = unifDist.nextIntFromTo(1, 3);
                    if (topInd==1) {
                        current.setSplitLength(p, -l);
                    }
                    else { // Crossed boundary -- perform NNI
                        Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(current, p);
                        q = nniSplit[topInd-2]; // topInd-2 = 0 or 1
                        current.remove(p);
                        current.add(q, -l);
                        splits.set(splitInd, q);
                    }
                }
            }
            catch (AlgorithmException err) {
                System.out.println("Error setting edge length in uniform ball sampler.");
            }
            
            /* Accept / reject */
            
            if (propBound<r) {
                // Accept
                bound = propBound;
                if (q!=null) crossedBoundary[0] = true;
                acc++;
            }
            else { // Approx bound says you might be outside ball, so...
                // Compute distance from x0
                if (g==null) {
                    try {
                        g = new Geodesic(current, x0);
                    } catch (AlgorithmError ex) {
                        System.out.println("Error making geodesic in uniform ball sampler.");
                    }
                }
                else {
                    g.shiftSingleEnd(false, current);
                }
                // Update bound wth computed value
                propBound = g.getInternalLength();
                if (propBound<r) {
                    // Accept
                    bound = propBound;
                    if (q!=null) crossedBoundary[0] = true;
                    acc++;
                }
                else {
                    // Reject -- reverse move on tree
                    if (q!=null) {
                        current.remove(q);
                        current.add(p, bl);
                        splits.set(splitInd, p);
                    }
                    else {
                        try {
                            current.setSplitLength(p, bl);
                        } catch (AlgorithmError ex) {
                            System.out.println("Error setting edge length in uniform ball sampler.");
                        }
                    }
                }
            }
            
            if (i>=burnits & (thin>0) & ((i-burnits+1)%thin)==0) {
                res.add(current.efficientClone());
            }
             
        } // End main loop
        
        if (thin<0) res.add(current);

        double accProp = ((double) acc)/((double) nits+burnits);
        System.out.println("Acceptance prop = "+String.format("%7.7f",accProp));
        
        return res;
    }
    
    
/** The following methods all sample from the exponential bump distribution
     centred at x0 and "variance" sigma.
     
     f(x)\propto exp -1/(2sigma^2)d(x,x_0)^2
     The methods all rely on an MH-MCMC scheme where the innovation is the 
     standard random walk single step (symmetric Guassian RW).
     Since acc/rej step is simple: the prop ratio cancels (symmetry),
     so the acc prob is just the density ratio, and the normalizing const cancels. 
*/
    
/** Run MCMC, but no need to pass in boundary cross detection.
     This is for use in applications other than the bridge algorithms. */
    static public ArrayList<TreeAsSplits> sampleExponentialBump(TreeAsSplits x0, double sigma, double sd, int burnits, int thin, int nits, NormalDistribution normDist, DoubleUniform unifDist) {
        boolean[] crossedBoundary = new boolean[1];
        return sampleExponentialBump(x0, sigma, sd, burnits, thin, nits, normDist, unifDist, crossedBoundary);
    }
    

    /** Run MH-MCMC algorithm to sample from exponential bump. 
     crossedBoundary[] is used to return boolean indicating whether any orthant boundaries were crossed. 
     This is important for the bridge algorithms. 
     
     Pass in thin = -1 if you only want the final tree
     */
    static public ArrayList<TreeAsSplits> sampleExponentialBump(TreeAsSplits x0, double sigma, double sd, int burnits, int thin, int nits, NormalDistribution normDist, DoubleUniform unifDist, boolean[] crossedBoundary) {
        
        /* The set of trees to return */
        ArrayList<TreeAsSplits> res = new ArrayList();
        
        TreeAsSplits current = x0.efficientClone();
        normDist.setParams(0.0, sd);
        
       /* Get the non-trivial splits in a fixed order */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(current.getNonTrivialSplits());
        // Sort the splits to obtain reproducible results -- but time consuming!
//        java.util.Collections.sort(splits);
        if (splits.size()<current.getNumTaxa()-3) {
            System.out.println("Warning: request made to randomly walk an unresolved tree in exponentialbump sampler. Trees must be randomly resolved first.");
        }
        crossedBoundary[0] = false;
        
        /* Main MCMC loop */
        int splitInd, topInd, acc = 0;
        double f = -0.5/sigma/sigma;
        Split p, q;
        double x, l, bl;
        double dxy=0.0, olddxy;
        Geodesic g = null;
        for (int i=0; i<(nits+burnits); i++) {
            
            olddxy = dxy; // store distance from x0
            
            /* Pick a split at random */
            splitInd = unifDist.nextIntFromTo(0, splits.size()-1);
            p = splits.get(splitInd);
            /* Perturb edge length */
            x = normDist.sample();
            bl = current.getSplitLength(p); // Store current edge length
            l = x+bl;
            q = null;

            try {
                if (l>0.0) { // No boundary crossed
                    current.setSplitLength(p, l);
                }
                else {
                    topInd = unifDist.nextIntFromTo(1, 3);
                    if (topInd==1) {
                        current.setSplitLength(p, -l);
                    }
                    else { // Crossed boundary -- perform NNI
                        Split[] nniSplit = treedatasets.OperationsOnTreeAsSplits.getNNISplits(current, p);
                        q = nniSplit[topInd-2]; // topInd-2 = 0 or 1
                        current.remove(p);
                        current.add(q, -l);
                        splits.set(splitInd, q);
                    }
                }
            }
            catch (AlgorithmException err) {
                System.out.println("Error setting edge length in exponential bump sampler.");
            }
            
            /* Accept / reject */
            
            // Compute distance from x0
            if (g==null) {
                try {
                    g = new Geodesic(current, x0);
                } catch (AlgorithmError ex) {
                    System.out.println("Error making geodesic in exponential bump sampler.");
                }
            }
            else {
                g.shiftSingleEnd(false, current);
            }
            
            // Calc log density ratio
            dxy = g.getInternalLength();
            double logDensityRatio = f*dxy*dxy-f*olddxy*olddxy;
            if (Math.log(unifDist.nextDouble())<logDensityRatio) {
                // Accept
                acc++;
            }
            else {
                // Reject
                dxy = olddxy;
                // reverse move on tree
                if (q!=null) {
                    current.remove(q);
                    current.add(p, bl);
                    splits.set(splitInd, p);
                }
                else {
                    try {
                        current.setSplitLength(p, bl);
                    } catch (AlgorithmError ex) {
                        System.out.println("Error setting edge length in exponential bump sampler.");
                    }
                }
            }
            
            if (i>=burnits & (thin>0) & ((i-burnits+1)%thin)==0) {
                res.add(current.efficientClone());
            }
             
        } // End main loop
        
        if (thin<0) res.add(current);

        double accProp = ((double) acc)/((double) nits+burnits);
        System.out.println("Acceptance prop = "+String.format("%7.7f",accProp));
        
        return res;
    }
    
       /** Run MH-MCMC algorithm to sample from exponential bump - this time using GGF as proposal**/
    
    static public ArrayList<TreeAsSplits> sampleExponentialBumpViaMVNStep(TreeAsSplits x0, double sigma, double sd, int burnits, int thin, int nits, NormalDistribution normDist, DoubleUniform unifDist, boolean[] crossedBoundary) {
        
        /* The set of trees to return */
        ArrayList<TreeAsSplits> res = new ArrayList();
        
        TreeAsSplits current = x0.efficientClone();
        normDist.setParams(0.0, sd);
        
       /* Get the non-trivial splits in a fixed order */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(current.getNonTrivialSplits());
        // Sort the splits to obtain reproducible results -- but time consuming!
//        java.util.Collections.sort(splits);
        if (splits.size()<current.getNumTaxa()-3) {
            System.out.println("Warning: request made to randomly walk an unresolved tree in exponentialbump sampler. Trees must be randomly resolved first.");
        }
        crossedBoundary[0] = false;
        
        /* Main MCMC loop */
        int splitInd, topInd, acc = 0;
        double f = -0.5/sigma/sigma;
        double dxy=0.0, olddxy;
        Geodesic g = null;
        double[] infor= new double[3];
        for (int i=0; i<(nits+burnits); i++) {
            
            olddxy = dxy; // store distance from x0
            TreeAsSplits old=current.clone();


            try {
                sampleDirectMVN(current, sd, normDist, unifDist, infor);
            }
            catch (AlgorithmException err) {
                System.out.println("Error setting edge length in exponential bump sampler.");
            }
            
            /* Accept / reject */
            
            // Compute distance from x0
            if (g==null) {
                try {
                    g = new Geodesic(current, x0);
                } catch (AlgorithmError ex) {
                    System.out.println("Error making geodesic in exponential bump sampler.");
                }
            }
            else {
                g.shiftSingleEnd(false, current);
            }
            
            // Calc log density ratio
            dxy = g.getInternalLength();
            double logDensityRatio = f*dxy*dxy-f*olddxy*olddxy;
            if (Math.log(unifDist.nextDouble())<logDensityRatio) {
                // Accept
                acc++;
            }
            else {
                // Reject
                dxy = olddxy;
                // reverse move on tree
               current=old;
            }
            
            if (i>=burnits & (thin>0) & ((i-burnits+1)%thin)==0) {
                res.add(current.efficientClone());
            }
             
        } // End main loop
        
        if (thin<0) res.add(current);

        double accProp = ((double) acc)/((double) nits+burnits);
        System.out.println("Acceptance prop = "+String.format("%7.7f",accProp));
        
        return res;
    }
    
    
    /* MVN sampler ---------------------------------------------------------- */
    
    /** Private class for sorting edges -- used by MVN sampler */
    static private class IntDoublePair implements Comparable {
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

    
    /** Perturb tree via "direct" MVN Step. --GGF
        Return log density.
        Argument sigma is SD of isotropic multivariate normal distro used to sample direction and distance.
    */
    static public void sampleDirectMVN(TreeAsSplits theTree, double sigma, NormalDistribution normDist, DoubleUniform unifDist, double[] info) throws AlgorithmException {

        normDist.setParams(0.0, sigma);
        int nprime = theTree.getNumInternalSplits();
        double logDens = 0.0;
        
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
        logDens = -nprime*Math.log(sigma)-0.5*lenSqu/(sigma*sigma); //  multivariate normal density

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
            if ((l[i]<0.0)&&(u[i]<0.0)) boundaries.add(new IntDoublePair(i,-(y*len)/u[i])); // -y*len/u[i] is the "time" when you hit the boundary. Order by time!
        }
        
        info[0] = logDens;
        info[1] = -0.693147*boundaries.size(); //log(1/2)
        info[2] = len;       
        
        if (boundaries.size()==0) {
            // Nothing to do -- return
            return;
        }
        
        /* We crossed at least one boundary */
        java.util.Collections.sort(boundaries);
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
             
        return;
    }
    
    
    
    
    /* might delete this
    
    public static TreeAsSplits sampleCentralGaussian(double sd, NormalDistribution normDist, DoubleUniform unifDist,double[] info) {
        StarTree s = StarTree.getInstance();
        TreeAsSplits theTree = s.getTreeAsSplits();
        double logDens=0.0;
        double lensqu=0.0;
        double len=0.0;
        double numOfOrthants=1;
        try {
            double l = TreeResolver.resolveTree(theTree, s.getTree(), unifDist);
        } catch (AlgorithmException ex) {
            System.out.println("Error resolving star tree.");
        }
        try {
            Iterator<Split> it = theTree.getSplits().iterator();
            while (it.hasNext()) {
                Split p = it.next();
                if (p.isTerminal()==null) {
                    theTree.setSplitLength(p, Math.abs(normDist.sample(0.0, sd)));
                    lensqu+=theTree.getSplitLength(p)*theTree.getSplitLength(p);
                    logDens+=Math.log(sd)-0.5*theTree.getSplitLength(p)/(sd*sd);
                    len+=theTree.getSplitLength(p);
                }
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error setting edge length when sampling central Gaussian.");
        }
        
        for(int i=1;i<2*theTree.getNumTaxa()-5;i++){
            if((double)i/2!=i/2){
            numOfOrthants=numOfOrthants*i;
        }
        }
        
        info[0] = logDens;
        info[1] = Math.log(1/numOfOrthants); 
        info[2] = len;  
        
        //HashMap<TreeAsSplits,double[]> Result = new HashMap();
        //Result.put(theTree,info);
        
        
        return theTree;
    }
    
    */
    
     
    
    /* Test area */
    public static void main(String[] args) throws AlgorithmError {
        
        // Initial tree x0 
        String initialTreeFilename = "/home/ntmwn/research/diffusion/bridge/sims/randomtree_48_gamma_2_20.txt";
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
        
        double r = 0.5;
        double sd = 1.5;
        int burnits = 0, thin = 10, nits = 20000;
        NormalDistribution normDist = new NormalDistribution(0.0, sd);
        DoubleUniform unifDist = new DoubleUniform(simulation.Random.getEngine());
       
        ArrayList<TreeAsSplits> sample = sampleExponentialBump(x0, r, sd, burnits, thin, nits, normDist, unifDist);

        for (int i=1; i<sample.size(); i++) {
            Geodesic h = new Geodesic(sample.get(i), x0);
            System.out.println(String.format("%7.7f", h.getInternalLength()));
        }

    }

}
