/*
BridgeWithNNIStepLike
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

import geodesics.Geodesic;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import simulation.CategoricalDistribution;
import simulation.NormalDistribution;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.TreeAsSplits;

/** Implementation of the forward bridge using NNI likelihood obtained from a 
    simple geodesic at each step.
 */

public class BridgeWithNNIStepLike extends ForwardStepBridge {
    
    public BridgeWithNNIStepLike(TreeAsSplits tA, TreeAsSplits tB, int n) {
        super(tA, tB, n);
    }
    
    @Override
    /** Make empty copy for use in MCMC proposal */
    public ForwardStepBridge makeEmptyCopy() {
        return new BridgeWithNNIStepLike(treeA, treeB, numSteps);
    }
    @Override
    public ForwardStepBridge makeFreshBridge(TreeAsSplits tA, TreeAsSplits tB) {
        return new BridgeWithNNIStepLike(tA, tB, numSteps);
    }
    @Override
    public ForwardStepBridge makeEmptyCopy(TreeAsSplits x, boolean whichEnd) {
        if (whichEnd) return new BridgeWithNNIStepLike(treeA, x, numSteps);
        else return new BridgeWithNNIStepLike(x, treeB, numSteps);
    }

    
    /** Compute the log likelihood associated with a particular step. 
     Return [0] = log like; [1] = log prop density. */
    protected double[] calcStepLogLike(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException {
        
        int nprime = tA.getNumTaxa()-3;
        double likeStepSD = Math.sqrt(t0/numSteps);
        if (!g.isSimple()) {
            throw new AlgorithmException("Non-simple segment in BridgeWithNNIStepLike");
        }

        /* Perform N' steps linking the two trees */
        HashMap<Split, Split> splitMap =  g.getSimpleSplitMap();
        /* Remove pendant edges */
        Iterator<Map.Entry<Split,Split>> it = splitMap.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<Split,Split> pair = (Map.Entry)it.next();
            if (pair.getKey().isTerminal()!=null)
                it.remove(); 
        }
            
        /* Sample from key set uniformly at random, N' times */
        
        double ll = 0.0; // log like
        double lpd = 0.0; // log prop density
            
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(splitMap.keySet());
        TreeAsSplits stepTree = tA.efficientClone();
        for (int i=0; i<(tA.getNumTaxa()-3); i++) {
            ArrayList<Split> compatible = new ArrayList();
            Iterator<Split> itP = splits.iterator();
            while (itP.hasNext()) {
                Split p, q;
                p = itP.next();
                q = splitMap.get(p);
              
                boolean good = true;
                if (!stepTree.contains(q)) {
                    /* Iterate thru splits in stepTree (except p) and ensure compatible with q */
                    Iterator<Split> itStep = stepTree.getSplitIterator();

                    while (itStep.hasNext()) {
                        Split r = itStep.next();
                        if (r!=p) {
                            if (!r.isCompatible(q)) {
                                good = false;
                                break;
                            }
                        }
                    }
                }
                if (good) compatible.add(p);    
            } // End loop thru splits we might replace
                
            Split p = (Split) CategoricalDistribution.sampleFromListWithReplacement(compatible, unifDist); // Doesn't matter if you replace or not. Replace is faster!
            lpd -= Math.log(compatible.size());
            Split q = splitMap.get(p);
            double pl = tA.getSplitLength(p);
            double ql = tB.getSplitLength(q);
            stepTree.remove(p);
 //           stepTree.addWithCheck(q, ql); // NB: CHANGE THIS TO add AFTER DEBUGGING
            stepTree.add(q, ql);
            splits.remove(p);
                
            /* Compute loglike */
            if (p.equals(q)) {
                ll += Math.log(NormalDistribution.pdf(ql-pl, 0.0, likeStepSD)+NormalDistribution.pdf(pl+ql, 0.0, likeStepSD)/3.0);
            }
            else {
                ll += NormalDistribution.logpdf(pl+ql, 0.0, likeStepSD)-Math.log(3.0); // change of splits
            }
        }        

        return new double[] {ll, lpd};    
    }
    
}
