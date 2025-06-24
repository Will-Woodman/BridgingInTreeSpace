/*
BridgeWithApproxMVNLike
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import simulation.NormalDistribution;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.TreeAsSplits;

/** Implementation of the forward bridge using approx MVN log like obtained from a 
    simple geodesic at each step (Also called Gaussian via Geodesic firing)
 */

public class BridgeWithApproxMVNLike extends ForwardStepBridge {
    
    public BridgeWithApproxMVNLike(TreeAsSplits tA, TreeAsSplits tB, int n) {
        super(tA, tB, n);
    }
    
    @Override
    /** Make empty copy for use in MCMC proposal */
    public ForwardStepBridge makeEmptyCopy() {
        return new BridgeWithApproxMVNLike(treeA, treeB, numSteps);
    }
    @Override
    public ForwardStepBridge makeFreshBridge(TreeAsSplits tA, TreeAsSplits tB) {
        return new BridgeWithApproxMVNLike(tA, tB, numSteps);
    }
    @Override
    public ForwardStepBridge makeEmptyCopy(TreeAsSplits x, boolean whichEnd) {
        if (whichEnd) return new BridgeWithApproxMVNLike(treeA, x, numSteps);
        else return new BridgeWithApproxMVNLike(x, treeB, numSteps);
    }

    /** Compute the log likelihood associated with a particular step. 
     Return [0] = log like; [1] = log prop density. 
     Assumes both trees fully resolved and linked by a simple geodesic. */
    protected double[] calcStepLogLike(TreeAsSplits tA, TreeAsSplits tB, Geodesic g, double t0) throws AlgorithmException {
        int nprime = tA.getNumTaxa()-3;
        double likeStepSD = Math.sqrt(t0/numSteps);
        if (!g.isSimple()) {
            throw new AlgorithmException("Non-simple segment in BridgeWithApproxMVNLike");
        }
        
        double ll = 0.0;
        HashMap<Split, Split> splitMap =  g.getSimpleSplitMap();
        
        Iterator<Split> it = tA.getSplitIterator();
        while (it.hasNext()) {
            Split p = it.next();
            if (p.isTerminal()==null) {
                double la = tA.getSplitLength(p);
                if (!tB.contains(p)) {
                    // This split contracts and expands  
                    double lb = tB.getSplitLength(splitMap.get(p));
                    //GGF likelihood
                    ll += NormalDistribution.logpdf(la+lb, 0.0, likeStepSD)-Math.log(2);
                    
                }
                else {
                    // Shared split
                    
                    double lb = tB.getSplitLength(p);
                    //GGF likelihood
                    ll += NormalDistribution.logpdf(la-lb, 0.0, likeStepSD);
                } 
            }
        }

        return new double[] {ll, 0.0};
        
    }
    
}
