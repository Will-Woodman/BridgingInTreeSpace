/*
 * GeodesicForBridging.java

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

import bridge.ForwardStepBridge.IntDoublePair;
import cern.jet.random.tdouble.DoubleUniform;
import geodesics.Geodesic;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.TreeAsSplits;

/**
    An extension of the basic BHV geodesic, with methods added specifically for bridging.
    NB: after debugging change all "addwithcheck" to "add".
 */

public class GeodesicForBridging extends Geodesic {
    
    public GeodesicForBridging(TreeAsSplits tA, TreeAsSplits tB) throws AlgorithmError {
        super(tA, tB);
    }

    public double[] getCriticalLambda() {
        int k = partitionA.size();
        double[] criticalLambda = new double[k];
        double na, nb;
        HashSet<Split> ai, bi;
        for (int i=0; i<k; i++) {
            ai = partitionA.get(i);
            bi = partitionB.get(i);
            na = Math.sqrt(edgeLenNorm(ai, treeA));
            nb = Math.sqrt(edgeLenNorm(bi, treeB));
            criticalLambda[i] = na/(na+nb);
        }
        return criticalLambda;
    }
    
    public int[][] getPartitionSizes() {
        int[][] sizes = new int[2][partitionA.size()];
        for (int i=0; i<partitionA.size(); i++) {
            sizes[0][i] = partitionA.get(i).size();
            sizes[1][i] = partitionB.get(i).size();
        }
        return sizes;
    }
    
    public TreeAsSplits getLimitTree(int k, boolean dir)  {
        
        TreeAsSplits theTree = new TreeAsSplits(treeA.getNumTaxa());
        double na = Math.sqrt(edgeLenNorm(partitionA.get(k), treeA));
        double nb = Math.sqrt(edgeLenNorm(partitionB.get(k), treeB));
        double s = na/(na+nb);
        
        try {
        /* Assemble topology (i.e. the set of splits):
           first "unshared" splits, which come from the partition A_i and B_i */
 
        /* Loop thru parts of the B partition which we need to add:
           add in B_1 up to B_(orthontIndex) */
        if (!dir) {
            for (Split p : partitionB.get(k)) {
                theTree.addWithCheck(p, 0.0);
            }
        }
        for (int i=0; i<k; i++) {
            HashSet<Split> bi = partitionB.get(i);
            // Loop thru Bi: work out length and add to treeMap
            double nai = Math.sqrt(edgeLenNorm(partitionA.get(i), treeA));
            double nbi = Math.sqrt(edgeLenNorm(bi, treeB));
            for (Split p : bi) {
                double e = treeB.getSplitLength(p);
                double u = (s*nbi-(1.0-s)*nai)*e/nbi;
                theTree.addWithCheck(p, u);
            }
        }
        /* Loop thru parts of the A partition which we need to add:
           add A_k with zero edge length */
        if (dir) {
            for (Split p : partitionA.get(k)) {
                theTree.addWithCheck(p, 0.0);
            }
        }
        /* Add the other A_i in */
        for (int i=k+1; i<partitionA.size(); i++) {
            HashSet<Split> ai = partitionA.get(i);
            // Loop thru Ai: work out length and add to treeMap
            double nai = Math.sqrt(edgeLenNorm(ai, treeA));
            double nbi = Math.sqrt(edgeLenNorm(partitionB.get(i), treeB));
            for (Split p : ai) {
                double e = treeA.getSplitLength(p);
                double u = ((1.0-s)*nai-s*nbi)*e/nai;
                theTree.addWithCheck(p, u);
            }
        }

        /* Remains to do the shared splits -- easier! */
        double u, v, e;
        for (Split p : sharedSplits) {
            if (treeA.contains(p)) u = treeA.getSplitLength(p);
            else u = 0.0;
            if (treeB.contains(p)) v = treeB.getSplitLength(p);
            else v = 0.0;
            e = (1.0-s)*u+s*v; // Linear scale on shared splits!
            theTree.addWithCheck(p, e);
        }
        
        }
        catch (AlgorithmException anEx) {
            System.out.println("Error in getLimitTree: "+anEx.getMessage());
        }
        
        return theTree;
    }
    

}
