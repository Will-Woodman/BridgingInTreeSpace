/*
    StarTree
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

/** Class representing the start tree for a fixed set of taxa. 
  Used so that the star tree can be rapidly cloned in order to produce central 
  gaussian tree distributions. 
 */

import java.util.Iterator;
import treebase.AlgorithmError;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;

    public class StarTree {
        
        protected Tree theTree;
        protected TreeAsSplits theTreeAsSplits;
    
        protected StarTree() {
            // Exists only to defeat instantiation.

        }

        private static class StarTreeSingletonHolder {
            public static final StarTree INSTANCE = new StarTree();
        }

        public static StarTree getInstance() {
            return StarTreeSingletonHolder.INSTANCE;
        }
        
        public void setTree(TreeAsSplits x) throws AlgorithmError {
            theTreeAsSplits = new TreeAsSplits(x.getNumTaxa());
            Iterator<Split> it = x.getSplitIterator();
            while (it.hasNext()) {
                Split p = it.next();
                if (p.isTerminal()!=null) {
                    // Triv split
                    theTreeAsSplits.add(p, 1.0);
                }
            }
            // Sanity check
            if (theTreeAsSplits.getNumSplits()!=theTreeAsSplits.getNumTaxa()) {
                throw new AlgorithmError("Strange problem generating a star tree: missing splits.");
            }
            
            theTree = theTreeAsSplits.getTree();
        }
        
        public Tree getTree() {
            return theTree.clone();
        }
        
        public TreeAsSplits getTreeAsSplits() {
            return theTreeAsSplits.efficientClone();
        }
    
    
}
