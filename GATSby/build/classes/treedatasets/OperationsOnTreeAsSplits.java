/*
 * OperationsOnTreeAsSplits.java

    Copyright (C) 2012  Tom M. W. Nye

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

package treedatasets;

/**
 * Set of static methods for performing operations on TreeAsSplits objects
 *
 * Operations:
 * 1. Randomly resolving
 * 2. Given fully resolved tree and an edge (split) find splits obtained by NNI on that edge
 */

import treebase.TreeAsSplits;
import java.util.*;
import treebase.Split;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.Graph;
import cern.jet.random.tdouble.DoubleUniform;
import simulation.Random;

public class OperationsOnTreeAsSplits {

    public static final boolean DEBUG = true;

    /** Resolve a tree by adding in splits from theSplits.
     There's no rationale over the order of addition: just keep adding until
     no more splits can go in.
     Method is based on a hashset and so is not deterministic. */
    public static void resolve(TreeAsSplits theTreeMap, HashSet<Split> theSplits) {
        HashSet<Split> addedSplits = new HashSet();
        if (theSplits!=null) {
            Iterator<Split> it = theSplits.iterator();
            for (int i=0; i<theSplits.size(); i++) {
                Split p = it.next();
                if ((theTreeMap.checkCompatibility(p))&&(!theTreeMap.contains(p))) {
                    theTreeMap.add(p, 1.0);
                    addedSplits.add(p);
                }
            }
        }

        // Set all added splits to have length zero
        try {
            Iterator<Split> it = addedSplits.iterator();
            for (int i=0; i<addedSplits.size(); i++) {
                theTreeMap.setSplitLength(it.next(), 0.0);
            }
        }
        catch (AlgorithmException anErr) {
            System.out.println("Error resolving a tree: seems like a split wasn't added."); // Shouldn't be possible
        }
    }

    public static void resolve(HashSet<Split> treeSplits, HashSet<Split> refSplits) {
        if (refSplits!=null) {
            Iterator<Split> it = refSplits.iterator();
            for (int i=0; i<refSplits.size(); i++) {
                Split p = it.next();
                if ((p.testCompatibilityWithSetofSplits(treeSplits))&&(!treeSplits.contains(p))) {
                    treeSplits.add(p);
                }
            }
        }
    }

    /** Resolve a tree by randomly resolving degree > 3 vertices
     Return log density value. 
     Pass back the set of splits created, and set the length of all
     splits you added to have zero length. 
     Pass in newSplits=null, if not needed
     This method is inefficient: relies on converting the TreeAsSplits to a Tree.
     NB: set up so that all lists are ordered, so should have reproducible results. */
    public static double randomlyResolve(TreeAsSplits theTreeMap) {
        return randomlyResolve(theTreeMap, null);
    }
    public static double randomlyResolve(TreeAsSplits theTreeMap, HashSet<Split> newSplits) {
        
        /* Keep track of which splits we've added */
        HashSet<Split> addedSplits = new HashSet();
        
        /* Store log density */
        double logDens = 0.0;
        
        if (theTreeMap.getNumSplits()!=(theTreeMap.getNumTaxa()-3)) {
            // Still not fully resolved: build a tree, find degree >3 vertices and resolve them
            Tree theTree = null;
            try {
                theTree = theTreeMap.getTree();
            }
            catch (AlgorithmException anErr) {
                System.out.println("Error: problem resolving a tree. This should not be possible. ");
            }

            // Set up random numbers
            DoubleUniform u = new DoubleUniform(Random.getEngine());

            // Find a degree>3 vertex
            Graph.Vertex v;
            do {
                v = null;
                ArrayList<Graph.Vertex> vertices = theTree.getOrderedListOfVertices();
                for (int i=0; i<theTree.numVertices(); i++) {
                    Graph.Vertex w = vertices.get(i);
                    if (w.degree()>3) v=w;
                }
                if (v!=null) {
                    // OK, found a vertex to resolve
                    try {
                        // Resolve this vertex.
                        logDens -= Math.log(0.5*v.degree()*(v.degree()-1));

                        // Start by adding two new vertices
                        Graph.Vertex vA = theTree.addNewVertex("");
                        Graph.Vertex vB = theTree.addNewVertex("");
                        theTree.connect(vA, vB, 1.0); // New branch
                        // Get the new edge as we'll need to find a split for it later
                        Graph.Edge newEdge  = vA.getEdge(vB);

                        // Get two random neighbours & join to vA:
                        // Start by getting an ordered list of neighbours
                        ArrayList<Graph.Vertex> orderedNeighbours = new ArrayList();
                        for (int i=0; i<vertices.size(); i++) {
                            Graph.Vertex w = vertices.get(i);
                            if (w.isConnectedTo(v)) orderedNeighbours.add(w);
                        }
                        // Sanity check
                        if (orderedNeighbours.size()!=v.degree()) System.out.println("Warning: algorithm for resolving a vertex randomly has failed. This should not be possible.");

                        // Get two integers from 0 to deg(v)-1
                        int[] k = new int[2];
                        k[0] = u.nextIntFromTo(0, (v.degree()-1));
                        k[1] = u.nextIntFromTo(0, (v.degree()-2));
                        k[1] = (k[1]<k[0]) ? k[1] : k[1]+1;

                        // Get two random neighbours & join to vA
                        for (int i=0; i<2; i++) {
                            Graph.Vertex w = orderedNeighbours.get(k[i]);
                            Graph.Edge e  = v.getEdge(w);
                            double z = e.getLength();
                            theTree.cut(e);
                            theTree.connect(vA, w, z);
                        }
                        // Join remaining two vertices to vB
                        for (int i=0; i<orderedNeighbours.size(); i++) {
                            if ((i!=k[0])&&(i!=k[1])) {
                                // This vertex should be joined to vB
                                Graph.Vertex w = orderedNeighbours.get(i);
                                Graph.Edge e  = v.getEdge(w);
                                double z = e.getLength();
                                theTree.cut(e);
                                theTree.connect(vB, w, z);
                            }
                        }
                        if ((vA.degree()<3)||(vB.degree()<3)) {
                            System.out.println("Warning: when resolving a tree new vertex has degree less than 3. Trouble ahead...");
                        }
                        theTree.remove(v);

                        // Now get the split.
                        theTree.generateSplits();
                        Split newSplit = theTree.getSplit(newEdge);
                        addedSplits.add(newSplit);
                        theTreeMap.add(newSplit, 1.0);
                    }
                    catch (AlgorithmException anErr) {
                        System.out.println("Error during tree resolution. "+anErr.getMessage());
                    }

                }
            }
            while (v!=null); // If you found an unresolved vertex, then loop thru' again

        }

        // Set all added splits to have length zero
        try {
            Iterator<Split> it = addedSplits.iterator();
            for (int i=0; i<addedSplits.size(); i++) {
                theTreeMap.setSplitLength(it.next(), 0.0);
            }
        }
        catch (AlgorithmException anErr) {
            System.out.println("Error resolving a tree: seems like a split wasn't added.");
        }

        if (newSplits!=null) {
            newSplits.clear();
            newSplits.addAll(addedSplits);
        }
        
        return logDens;
    }

    /** Resolve a tree by randomly resolving degree > 3 vertices
     Return the set of edges created, and set the length of all
     edges you added to have zero length.
    */
    public static HashSet<Graph.Edge> randomlyResolve(Tree theTree) {

        // Set up random numbers
        DoubleUniform u = new DoubleUniform(Random.getEngine());

        // Set up output
        HashSet<Graph.Edge> newEdges = new HashSet();

        // Find a degree>3 vertex
        Graph.Vertex v;
        do {
            v = null;
            ArrayList<Graph.Vertex> vertices = theTree.getOrderedListOfVertices();
            for (int i=0; i<theTree.numVertices(); i++) {
                Graph.Vertex w = vertices.get(i);
                if (w.degree()>3) v=w;
            }
            if (v!=null) {
                // OK, found a vertex to resolve
                try {
                    // Resolve this vertex.

                    // Start by adding two new vertices
                    Graph.Vertex vA = theTree.addNewVertex("");
                    Graph.Vertex vB = theTree.addNewVertex("");
                    theTree.connect(vA, vB, 1.0); // New branch
                    newEdges.add(vA.getEdge(vB));

                    // Get two random neighbours & join to vA:
                    // Start by getting an ordered list of neighbours
                    ArrayList<Graph.Vertex> orderedNeighbours = new ArrayList();
                    for (int i=0; i<vertices.size(); i++) {
                        Graph.Vertex w = vertices.get(i);
                        if (w.isConnectedTo(v)) orderedNeighbours.add(w);
                    }
                    // Sanity check
                    if (orderedNeighbours.size()!=v.degree()) System.out.println("Warning: algorithm for resolving a vertex randomly has failed. This should not be possible.");

                    // Get two integers from 0 to deg(v)-1
                    int[] k = new int[2];
                    k[0] = u.nextIntFromTo(0, (v.degree()-1));
                    k[1] = u.nextIntFromTo(0, (v.degree()-2));
                    k[1] = (k[1]<k[0]) ? k[1] : k[1]+1;

                    // Get two random neighbours & join to vA
                    for (int i=0; i<2; i++) {
                        Graph.Vertex w = orderedNeighbours.get(k[i]);
                        Graph.Edge e  = v.getEdge(w);
                        double z = e.getLength();
                        theTree.cut(e);
                        Graph.Edge f = theTree.connect(vA, w, z);
                        if (newEdges.contains(e)) {
                            newEdges.remove(e);
                            newEdges.add(f);
                        }
                    }
                    // Join remaining two vertices to vB
                    for (int i=0; i<orderedNeighbours.size(); i++) {
                        if ((i!=k[0])&&(i!=k[1])) {
                            // This vertex should be joined to vB
                            Graph.Vertex w = orderedNeighbours.get(i);
                            Graph.Edge e  = v.getEdge(w);
                            double z = e.getLength();
                            theTree.cut(e);
                            Graph.Edge f = theTree.connect(vB, w, z);
                            if (newEdges.contains(e)) {
                                newEdges.remove(e);
                                newEdges.add(f);
                            }
                        }
                    }
                    if ((vA.degree()<3)||(vB.degree()<3)) {
                        System.out.println("Warning: when resolving a tree new vertex has degree less than 3. Trouble ahead...");
                    }
                    theTree.remove(v);

                }
                catch (AlgorithmException anErr) {
                    System.out.println("Error during tree resolution. "+anErr.getMessage());
                }

            }
        }
        while (v!=null); // If you found an unresolved vertex, then loop thru' again

        // Set edge lengths to zero
        Iterator<Graph.Edge> itE = newEdges.iterator();
        for (int i=0; i<newEdges.size(); i++) {
            itE.next().setLength(0.0);
        }

        return newEdges;
    }


    /** Get NNI splits of a particular split.
     * NB: This assumes the tree is fully resolved, so that there are exactly two possibilities.
     The return array contains the two splits. */
    public static Split[] getNNISplits(TreeAsSplits theTreeMap, Split p) {
        if (theTreeMap.getNumSplits()!=(2*theTreeMap.getNumTaxa()-3)) {
            System.out.println("Warning: NNI operation requested on a tree which is not fully resolved. ");
        }

        // Sanity check: is p in the tree?
        if (!theTreeMap.contains(p)) {
            System.out.println("Algorithm error: bad split preparing an NNI move.");
        }
        // Comment out above later

        Split p1, p2;
        HashSet<String> pA1=null, pA2=null, pB1=null, pB2=null; // Sets of species on the 4 branches hanging from p
        HashSet<String> pA = p.getTaxonSubset(null); // Species on one side of p
        HashSet<String> pB = p.allTaxa(); // Species on the other side
        pB.removeAll(pA);
        HashSet<Split> splits = theTreeMap.getSplits();

        // Loop thru' p1 in the tree
        Iterator<Split> itA = splits.iterator();
        boolean done = false;
        for (int i=0; i<splits.size(); i++) {
            if (!done) {
                p1 = itA.next();
                if (!p.equals(p1)) {
                    pA1 = p1.getTaxonSubset(pA);
                    if (pA1!=null) {
                        pA2 = new HashSet();
                        pA2.addAll(pA);
                        pA2.removeAll(pA1);
                        // Is the split defined by pA2 in the tree?
                        HashSet<String> q = p.allTaxa();
                        q.removeAll(pA2);
                        try {
                            Split c = new Split(pA2,q);

                            if (theTreeMap.contains(c)) {
                                done = true; // Got pA1, pA2
                            }
                        }
                        catch (AlgorithmException anErr) {
                            System.out.println("Logic error making split for nni. "+anErr.getMessage());
                        }
                    }
                }
            }
        }
        done = false;
        Iterator<Split> itB = splits.iterator();
        for (int i=0; i<splits.size(); i++) {
            if (!done) {
                p1 = itB.next();
                if (!p.equals(p1)) {
                    pB1 = p1.getTaxonSubset(pB);
                    if (pB1!=null) {
                        pB2 = new HashSet();
                        pB2.addAll(pB);
                        pB2.removeAll(pB1);
                        // Is the split defined by pB2 in the tree?
                        HashSet<String> q = p.allTaxa();
                        q.removeAll(pB2);
                        try {
                            Split c = new Split(pB2,q);

                            if (theTreeMap.contains(c)) {
                                done = true; // Got pB1, pB2
                            }
                        }
                        catch (AlgorithmException anErr) {
                            System.out.println("Logic error making split for nni. "+anErr.getMessage());
                        }
                    }
                }
            }
        }

        // Sanity check: make sure union of hashSets in x covers the universe.
        if (DEBUG) {
            if ((pA1==null)|(pA2==null)|(pB1==null)|(pB2==null)) System.out.println("Bad NNI algorithm: null set of splits");
            HashSet<String> u = new HashSet();
            u.addAll(pA1);
            u.addAll(pA2);
            u.addAll(pB1);
            u.addAll(pB2);
            if (u.size()!=theTreeMap.getNumTaxa()) System.out.println("Bad NNI algorithm: universe not covered");
        }

        // Finally build the two splits
        Split q1=null, q2=null;
        try {
            HashSet<String> x = new HashSet();
            HashSet<String> y = new HashSet();
            x.addAll(pA1); x.addAll(pB1);
            y.addAll(pA2); y.addAll(pB2);
            q1 = new Split(x,y);
            x = new HashSet();
            y = new HashSet();
            x.addAll(pA1); x.addAll(pB2);
            y.addAll(pA2); y.addAll(pB1);
            q2 = new Split(x,y);
        }
        catch (AlgorithmException anErr) {
            System.out.println("Error creating NNI splits.");
        }

        // Another sanity check
        if (DEBUG) {
            Iterator<Split> it = theTreeMap.getSplitIterator();
            for (int i=0; i<theTreeMap.getNumSplits(); i++) {
                Split r = it.next();
                if (r!=p) {
                    if (!r.isCompatible(q1)) {
                        System.out.println("Bad sanity check finding nni splits");
                    }
                     if (!r.isCompatible(q2)) {
                        System.out.println("Bad sanity check finding nni splits");
                    }
                }
            }
        }

        // Return the splits
        Split[] res = new Split[2];
        res[0] = q1;
        res[1] = q2;
        return res;
    }
    
}