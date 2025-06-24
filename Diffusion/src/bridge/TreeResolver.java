/*
   TreeResolver
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
import java.util.HashSet;
import java.util.Iterator;

import simulation.CategoricalDistribution;
import simulation.Random;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Graph.Edge;
import treebase.Graph.Vertex;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;

/**
    
 */

public class TreeResolver {
    
    public static double resolveVertex(Tree theTree, Vertex unresolvedVertex, DoubleUniform unifDist, HashSet<Split> addedSplits) throws AlgorithmException {
        if (unresolvedVertex.degree()<4) return 0.0;
        
        if (unresolvedVertex.degree()==4) {
            /* Do something fast */
            // Sample a pair of vertices
            HashSet<Vertex> neighbours = unresolvedVertex.getNeighbours();
            Vertex w1 = (Vertex) CategoricalDistribution.sampleFromSetWithoutReplacement(neighbours, unifDist);
            Vertex w2 = (Vertex) CategoricalDistribution.sampleFromSetWithoutReplacement(neighbours, unifDist);
            Iterator<Vertex> itV = neighbours.iterator();
            Vertex u1 = itV.next();
            Vertex u2 = itV.next();
            
            try{
                // Cut edges
                Edge ew1 = w1.getEdge(unresolvedVertex); 
                double lw1 = ew1.getLength();
                theTree.cut(ew1);
                Edge ew2 = w2.getEdge(unresolvedVertex); 
                double lw2 = ew2.getLength();
                theTree.cut(ew2);
                Edge eu1 = u1.getEdge(unresolvedVertex); 
                double lu1 = eu1.getLength();
                theTree.cut(eu1);
                Edge eu2 = u2.getEdge(unresolvedVertex); 
                double lu2 = eu2.getLength();
                theTree.cut(eu2);
                theTree.remove(unresolvedVertex);

                // Join w1 & w2, u1 & u2
                Vertex w = theTree.addNewVertex("");
                Vertex u = theTree.addNewVertex("");
                theTree.connect(w1, w, lw1);
                theTree.connect(w2, w, lw2);
                theTree.connect(u1, u, lu1);
                theTree.connect(u2, u, lu2);
                Edge e = theTree.connect(u,w,0.0);
                if (addedSplits!=null) addedSplits.add(theTree.getSingleSplit(e));
            }
            catch (AlgorithmException anErr) {
                throw new AlgorithmError("Error resolving degree 4 vertex. "+anErr.getMessage());
            }
            
            /* return log 1/3 */
            return -1.0986123;
        }
        
        /* Sample unrooted topology uniformly at random */
        int n = unresolvedVertex.degree();
        HashSet<Vertex> stems = unresolvedVertex.getNeighbours();
        Vertex r = theTree.addNewVertex("");
        HashSet<Edge> edges = new HashSet();
        /* Join 3 stems to r */
        for (int i=1; i<=3; i++) {
            Vertex u = (Vertex) CategoricalDistribution.sampleFromSetWithoutReplacement(stems, unifDist);
            Edge e = theTree.connect(u, r, 0.0);
            edges.add(e);
        }
        
        /* Repeat n-3 times:
            Select a stem at random
            Subdivide an edge
            Join stem to edge
        */
        for (int i=0; i<(n-3); i++) {
            Vertex u = (Vertex) CategoricalDistribution.sampleFromSetWithoutReplacement(stems, unifDist);
            Edge e = (Edge) CategoricalDistribution.sampleFromSetWithoutReplacement(edges, unifDist);
            // Divide e
            Vertex deg2v = theTree.divide(e, "");
            Iterator<Edge> itE = deg2v.edgeIterator();
            edges.add(itE.next()); // Add edge 1
            edges.add(itE.next()); // Add edge 2
            // Join stem to new vertex
            Edge f = theTree.connect(u, deg2v, 0.0);
            edges.add(f);
        }
        
        /* Sort out stem edge lengths and remove unresolved vertex. */
        int count = 0;
        for (Edge e : edges) {
            Vertex[] w = e.getVertices();
            if ((w[0].isConnectedTo(unresolvedVertex))&&(w[1].isConnectedTo(unresolvedVertex))) {
                throw new AlgorithmError("Bad graph logic resolving vertex.");
            }
            if (w[0].isConnectedTo(unresolvedVertex)) {
                Edge eu = w[0].getEdge(unresolvedVertex);
                e.setLength(eu.getLength());
                theTree.cut(eu);
                count++;
            }
            if (w[1].isConnectedTo(unresolvedVertex)) {
                Edge eu = w[1].getEdge(unresolvedVertex);
                e.setLength(eu.getLength());
                theTree.cut(eu);
                count++;
            }
        }
        theTree.remove(unresolvedVertex);
        if (count!=n) {
            throw new AlgorithmError("Failed sanity check resolving vertex.");
        }
        
        /* Get the splits for the new edges */
        if (addedSplits!=null) {
            theTree.generateSplits();
            for (Edge e : edges) {
                addedSplits.add(theTree.getSplit(e));
            }
        }
        
        /* Compute log density: -log (2n-5)!! */
        double logDensity = 0.0;
        for (int i=1; i<=(n-2); i++) {
            logDensity -= Math.log(2*i-1);
        }
        
        return logDensity;
    }
    
    public static double resolveTree(Tree theTree, DoubleUniform unifDist) throws AlgorithmException {
        
        double logDens = 0.0;
        Vertex unresolved;
        do {
            unresolved = null;
            Iterator<Vertex> it = theTree.getVertexIterator();
            while (it.hasNext()) {
                Vertex v = it.next();
                if (v.degree()>3) {
                    unresolved = v;
                }
            }
            if (unresolved!=null) logDens += resolveVertex(theTree, unresolved, unifDist, null);
        }
        while (unresolved!=null);
        
        return logDens;
    }
 
    public static double resolveTree(TreeAsSplits t, DoubleUniform unifDist) throws AlgorithmException {
        
        if (t.fullyResolved()) return 0.0;
        
        Tree theTree = null;
        try {
            theTree = t.getTree();
        } catch (AlgorithmError ex) {
            System.out.println("Error converting splts to tree when resolving.");
        }
        
        double logDens = 0.0;
        Vertex unresolved;
        HashSet<Split> splits = new HashSet();
        do {
            unresolved = null;
            Iterator<Vertex> it = theTree.getVertexIterator();
            while (it.hasNext()) {
                Vertex v = it.next();
                if (v.degree()>3) {
                    unresolved = v;
                }
            }
            if (unresolved!=null) logDens += resolveVertex(theTree, unresolved, unifDist, splits);
        }
        while (unresolved!=null);
        
        /* Add splits */
        for (Split s : splits) {
            if (!t.contains(s)) {
                t.add(s, 0.0);
            }
        }
        if (!t.fullyResolved()) {
            throw new AlgorithmError("Sanity check failed resolving tree: tree is unresolved!");
        }
        
        return logDens;
    }
    
    
    /** Pass in a treeAsSplits object plus the Tree it is equivalent to.  */
    public static double resolveTree(TreeAsSplits t, Tree theTree, DoubleUniform unifDist) throws AlgorithmException {
        
        if (t.fullyResolved()) return 0.0;
        
        double logDens = 0.0;
        Vertex unresolved;
        HashSet<Split> splits = new HashSet();
        do {
            unresolved = null;
            Iterator<Vertex> it = theTree.getVertexIterator();
            while (it.hasNext()) {
                Vertex v = it.next();
                if (v.degree()>3) {
                    unresolved = v;
                }
            }
            if (unresolved!=null) logDens += resolveVertex(theTree, unresolved, unifDist, splits);
        }
        while (unresolved!=null);
        
        /* Add splits */
        for (Split s : splits) {
            if (!t.contains(s)) {
                t.add(s, 0.0);
            }
        }
        if (!t.fullyResolved()) {
            throw new AlgorithmError("Sanity check failed resolving tree: tree is unresolved!");
        }
        
        return logDens;
    }
    
    
    /* TESTING AREA */

    public static void main(String[] args) throws AlgorithmException {        
        DoubleUniform unifDist = new DoubleUniform(Random.getEngine());
        
        Tree x = new Tree("(A:1,B:1,C:1,D:1,E:1,(F:1,G:1):1);");
        Iterator<Vertex> it = x.getVertexIterator();
        Vertex v = null;
        while (it.hasNext()) {
            v = it.next();
            if (v.degree()!=1) break;
        }
        double l = resolveVertex(x, v, unifDist, null);
        System.out.println("\n");
        System.out.println(x.toString());
        System.out.println(String.format("%7.7f",l));
        
        System.out.println("\n");
        Tree y = new Tree("(A:1,B:1,C:1,D:1,E:1,(F:1,G:1,H:1):1);");
        l = resolveTree(y, unifDist);
        System.out.println(y.toString());
        System.out.println(String.format("%7.7f",l));
    }

}
