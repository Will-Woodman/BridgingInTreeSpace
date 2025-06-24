/*
 * TreeWithTopologicalOperations
   Copyright (C) 2013  Sarah E. Heaps and Tom M. W.Nye

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

    Contact the authors at:  <sarah.heaps@ncl.ac.uk>, <tom.nye@ncl.ac.uk>
 */

package treebase;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Extends RootedTree to incorporate topological operations
 * NB: always recalculate splits after any topological operation!
 */

public class RootedTreeWithTopologicalOperations extends RootedTree implements TopologicalOperations {

    /* Constructors  */

    /** No arg constructor for use by extending classes */
    protected RootedTreeWithTopologicalOperations() {
        super();
    }

    /** Constructors based on Newick strings */
    public RootedTreeWithTopologicalOperations(String theString) throws AlgorithmException {
        super(theString);
    }

    public RootedTreeWithTopologicalOperations(java.io.File theFile) throws AlgorithmException {
        super(theFile);
    }

    /** Constructors based on rooting an unrooted tree */
    public RootedTreeWithTopologicalOperations(Tree t, Vertex r, HashMap<Vertex,Vertex> corresp) throws AlgorithmError {
        super(t, r, corresp);
    }
    public RootedTreeWithTopologicalOperations(Tree t, Vertex r) throws AlgorithmError {
        super(t, r);
     }
    public RootedTreeWithTopologicalOperations(Tree t, Edge e, HashMap<Vertex,Vertex> corresp) throws AlgorithmError {
        super(t, e, corresp);
    }
    public RootedTreeWithTopologicalOperations(Tree t, Edge e) throws AlgorithmError {
        super(t, e);
     }

    /** Constructor from a hashmap from splits to lengths */
    public RootedTreeWithTopologicalOperations(HashMap<Split,Double> splitsToLengths, Split rootSplit) throws AlgorithmError {
        super(splitsToLengths, rootSplit);
    }

    /** Constructor from a RootedTree  */
    public RootedTreeWithTopologicalOperations(RootedTree t) {
        super();
        this.buildFromTemplate(t);
    }

    /** Clone  */
    public RootedTreeWithTopologicalOperations clone() {
        return this.clone(null);
    }

    /** Clone a tree but also return the correspondence from old to new vertices.
     Pass in null if this information is not needed, otherwise an empty map */
    public RootedTreeWithTopologicalOperations clone(HashMap<Vertex,Vertex> corresp) {
        RootedTreeWithTopologicalOperations t = new RootedTreeWithTopologicalOperations();
        t.buildFromTemplate(this, corresp);
        return t;
    }


    /* Topological moves */

    /** Perform an NNI on an edge.
     NB: THIS ASSUMES FULLY RESOLVED TREE!
     Pass in the array of vertices obtained from getVerticesAdjacentToEdge.
     This array is modified by the method so it can be passed back to reverseNNI.
    */
    @Override
    public Edge performNNI(Edge e, Vertex[][] v, double newLength, boolean choice) throws AlgorithmException {
        // Check whether edge is connected to the root
        int rootInd = 0;
        if(e.vertex1==getRoot()) rootInd = 1;
        if(e.vertex2==getRoot()) rootInd = 2;

        Vertex top=null;
        if (rootInd==0) {
            // Set up info for ancestry calcs
            Vertex anc1 = ((VertexWithAncestor)e.vertex1).getAncestor();
            if (anc1!=e.vertex2) {
                top = anc1;
            }
            else {
                Vertex anc2 = ((VertexWithAncestor)e.vertex2).getAncestor();
                if (anc2!=e.vertex1) {
                    top = anc2;
                }
                else throw new AlgorithmError("Error with ancestry calcs following NNI.");
            }
        }

        if (!edges.contains(e)) throw new AlgorithmError("Bad edge sent to request to perform NNI.");

        // If the vertex at either end of e is of degree two, throw exception
        if(e.vertex1.degree()==2 || e.vertex2.degree()==2) throw new AlgorithmException("Cannot perform NNI on edge connected to a degree 2 vertex.");

        // Cut the edge
        this.cut(e);
        // Add two new vertices
        Vertex w0 = this.addNewVertex("");
        Vertex w1 = this.addNewVertex("");
        Edge newEdge = this.connect(w0, w1, newLength);

        /* Now connect up the vertices */
        Edge f;
        if (choice) {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w1);
        }
        else {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w1);
        }
        this.remove(v[0][0]);
        this.remove(v[1][0]);
        v[0][0] = w0;
        v[1][0] = w1;
        if (choice) {
            Graph.Vertex w = v[1][1];
            v[1][1] = v[0][2];
            v[0][2] = w;
        } else {
            Graph.Vertex w = v[1][2];
            v[1][2] = v[0][2];
            v[0][2] = w;
        }

        // If NNI was on an edge connected to a degree 3 root, set root to vertex on new edge
        if(rootInd==1) setRoot((VertexWithAncestor)v[0][0]);
        else if(rootInd==2) setRoot((VertexWithAncestor)v[1][0]);

        // Recalculate ancestry
        if(rootInd==0) {

            /* Small local change
             Four cases for the subtree containing the root, and for each case five ancestors to set */
            if (top==v[0][1]) {
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[0][1]);
                ((VertexWithAncestor)v[0][2]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][1]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[1][2]).setAncestor((VertexWithAncestor)v[1][0]);
            }
            else if (top==v[0][2]) {
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[0][2]);
                ((VertexWithAncestor)v[0][1]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][1]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[1][2]).setAncestor((VertexWithAncestor)v[1][0]);
            }
            else if (top==v[1][1]) {
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[1][1]);
                ((VertexWithAncestor)v[1][2]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][1]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[0][2]).setAncestor((VertexWithAncestor)v[0][0]);
            }
            else if (top==v[1][2]) {
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[1][2]);
                ((VertexWithAncestor)v[1][1]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][1]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[0][2]).setAncestor((VertexWithAncestor)v[0][0]);
            }
            else throw new AlgorithmError("Error with ancestry calcs following NNI.");


        } else {
            recalculateAncestry();
        }

        return newEdge;
    }


    /** Reverse the NNI on an edge.
     NB: THIS ASSUMES FULLY RESOLVED TREE!
     The oldEdge and NewEdge should be the 1st and 2nd entries in the array returned by performNNI. The vertex[][] should
     be that passed to performNNI (and modified by that method). The boolean should be the same as that passed to
     performNNI in order to reverse the move.
     */
    @Override
    public Edge reverseNNI(Edge oldEdge, Edge newEdge, Vertex[][] v, boolean choice) throws AlgorithmException {
        if (!edges.contains(newEdge)) throw new AlgorithmError("Bad edge sent to request to reverse NNI.");

        // If the vertex at either end of e is of degree two, throw exception
        if(newEdge.vertex1.degree()==2 || newEdge.vertex2.degree()==2) throw new AlgorithmException("Cannot reverse NNI on edge connected to a degree 2 vertex. This should not be possible.");

        // Check whether edge is connected to the root
        int rootInd = 0;
        if(newEdge.vertex1==getRoot()) rootInd = 1;
        if(newEdge.vertex2==getRoot()) rootInd = 2;

        Vertex top=null;
        if (rootInd==0) {
            // Set up info for ancestry calcs
            Vertex anc1 = ((VertexWithAncestor)newEdge.vertex1).getAncestor();
            if (anc1!=newEdge.vertex2) {
                top = anc1;
            }
            else {
                Vertex anc2 = ((VertexWithAncestor)newEdge.vertex2).getAncestor();
                if (anc2!=newEdge.vertex1) {
                    top = anc2;
                }
                else throw new AlgorithmError("Error with ancestry calcs reversing NNI.");
            }
        }

        // Cut the edge
        this.cut(newEdge);
        // Add two new vertices
        Vertex w0 = this.addNewVertex(""); //when I come to override, addNewVertex can contain the sequences which
        Vertex w1 = this.addNewVertex(""); //I'll be able to get from oldEdge
        Edge oldEdgeCopy = this.connect(w0, w1, oldEdge.getLength()); //when I come to override, connect can contain
                                                   //the subs history which I can get from oldEdge

        // Now connect up the vertices
        Edge f;
        if (choice) {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w1);
        }
        else {
            f = v[0][0].getEdge(v[0][1]); //fine
            f.changeVertex(v[0][0], w0); //fine
            f = v[1][0].getEdge(v[1][2]); //breaks
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w1);
        }
        this.remove(v[0][0]);
        this.remove(v[1][0]);
        v[0][0] = w0;
        v[1][0] = w1;
        if (choice) {
            Graph.Vertex w = v[1][1];
            v[1][1] = v[0][2];
            v[0][2] = w;
        } else {
            Graph.Vertex w = v[1][2];
            v[1][2] = v[0][2];
            v[0][2] = w;
        }

        // If NNI was on an edge connected to a degree 3 root, set root to vertex on new edge
        if(rootInd==1) setRoot((VertexWithAncestor)w0);
        else if(rootInd==2) setRoot((VertexWithAncestor)w1);

        // Recalculate ancestry
        if(rootInd==0) {

            /* Small local change
             Four cases for the subtree containing the root, and for each case five ancestors to set */
            if (top==v[0][1]) {
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[0][1]);
                ((VertexWithAncestor)v[0][2]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][1]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[1][2]).setAncestor((VertexWithAncestor)v[1][0]);
            }
            else if (top==v[0][2]) {
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[0][2]);
                ((VertexWithAncestor)v[0][1]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[1][1]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[1][2]).setAncestor((VertexWithAncestor)v[1][0]);
            }
            else if (top==v[1][1]) {
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[1][1]);
                ((VertexWithAncestor)v[1][2]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][1]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[0][2]).setAncestor((VertexWithAncestor)v[0][0]);
            }
            else if (top==v[1][2]) {
                ((VertexWithAncestor)v[1][0]).setAncestor((VertexWithAncestor)v[1][2]);
                ((VertexWithAncestor)v[1][1]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][0]).setAncestor((VertexWithAncestor)v[1][0]);
                ((VertexWithAncestor)v[0][1]).setAncestor((VertexWithAncestor)v[0][0]);
                ((VertexWithAncestor)v[0][2]).setAncestor((VertexWithAncestor)v[0][0]);
            }
            else throw new AlgorithmError("Error with ancestry calcs following NNI.");

        } else {
            recalculateAncestry();
        }

        return oldEdgeCopy;
    }


    /** Carry out an SPR move: override unrooted version */
    public void performTopologicalSPR(Edge pruneEdge, Edge graftEdge) {

        try {
            cut(pruneEdge);

            /* Work out which vertex on the prune edge attaches to the "main tree"
             (ie. the component of the graph containing the graft edge)
             and which vertex is on the pruned sub-tree.
             */

            Vertex[] v = pruneEdge.getVertices();
            Vertex subtreeVertex = null;
            Vertex maintreeVertex = null;
            if (this.subtreeContainsEdge(v[0], v[1], graftEdge)) {
                subtreeVertex = v[1];
                maintreeVertex = v[0];
            }
            else if (this.subtreeContainsEdge(v[1], v[0], graftEdge)) {
                subtreeVertex = v[0];
                maintreeVertex = v[1];
            }
            else {
                javax.swing.JOptionPane.showMessageDialog(null,"Serious problem performing SPR: subtree not found.","Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
                System.exit(1);
            }

            // OK got the vertices

            // Subdivide the graft edge
            Vertex attachVertex = this.divide(graftEdge, "");
            // Connect the new vertex to the pruned tree
            connect(subtreeVertex, attachVertex, pruneEdge.getLength());
            // Absorb degree two vertices if necessary
            if (maintreeVertex.degree()==2) {
                removeDegreeTwoVertex(maintreeVertex);
            }

        }
        catch (AlgorithmException anErr) {
                javax.swing.JOptionPane.showMessageDialog(null,"Unable to perform SPR: "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }
        recalculateAncestry();

    }

    /** Carry out a reversible SPR move on a fully resolved tree. You can specify new lengths
        for the four new edges created during the move. Note that regrafting the subtree on
        edge E creates two new edges on either side of the attached subtree. The edge which defines
        the same split as the original edge, E, is given length lengthSplitEdge1. The other is given
        the length lengthSplitEdge2. NOTE: slight inconsistency in the case that the prune and graft
        edges are adjacent. The four new branch lengths are not well defined: lengthMergedEdge and
        lengthSplitEdge2 refer to the same edge. The former length is used.
        The method return an array of edges of dimension 4 x 2. Components 0-3 in the first column contain
        edges from the original tree: the two edges merged into one where the subtree is plucked
        (entry 0 is the edge defining the same split as the merged edge in the new tree, entry 1 is the
        other); the edge which is split in to two (due to regrafting) and the edge which is pruned.
        Components 0-3 in the second column contain edges from the new tree: the new edges on either side
        of the newly attached subtree (in the order above), and the graft edge and the prune edges for the
        reverse move.*/
    @Override
    public Edge[][] performSPR(Edge pruneEdge, Edge graftEdge, double lengthMergedEdge, double lengthSplitEdge1, double lengthSplitEdge2, double lengthPrunedEdge) throws AlgorithmError {
        boolean recalcAncestry = false; // For a rooted tree with degree 3 root, does the root "move"?
        boolean graftDescendantOfPrune; // Is the graft edge a descendant of the prune edge?
        Edge[][] edgeArray = new Edge[4][2];
        // Record the edges which are plucked off and split in to two
        edgeArray[2][0] = graftEdge;
        edgeArray[3][0] = pruneEdge;
        // Are the edges adjacent?
        boolean areAdjacent = pruneEdge.getNeighbouringEdges().contains(graftEdge);
        try {
            // Work out which vertex on the prune edge attaches to the "main tree"
            // (ie. the component of the graph containing the graft edge)
            // and which vertex is on the pruned sub-tree.
            VertexWithAncestor subtreeVertex=null, maintreeVertex=null;
            VertexWithAncestor graftUp = getAncestralVertex(graftEdge);
            VertexWithAncestor graftDown = getDescendantVertex(graftEdge);

            Vertex[] v = pruneEdge.getVertices();
            if (this.subtreeContainsEdge(v[0], v[1], graftEdge)) {
                subtreeVertex = (VertexWithAncestor)v[1];
                maintreeVertex = (VertexWithAncestor)v[0];
            }
            else if (this.subtreeContainsEdge(v[1], v[0], graftEdge)) {
                subtreeVertex = (VertexWithAncestor)v[0];
                maintreeVertex = (VertexWithAncestor)v[1];
            }
            else {
                throw new AlgorithmError("Serious problem performing SPR: subtree not found.");
            }
            if((subtreeVertex==getRoot() || maintreeVertex==getRoot()) && getRoot().degree()==2) throw new AlgorithmException("Cannot prune from "
                    + "edge connected to the root.");
            else {
                if(maintreeVertex.getAncestor()==subtreeVertex) graftDescendantOfPrune = true;
                else if(subtreeVertex.getAncestor()==maintreeVertex) graftDescendantOfPrune = false;
                else throw new AlgorithmException("Problem finding direction of ancestry on prune edge. This should "
                        + "not be possible.");
            }

            Iterator<Edge> itE = maintreeVertex.getEdges().iterator();
            Edge eM1 = itE.next();
            Edge eM3 = itE.next();
            if(eM1==pruneEdge) eM1 = itE.next();
            if(eM3==pruneEdge) eM3 = itE.next();
            // Find splits on the edges which will be merged and the split defined by the grafting egde
            Split originalSplitM1 = this.getSingleSplit(eM1);
            Split originalSplitM3 = this.getSingleSplit(eM3);
            Split originalSplitS = this.getSingleSplit(graftEdge);

            // Pluck off the subtree
            cut(pruneEdge);

            // OK got the vertices

            // Subdivide the graft edge
            VertexWithAncestor attachVertex = (VertexWithAncestor)this.divide(graftEdge, "");
            itE = attachVertex.getEdges().iterator();
            Edge eS1 = itE.next();
            Edge eS3 = itE.next();

            // Connect the new vertex to the pruned tree
            connect(subtreeVertex, attachVertex, lengthPrunedEdge);

            // Record the prune edge for the reverse move:
            edgeArray[3][1] = attachVertex.getEdge(subtreeVertex);

            // Absorb degree two vertex
            if(maintreeVertex==getRoot()) {
                setRoot(attachVertex);
                recalcAncestry = true;
            }
            Iterator<Vertex> it = maintreeVertex.neighboursIterator();
            VertexWithAncestor v1 = (VertexWithAncestor)it.next();
            VertexWithAncestor v3 = (VertexWithAncestor)it.next();
            // Remove the degree two vertex and recale the resulting merged edge
            removeDegreeTwoVertex(maintreeVertex);
            v1.getEdge(v3).setLength(lengthMergedEdge);
            // Record the graft edge for the reverse move
            edgeArray[2][1] = v1.getEdge(v3);

            // Now that the new tree has been constructed, find the split on new
            // merged edge and the splits on edges eS1 and eS3
            Split newSplitS1, newSplitS3;
            if(areAdjacent) {
                if(edges.contains(eS1)) newSplitS1 = this.getSingleSplit(eS1);
                else newSplitS1 = this.getSingleSplit(v1.getEdge(v3));
                if(edges.contains(eS3)) newSplitS3 = this.getSingleSplit(eS3);
                else newSplitS3 = this.getSingleSplit(v1.getEdge(v3));
            } else {
                newSplitS1 = this.getSingleSplit(eS1);
                newSplitS3 = this.getSingleSplit(eS3);
            }
            if(originalSplitS.equals(newSplitS1)) {
                // If true, the edge eS1 defines the same split as the original
                // edge prior to being split in to two. Otherwise, this should
                // be true of the edge eS3
                eS1.setLength(lengthSplitEdge1);
                eS3.setLength(lengthSplitEdge2);
                edgeArray[0][1] = eS1;
                edgeArray[1][1] = eS3;
            } else if(originalSplitS.equals(newSplitS3)) {
                eS1.setLength(lengthSplitEdge2);
                eS3.setLength(lengthSplitEdge1); //lengthSplitEdge1 applies to the edge defining the same split
                edgeArray[0][1] = eS3;
                edgeArray[1][1] = eS1;
            } else throw new AlgorithmError("Problem performing SPR: matching splits not found. This should not be possible.");

            Split newSplitM = this.getSingleSplit(v1.getEdge(v3));
            if(originalSplitM1.equals(newSplitM)) {
                // If true, the merged edge on the new tree defines the same
                // split as the original edge eM1. Otherwise,  this is true of
                // the original edge eM3.
                edgeArray[0][0] = eM1; //edgeArray[0][0] applies to the edge defining the same split
                edgeArray[1][0] = eM3;
            } else if(originalSplitM3.equals(newSplitM)) {
                edgeArray[0][0] = eM3;
                edgeArray[1][0] = eM1;
            } else throw new AlgorithmException("Problem performing SPR: matching splits not found. This should not be possible.");

            // Recalculate ancestry
            if(recalcAncestry || areAdjacent) this.recalculateAncestry();
            else {
                if(!graftDescendantOfPrune) {
                    subtreeVertex.setAncestor(attachVertex);
                    attachVertex.setAncestor(graftUp);
                    graftDown.setAncestor(attachVertex);
                } else {
                    attachVertex.setAncestor(subtreeVertex);
                    graftDown.setAncestor(attachVertex);
                    // Reverse ancestry on the spine of the tree between graftUp
                    // and the vertex on the end of the edge created by merging
                    VertexWithAncestor u = reverseAncestryInSPR(graftUp, graftUp.getAncestor());
                    Iterator<Vertex> itN = u.neighboursIterator();
                    for(int i=0; i<u.degree(); i++) {
                        VertexWithAncestor x = (VertexWithAncestor) itN.next();
                        if(x.getAncestor()==null) {
                            x.setAncestor(u);
                            break;
                        }
                    }
                    graftUp.setAncestor(attachVertex);
                    //recalculateAncestry();
                }

            }

        }
        catch (AlgorithmException anErr) {
                throw new AlgorithmError("Algorithm Error: Unable to perform SPR: "+anErr.getMessage());
        }

        return edgeArray;

    }


    protected VertexWithAncestor reverseAncestryInSPR(VertexWithAncestor v, VertexWithAncestor fromV) throws AlgorithmException {

        if(fromV==null) return v;
        VertexWithAncestor w = reverseAncestryInSPR(fromV, fromV.getAncestor());
        fromV.setAncestor(v);
        return w;
    }


    /** Reverses the SPR move performSPR. Pass in the array of edges returned by this method. The array returned here
        contains the four edges created during the move.*/
    @Override
    public Edge[] reverseSPR(Edge[][] edgeArray) throws AlgorithmError {
        boolean recalcAncestry = false; // For a rooted tree with degree 3 root, does the root "move"?
        boolean graftDescendantOfPrune; // Is the graft edge a descendant of the prune edge?

        Edge mergedEdge1 = edgeArray[0][0];
        Edge mergedEdge2 = edgeArray[1][0];
        Edge splitEdge = edgeArray[2][0];
        Edge prunedEdge = edgeArray[3][0];
        Edge reverseGraftEdge = edgeArray[2][1];
        Edge reversePruneEdge = edgeArray[3][1];

        boolean areAdjacent = reversePruneEdge.getNeighbouringEdges().contains(reverseGraftEdge);

        Edge[] originalEdges = new Edge[4];

        try {
            // Work out which vertex on the prune edge attaches to the "main tree"
            // (ie. the component of the graph containing the graft edge)
            // and which vertex is on the pruned sub-tree.
            // Simple on a rooted tree!
            VertexWithAncestor subtreeVertex=null, maintreeVertex=null;
            VertexWithAncestor graftUp = getAncestralVertex(reverseGraftEdge);
            VertexWithAncestor graftDown = getDescendantVertex(reverseGraftEdge);

            Vertex[] v = reversePruneEdge.getVertices();
            if (this.subtreeContainsEdge(v[0], v[1], reverseGraftEdge)) {
                subtreeVertex = (VertexWithAncestor)v[1];
                maintreeVertex = (VertexWithAncestor)v[0];
            }
            else if (this.subtreeContainsEdge(v[1], v[0], reverseGraftEdge)) {
                subtreeVertex = (VertexWithAncestor)v[0];
                maintreeVertex = (VertexWithAncestor)v[1];
            }
            else {
                throw new AlgorithmError("Serious problem performing SPR: subtree not found.");
            }
            if((subtreeVertex==getRoot() || maintreeVertex==getRoot()) && getRoot().degree()==2) throw new AlgorithmException("Cannot prune from "
                    + "edge connected to the root.");
            else {
                if(maintreeVertex.getAncestor()==subtreeVertex) graftDescendantOfPrune = true;
                else if(subtreeVertex.getAncestor()==maintreeVertex) graftDescendantOfPrune = false;
                else throw new AlgorithmException("Problem finding direction of ancestry on prune edge. This should "
                        + "not be possible.");
            }

            // Find split defined by the grafting egde
            Split originalSplitS = this.getSingleSplit(reverseGraftEdge);

            // Pluck off the subtree
            cut(reversePruneEdge);

            // OK got the vertices

            // Subdivide the graft edge
            VertexWithAncestor attachVertex = (VertexWithAncestor)this.divide(reverseGraftEdge, "");
            Iterator<Edge> itE = attachVertex.getEdges().iterator();
            Edge eS1 = itE.next();
            Edge eS3 = itE.next();

            // Connect the new vertex to the pruned tree
            originalEdges[3] = connect(subtreeVertex, attachVertex, prunedEdge.getLength());

            // Now the pruned subtree is attached, find which of the two edges
            // eS1 and eS3 defines the same split as the original edge
            Split newSplitS1 = this.getSingleSplit(eS1);
            Split newSplitS3 = this.getSingleSplit(eS3);
            if(originalSplitS.equals(newSplitS1)) {
                // If true, the edge eS1 defines the same split as the original
                // edge. Otherwise, this should true of the edge eS3
                eS1.setLength(mergedEdge1.getLength());
                eS3.setLength(mergedEdge2.getLength());
                originalEdges[0] = eS1;
                originalEdges[1] = eS3;
            } else if(originalSplitS.equals(newSplitS3)) {
                eS1.setLength(mergedEdge2.getLength());
                eS3.setLength(mergedEdge1.getLength()); //mergedEdge1 applies to the edge defining the same split
                originalEdges[0] = eS3;
                originalEdges[1] = eS1;
            } else throw new AlgorithmError("Problem performing SPR: matching splits not found. This should not be possible.");

            // Absorb degree two vertices if necessary
            if(maintreeVertex==getRoot()) {
                setRoot(attachVertex);
                recalcAncestry = true;
            }
            Iterator<Vertex> it = maintreeVertex.neighboursIterator();
            VertexWithAncestor v1 = (VertexWithAncestor)it.next();
            VertexWithAncestor v3 = (VertexWithAncestor)it.next();
            // Remove the degree two vertex and recale the resulting merged edge
            removeDegreeTwoVertex(maintreeVertex);
            originalEdges[2] = v1.getEdge(v3);
            originalEdges[2].setLength(splitEdge.getLength());

            // Recalculate ancestry
            if(recalcAncestry || areAdjacent) this.recalculateAncestry();
            else {
                if(!graftDescendantOfPrune) {
                    subtreeVertex.setAncestor(attachVertex);
                    attachVertex.setAncestor(graftUp);
                    graftDown.setAncestor(attachVertex);
                } else {
                    attachVertex.setAncestor(subtreeVertex);
                    graftDown.setAncestor(attachVertex);
                    // Reverse ancestry on the spine of the tree between graftUp
                    // and the vertex on the end of the edge created by merging
                    VertexWithAncestor u = reverseAncestryInSPR(graftUp, graftUp.getAncestor());
                    Iterator<Vertex> itN = u.neighboursIterator();
                    for(int i=0; i<u.degree(); i++) {
                        VertexWithAncestor x = (VertexWithAncestor) itN.next();
                        if(x.getAncestor()==null) {
                            x.setAncestor(u);
                            break;
                        }
                    }
                    graftUp.setAncestor(attachVertex);
                    //recalculateAncestry();
                }

            }

        }
        catch (AlgorithmException anErr) {
                throw new AlgorithmError("Algorithm Error: Unable to reverse SPR: "+anErr.getMessage());
        }

        return originalEdges;
    }


    /** Move the (degree two) root. You can specify new lengths for the three new edges created during
        the move. Note that removing a degree two root creates a new single edge whose length is given
        by lengthMergedEdge. Inserting a degree two root on the newRootEdge creates two new edges on
        either side. Of these two edges, the one which is ancestral to the new merged edge is assigned
        the length lengthSplitEdge1. The other is assigned the length lengthSplitEdge2.
        The method returns an array of edges of dimension 3 x 2. In the first column, components 0-2 contain
        edges from the original tree: the edges on either side of the old root, ordered as above in
        relation to the edge chosen for rooting and the new root edge before it is divided into two.
        In the second column, components 0-2 contain edges from the new tree, namely the edges on either
        side of the root assigned lengths lengthSplitEdge1 and lengthSplitEdge1 respectively, and the
        newRootEdge for the reverse move (ie the edge created by merging during the root move).*/
    public Edge[][] performRootMove(Edge newRootEdge, double lengthMergedEdge, double lengthSplitEdge1, double lengthSplitEdge2) throws AlgorithmError {
        if(getRoot().degree()!=2) throw new AlgorithmError("Algorithm Error: root moves are only intended for rooted trees whose root "
                + "is of degree 2.");
        Edge[][] edgeArray = new Edge[3][2];
        try {
            edgeArray[2][0] = newRootEdge;

            Iterator<Edge> itE = getRoot().edgeIterator();
            Edge eA = itE.next();
            HashSet<Edge> eChildren = new HashSet<Edge>();
            this.getAllDescendants(eA, null, eChildren);
            if(eChildren.contains(newRootEdge)) {
                edgeArray[0][0] = eA;
                edgeArray[1][0] = itE.next();
            } else {
                Edge eB = itE.next();
                eChildren.clear();
                this.getAllDescendants(eB, null, eChildren);
                if(eChildren.contains(newRootEdge)) {
                    edgeArray[0][0] = eB;
                    edgeArray[1][0] = eA;
                } else throw new AlgorithmException("Problem labelling edges during root move.");
            }

            // Remove degree two root and rescale newly created edge
            Iterator<Vertex> itN = getRoot().neighboursIterator();
            VertexWithAncestor v1 = (VertexWithAncestor)itN.next(); VertexWithAncestor v2 = (VertexWithAncestor)itN.next();
            cut(edgeArray[0][0]); cut(edgeArray[1][0]);
            vertices.remove(getRoot());
            edgeArray[2][1] = connect(v1, v2, lengthMergedEdge);

            // Subdivide the new root edge and label the new vertex as the root
            VertexWithAncestor newRoot = (VertexWithAncestor)this.divide(newRootEdge, "");
            this.setRoot(newRoot);

            // Recalculate ancestry
            this.recalculateAncestry();

            // Assign labels and set lengths accordingly
            itE = getRoot().edgeIterator();
            Edge eeA = itE.next();
            eChildren.clear();
            this.getAllDescendants(eeA, null, eChildren);
            if(eChildren.contains(edgeArray[2][1])) {
                eeA.setLength(lengthSplitEdge1);
                edgeArray[0][1] = eeA;
                edgeArray[1][1] = itE.next();
                edgeArray[1][1].setLength(lengthSplitEdge2);
            } else {
                Edge eeB = itE.next();
                eChildren.clear();
                this.getAllDescendants(eeB, null, eChildren);
                if(eChildren.contains(edgeArray[2][1])) {
                    eeB.setLength(lengthSplitEdge1);
                    eeA.setLength(lengthSplitEdge2);
                    edgeArray[0][1] = eeB;
                    edgeArray[1][1] = eeA;
                } else throw new AlgorithmException("Problem labelling edges during root move.");
            }

        } catch (AlgorithmException anErr) {
            edgeArray = null;
            throw new AlgorithmError("Algorithm Error: Unable to perform root move: "+anErr.getMessage());
        }

        return edgeArray;
    }


    /** Reverses the move-root method. Pass the array returned by performRootMove into
       reverseRootMove. The edges returned by this method are the three edges created during
       the move.*/
    public Edge[] reverseRootMove(Edge[][] edgeArray) throws AlgorithmError {
        if(getRoot().degree()!=2) throw new AlgorithmError("Algorithm Error: root moves are only intended for rooted trees whose root "
                + "is of degree 2.");
        Edge mergedEdge1 = edgeArray[0][0];
        Edge mergedEdge2 = edgeArray[1][0];
        Edge rootEdge = edgeArray[2][0];
        Edge reverseRootEdge = edgeArray[2][1];

        Edge[] originalEdges = new Edge[3];

        try {
            // Remove degree two root and rescale newly created edge
            Iterator<Vertex> itN = getRoot().neighboursIterator();
            VertexWithAncestor v1 = (VertexWithAncestor)itN.next(); VertexWithAncestor v2 = (VertexWithAncestor)itN.next();
            Edge e1 = getRoot().getEdge(v1); Edge e2 = getRoot().getEdge(v2);
            cut(e1); cut(e2);
            vertices.remove(getRoot());
            originalEdges[2] = connect(v1, v2, rootEdge.getLength());

            // Subdivide the new root edge and label the new vertex as the root
            Vertex[] vcs1 = mergedEdge1.getVertices();
            VertexWithAncestor newRoot = (VertexWithAncestor)this.divide(reverseRootEdge, "");
            this.setRoot(newRoot);

            // Recalculate ancestry
            this.recalculateAncestry();

            // Assign labels and set lengths accordingly
            Iterator<Edge> itE = getRoot().edgeIterator();
            Edge eeA = itE.next();
            HashSet<Edge> eChildren = new HashSet<Edge>();
            this.getAllDescendants(eeA, null, eChildren);
            if(eChildren.contains(originalEdges[2])) {
                eeA.setLength(mergedEdge1.getLength());
                originalEdges[0] = eeA;
                originalEdges[1] = itE.next();
                originalEdges[1].setLength(mergedEdge2.getLength());
            } else {
                Edge eeB = itE.next();
                eChildren.clear();
                this.getAllDescendants(eeB, null, eChildren);
                if(eChildren.contains(originalEdges[2])) {
                    eeB.setLength(mergedEdge1.getLength());
                    eeA.setLength(mergedEdge2.getLength());
                    originalEdges[0] = eeB;
                    originalEdges[1] = eeA;
                } else throw new AlgorithmException("Problem labelling edges during root move.");
            }

        } catch (AlgorithmException anErr) {
            throw new AlgorithmError("Algorithm Error: Unable to reverse root move: "+anErr.getMessage());
        }

        return originalEdges;
    }


}
