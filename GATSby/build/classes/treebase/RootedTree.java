/**
    RootedTree

    Copyright (C) 2011  Tom M. W. Nye

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

package treebase;

/**
 * RootedTree extends the Tree class by adding in a vertex labelled as the root,
 * and a map from vertices to immediate ancestors
 */

import java.util.*;

public class RootedTree extends Tree {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    private VertexWithAncestor root;

    /* ------------------------------------------*/
    // New Vertex inner class

    /* Inner class provides link to ancestor */

    public class VertexWithAncestor extends Vertex {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        protected VertexWithAncestor ancestor;

        /** Trivial constructor */
        public VertexWithAncestor(String s) {
            super(s);
            ancestor = null;
        }

        protected void setAncestor(VertexWithAncestor v) {
            ancestor = v;
        }
        public VertexWithAncestor getAncestor() {
            return ancestor;
        }

    }


    @Override
    public Vertex addNewVertex(String theLabel) {
        VertexWithAncestor v = new VertexWithAncestor(theLabel);
        vertices.add(v);
        if (theLabel!=null) {
            if (!theLabel.isEmpty()) {
                taxa.add(theLabel);
            }
        }        return v;
    }

    // End inner class
    /* ------------------------------------------*/


    protected RootedTree() {
        super();
        root = null;
    }

    public RootedTree(String theString) throws AlgorithmException {
        super();

        theString = NewickStringUtils.removeTrailingBrackets(theString);

        // Create the first vertex
        VertexWithAncestor v = (VertexWithAncestor) addNewVertex("");
        root = v;
        root.setAncestor(null);

        // Now parse the string
        generateConnectionsFromNewick(v, theString);
    }

    public RootedTree(java.io.File theFile) throws AlgorithmException {
        super();

        // Read in a string from a file
        String theString = new String();
        String s;
        try {
            java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(theFile));
            // Read a string
            while ( (s=br.readLine()) != null ) {
                theString += s;
            }
        }
        catch (java.io.IOException ioex) {
            throw new AlgorithmException("Error reading newick file.");
        }

        theString = NewickStringUtils.removeTrailingBrackets(theString);

        // Create the first vertex
        VertexWithAncestor v = (VertexWithAncestor) addNewVertex("");
        root = v;
        root.setAncestor(null);

        // Now parse the string
        generateConnectionsFromNewick(v, theString);

    }

    /** Generate vertices and edges from a newick string.
     Same as unrooted version, but ancestry added. */
    public void generateConnectionsFromNewick(VertexWithAncestor v, String theString) throws AlgorithmException {
        /* Split the string up into comma separated blocks, each with a branch length.
         For each block generate a new vertex and edges, and recursively parse each block */
         HashMap<String,Double> blocks = NewickStringUtils.parseToSubstrings(theString);

         if (blocks.size()==0) {
             // The string is already minimal -- i.e. just a taxon name
             v.label = new String(theString);
             taxa.add(v.label);
             return;
         }

         // Loop thru' blocks
         HashSet<String> vertexStrings = new HashSet<String>(blocks.keySet());
         Iterator<String> it = vertexStrings.iterator();
         for (int i=0; i<vertexStrings.size(); i++) {
             String s = it.next();
             double x = (blocks.get(s)).doubleValue();
             VertexWithAncestor w = (VertexWithAncestor) addNewVertex("");
             connect(v, w, x);
             w.setAncestor(v);
             generateConnectionsFromNewick(w, s);
         }
    }


    /** Clone.
     Pass in corresp=empty map get get info on how vertices match between original and clone,
     otherwise a null */
    public RootedTree clone(HashMap<Vertex,Vertex> corresp) {
        RootedTree t = new RootedTree();
        t.buildFromTemplate(this, corresp);
        return t;
    }
    /** Regular clone method */
    public RootedTree clone() {
        return this.clone(null);
    }

    /** Copy from a template tree */
    protected void buildFromTemplate(RootedTree t, HashMap<Vertex,Vertex> corresp) {
        if (corresp == null) {
            corresp = new HashMap<Vertex,Vertex>();
        }
        else {
            corresp.clear();
        }
        buildFromTemplate((Tree)t, corresp);
        root = (VertexWithAncestor) corresp.get(t.getRoot());
        this.recalculateAncestry();
    }
    protected void buildFromTemplate(RootedTree t) {
        buildFromTemplate(t, null);
    }
    /** Direct copy method: use with caution!!! */
    private void setVerticesEdgesAndRoot(HashSet<Vertex> v, HashSet<Edge> e, HashSet<String> t, VertexWithAncestor r) {
        vertices = v;
        edges = e;
        taxa = t;
        root = r;
    }
    @Override
    public void forceCopy(Tree t) {
        ((RootedTree)t).setVerticesEdgesAndRoot(vertices, edges, taxa, root);
    }

    /** Constructors based on rooted an unrooted tree */
    public RootedTree(Tree t, Vertex r, HashMap<Vertex,Vertex> corresp) throws AlgorithmError {
        buildFromTemplate((Tree)t, corresp);
        root = (VertexWithAncestor) corresp.get(r);
        if (root==null) throw new AlgorithmError("Root not found: unable to build rooted tree from unrooted tree");
        this.recalculateAncestry();
    }
    public RootedTree(Tree t, Vertex r) throws AlgorithmError {
        HashMap<Vertex,Vertex> corresp = new HashMap<Vertex,Vertex>();
        buildFromTemplate((Tree)t, corresp);
        root = (VertexWithAncestor) corresp.get(r);
        if (root==null) throw new AlgorithmError("Root not found: unable to build rooted tree from unrooted tree");
        this.recalculateAncestry();
     }
    public RootedTree(Tree t, Edge e, HashMap<Vertex,Vertex> corresp) throws AlgorithmError {
        Vertex r = t.divide(e, "");
        buildFromTemplate((Tree)t, corresp);
        root = (VertexWithAncestor) corresp.get(r);
        if (root==null) throw new AlgorithmError("Root not found: unable to build rooted tree from unrooted tree");
        this.recalculateAncestry();
    }
    public RootedTree(Tree t, Edge e) throws AlgorithmError {
        Vertex r = t.divide(e, "");
        HashMap<Vertex,Vertex> corresp = new HashMap<Vertex,Vertex>();
        buildFromTemplate((Tree)t, corresp);
        root = (VertexWithAncestor) corresp.get(r);
        if (root==null) throw new AlgorithmError("Root not found: unable to build rooted tree from unrooted tree");
        this.recalculateAncestry();
     }

    /** Constructor from a hashmap from splits to lengths */
    public RootedTree(HashMap<Split,Double> splitsToLengths, Split rootSplit) throws AlgorithmError {
        super();
        taxa = new HashSet<String>();

        if (splitsToLengths.size()==0) {
            VertexWithAncestor v = (VertexWithAncestor) addNewVertex("");
            return;
        }

        if (!(splitsToLengths.containsKey(rootSplit))) {
            throw new AlgorithmError("Attempt to construct rooted tree from splits failed -- root not included in split set.");
        }

        Edge rootEdge = this.buildFromSplits(splitsToLengths, rootSplit);

        // Finally add in the root
        root = (VertexWithAncestor)this.divide(rootEdge, "");
        recalculateAncestry();
        
    }

    /** Output as a string. */
    public String toOrderedString() {
        return toString(root, null, new CompareVerticesByOrdering());
    }

    /** Output as a string */
    public String toString() {
        return toString(root, null,null);
    }

    /** Output as a string, ignoring branch lengths.
     This will produce the same string every time, and so can be used to compare topologies. */
    public String toTopologyString() {
        return toString(root, null, false, new CompareVerticesByOrdering());
    }




/* -------------------------------------------------------------------------- */

    /* Methods relating to the root */

    public VertexWithAncestor getRoot() {
        return root;
    }

    /** Remove degree 2 vertex and update ancestry relations */
    public void removeDegreeTwoVertex(VertexWithAncestor v) throws AlgorithmException {
        if (v==root) return;
        Iterator<Vertex> it = v.neighboursIterator();
        VertexWithAncestor a = v.getAncestor();
        VertexWithAncestor d = null;
        VertexWithAncestor vTemp1 = (VertexWithAncestor)it.next();
        VertexWithAncestor vTemp2 = (VertexWithAncestor)it.next();
        if(vTemp1==a) d = vTemp2;
        else if(vTemp2==a) d = vTemp1;
        ((Graph)this).removeDegreeTwoVertex((Vertex)v);
        if(d!=null) d.setAncestor(a);
        else {
            /* Deal with case when the ancestry cannot be deduced. e.g. this occurs
               during SPR when the prune edge is an ancestor of the graft edge meaning
               the vertices on either side of the degree two vertex had been siblings
               on the original tree.
               This code will never be executed if the method is called on a
               connected tree - only during special circumstances such as mid-SPR.
            */
            vTemp1.setAncestor(null);
            vTemp2.setAncestor(null);
        }
    }

    /* Set the root to be vertex v and recalculate ancestry etc in correspondance */
    protected void setRoot(VertexWithAncestor v) throws AlgorithmException {
        if(!vertices.contains(v)) {
            throw new AlgorithmException("Algorithm Exception: attempting to specify root"
                    + "outside of set of vertices.");
        } else {
            root=v;
            root.setAncestor(null);
            recalculateAncestry();
        }
    }

    /** Find out whether an edge is connected to the root */
    public boolean isEdgeConnectedToRoot(Edge e) {
        Vertex[] w = e.getVertices();
        if (w[0]==root) return true;
        if (w[1]==root) return true;
        return false;
    }
    
    
   /** Find a vertex in the middle of the tree.
     Override Tree method. */
    @Override
    public Vertex pickCentralVertex(boolean topologyFlag) {
        return root;
    }
    
/* -------------------------------------------------------------------------- */

    /* Methods relating to ancestry */

    /** Recalculate ancestry map */
    public void recalculateAncestry() {
        recalculateAncestry(root, null);
    }
    private void recalculateAncestry(VertexWithAncestor v, VertexWithAncestor fromV) {
        Iterator<Vertex> it = v.neighboursIterator();
        VertexWithAncestor w;
        for (int i=0; i<v.degree(); i++) {
            w = (VertexWithAncestor) it.next();
            if (!(w.equals(fromV))) {
                w.setAncestor(v);
                recalculateAncestry(w, v);
           }
        }
    }

    /** Given an edge get the earlier vertex */
    public VertexWithAncestor getAncestralVertex(Edge e) {
        Vertex[] v = e.getVertices();
        VertexWithAncestor v1, v2;
        v1 = (VertexWithAncestor) v[0];
        v2 = (VertexWithAncestor) v[1];
        VertexWithAncestor anc =  v1.getAncestor();
        if (v2.equals(anc)) return v2;
        anc =  v2.getAncestor();
        if (v1.equals(anc)) return v1;
        javax.swing.JOptionPane.showMessageDialog(null,"Error getting ancestral vertex for an edge in a rooted tree. This shouldn't be possible. ","Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        return null;
    }

    /** Given an edge get the later vertex */
    public VertexWithAncestor getDescendantVertex(Edge e) {
        Vertex[] v = e.getVertices();
        VertexWithAncestor v1, v2;
        v1 = (VertexWithAncestor) v[0];
        v2 = (VertexWithAncestor) v[1];
        VertexWithAncestor anc = v1.getAncestor();
        if (v2.equals(anc)) return v1;
        anc = v2.getAncestor();
        if (v1.equals(anc)) return v2;
        javax.swing.JOptionPane.showMessageDialog(null,"Error getting descendant vertex for an edge in a rooted tree. This shouldn't be possible. ","Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        return null;
    }

    /** Given a vertex find ancestral edge */
    public Edge getAncestralEdge(VertexWithAncestor v) {
        VertexWithAncestor w = v.getAncestor();
        if ((w==null)&&(!(v==root))) {
            javax.swing.JOptionPane.showMessageDialog(null,"Error getting ancestor for a vertex in a rooted tree. This shouldn't be possible. ","Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
            return null;
        }
        if (v==root) return null;
        Edge e = null;
        try {
            e = v.getEdge(w);
        }
        catch (AlgorithmException anEx) {
            javax.swing.JOptionPane.showMessageDialog(null,"Error accessing edge between vertex and its ancestor. This shouldn't be possible. ","Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }
        return e;
    }

    /** Get Vertices and edges descended from a vertex */
    public void getAllDescendants(VertexWithAncestor v, HashSet<VertexWithAncestor> allV, HashSet<Edge> allE) {
        VertexWithAncestor anc = v.getAncestor();
        // Loop thru' descendants of v
        HashSet<Vertex> h = v.getNeighbours();
        Iterator<Vertex> it = h.iterator();
        VertexWithAncestor w;
        for (int i=0; i<h.size(); i++) {
            w = (VertexWithAncestor) it.next();
            if (!(w.equals(anc))) {
                /* Add w to allV and corresp edge to allE,
                   then continue recursively. */
                if (allV!=null) allV.add(w);
                if (allE!=null) {
                    try {
                        allE.add(v.getEdge(w));
                    }
                    catch (AlgorithmException anErr) {
                        javax.swing.JOptionPane.showMessageDialog(null,"Error getting descendants in a rooted tree. This shouldn't be possible. "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
                    }
                }
                getAllDescendants(w, allV, allE);
            }
        }
    }

    /** Get Vertices and edges descended from an edge */
    public void getAllDescendants(Edge e, HashSet<VertexWithAncestor> allV, HashSet<Edge> allE) {
        VertexWithAncestor v = getDescendantVertex(e);
        if (allV!=null) allV.add(v);
        getAllDescendants(v, allV, allE);
    }

    /** Get vertices and edges that are ancestral to a vertex */
    public void getAllAncestors(VertexWithAncestor v, HashSet<VertexWithAncestor> allV, HashSet<Edge> allE) {
        if (v==root) return;
        VertexWithAncestor anc = v.getAncestor();
        if (allV!=null) {
            allV.add(anc);
        }
        if (allE!=null) {
            try {
                allE.add(v.getEdge(anc));
            }
            catch (AlgorithmException anErr) {
                javax.swing.JOptionPane.showMessageDialog(null,"Error getting ancestors in a rooted tree. This shouldn't be possible. "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
            }
        }
        getAllAncestors(anc, allV, allE);
    }

    /** Get vertices and edges that are ancestral to an edge */
    public void getAllAncestors(Edge e, HashSet<VertexWithAncestor> allV, HashSet<Edge> allE) {
        VertexWithAncestor v = getAncestralVertex(e);
        if (allV!=null) allV.add(v);
        getAllAncestors(v, allV, allE);
    }

    /** Get the immediate descendants of a vertex.
     Unfortunately these are returned as Vertex not VertexAsAncestor,
     so additional casting might be necessary. */
    public HashSet<Vertex> getChildren(VertexWithAncestor v) {
        HashSet<Vertex> neighbours = v.getNeighbours();
        VertexWithAncestor a = v.getAncestor();
        if (a!=null) neighbours.remove(a);
        return neighbours;
    }


    /** Get descendant taxa of an edge */
    public HashSet<String> getDescendantTaxa(Edge e) {
        HashSet<String> dTaxa = new HashSet<String>();
        HashSet<VertexWithAncestor> allV = new HashSet<VertexWithAncestor>();
        getAllDescendants(e, allV, null);
        Iterator<VertexWithAncestor> it = allV.iterator();
        while(it.hasNext()) {
            VertexWithAncestor v = it.next();
            if(v.degree()==1) dTaxa.add(v.label);
        }
        return dTaxa;
    }

    /** Get descendant taxa of a vertex */
    public HashSet<String> getDescendantTaxa(VertexWithAncestor w) {
        HashSet<String> dTaxa = new HashSet<String>();
        HashSet<VertexWithAncestor> allV = new HashSet<VertexWithAncestor>();
        getAllDescendants(w, allV, null);
        if (w.degree()==1) allV.add(w);
        Iterator<VertexWithAncestor> it = allV.iterator();
        while(it.hasNext()) {
            VertexWithAncestor v = it.next();
            if(v.degree()==1) dTaxa.add(v.label);
        }
        return dTaxa;
    }


}
