/**
    Graph
    Representation of a graph as a set of vertices and edges

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
 * Representation of a graph as a set of vertices and edges
 */

import java.util.*;

public class Graph implements java.io.Serializable {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    protected HashSet<Vertex> vertices;
    protected HashSet<Edge> edges;

    protected Graph() {
        vertices = new HashSet<Vertex>();
        edges = new HashSet<Edge>();
    }

    /** Clone a graph
     Corresp is a hashmap from old vertices to new.
     Pass in null if this information is not needed, otherwise an empty map */
    protected Graph clone(HashMap<Vertex,Vertex> corresp) {
        Graph g = new Graph();
        g.buildFromTemplate(this, corresp);
        return g;
    }
     protected Graph clone() {
         return this.clone(null);
     }


    /** Copy from a template graph.
     Corresp is a hashmap from old vertices to new.
     Pass in null if this information is not needed, otherwise an empty map */
    protected void buildFromTemplate(Graph g, HashMap<Vertex,Vertex> corresp) {
        vertices.clear();
        edges.clear();

        // Set up the correspondence
        if (corresp==null) {
            corresp = new HashMap<Vertex,Vertex>();
        } else {
            corresp.clear();
        }

        // Create new vertices
        Vertex v;
        Iterator<Vertex> itV = g.getVertexIterator();
        for (int i=0; i<g.numVertices(); i++) {
            v = itV.next();
            Vertex vNew = addNewVertex(v.label);
            corresp.put(v, vNew);
         }

        // Deal with edges
        Edge e;
        try {
            Iterator<Edge> itE = g.getEdgeIterator();
            for (int i=0; i<g.numEdges(); i++) {
                e = itE.next();
                e.clone(corresp, this);
            }
        }
        catch (AlgorithmError anErr) {
            javax.swing.JOptionPane.showMessageDialog(null,"Error cloning an edge. This shouldn't be possible. "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }
    }
    protected void buildFromTemplate(Graph g) {
        buildFromTemplate(g, null);
    }


    /** Util methods */
    public int numVertices() {
        return vertices.size();
    }
    public int numEdges() {
        return edges.size();
    }
    public HashSet<Vertex> getVertices() {
        HashSet h = new HashSet();
        h.addAll(vertices);
        return h;
    }
    public HashSet<Edge> getEdges() {
        HashSet h = new HashSet();
        h.addAll(edges);
        return h;
    }
    public Iterator<Vertex> getVertexIterator() {
        return vertices.iterator();
    }
    public Iterator<Edge> getEdgeIterator() {
        return edges.iterator();
    }

    /** Use the following for getting the complement of some edge or vertex set */
    public HashSet<Edge> getEdgeSetComplement(HashSet<Edge> s) {
        HashSet<Edge> c = getEdges();
        c.removeAll(s);
        return c;
    }
    public HashSet<Vertex> getVertexSetComplement(HashSet<Vertex> s) {
        HashSet<Vertex> c = getVertices();
        c.removeAll(s);
        return c;
    }

    public Vertex findVertexByLabel(String theLabel) {
        Iterator<Vertex> it = vertices.iterator();
        Vertex v = null;
        while (it.hasNext()) {
            v = it.next();
            if (v.label.equals(theLabel)) return v;
        }
        return null;
    }

    public double getTotalEdgeLength() {
        double s = 0;
        Edge e;
        Iterator<Edge> it = edges.iterator();
        for (int i=0; i<edges.size(); i++) {
            e = it.next();
            s += e.getLength();
        }
        return s;
    }

    /** Scale all edges in the tree by a certain factor */
    public void scale(double fac) {
        Edge e;
        double s;
        Iterator<Edge> it = edges.iterator();
        for (int i=0; i<edges.size(); i++) {
            e = it.next();
            s = e.getLength()*fac;
            e.setLength(s);
        }
    }

    /* -------------------------------------------------------------------------- */

    /* Topological operations! */

    /** Add a vertex
     This allows you to construct graphs vertex-by-vertex */
    public Vertex addNewVertex(String theLabel)  {
        Vertex v = new Vertex(theLabel);
        vertices.add(v);
        return v;
    }


    /** Connect two vertices with an edge!
     Override this for extending classes. */
    public Edge connect(Vertex v, Vertex w, double x) throws AlgorithmError {
        if ((!vertices.contains(v))||(!vertices.contains(w))) {
            throw new AlgorithmError("Cannot join vertices from outside a graph.");
        }
        Edge e = new Edge(v,w,x);
        addEdge(e);
        return e;
    }
    /** Add an edge */
    protected void addEdge(Edge e) {
        e.vertex1.addConnection(e.vertex2,e);
        e.vertex2.addConnection(e.vertex1,e);
        edges.add(e);
    }

    /** Delete an edge */
    public void cut(Edge e) throws AlgorithmError {
        if (!edges.contains(e)) {
            throw new AlgorithmError("Cannot cut an edge not included in graph.");
        }
        e.disconnect();
    }

    /** Remove a vertex */
    public void remove(Vertex v) throws AlgorithmError {
        if (!vertices.contains(v)) {
            throw new AlgorithmError("Cannot remove a vertex not included in graph.");
        }
        HashSet<Edge> nbh = v.getEdges();
        Iterator<Edge> it = nbh.iterator();
        for (int i=0; i<nbh.size(); i++) {
            Edge e = it.next();
            e.disconnect();
        }
        vertices.remove(v);
    }

    /** Absorb a single degree 2 vertex */
    public void removeDegreeTwoVertex(Vertex v) throws AlgorithmException {
        if (!vertices.contains(v)) {
            throw new AlgorithmException("Attempt to absorb a vertex that isn't in the enveloping tree.");
        }
        if (v.degree()!=2) {
            throw new AlgorithmException("Attempt to absorb a vertex that doesn't have degree 2.");
        }

        Vertex u,w;
        Edge eu, ew;
        HashSet<Vertex> nb;
        Iterator<Vertex> itN;
        double x;

        // Get the two neighbours of v
        nb = v.getNeighbours();
        itN = nb.iterator();
        u = itN.next();
        w = itN.next();
        // Get corresponding edges
        eu = v.getEdge(u);
        ew = v.getEdge(w);
        x = eu.getLength()+ew.getLength();
        cut(eu);
        cut(ew);
        vertices.remove(v);
        connect(u, w, x);
    }

    /** Remove all vertices with degree 2 */
    public void removeDegreeTwoVertices() {
        // Find deg 2 vertices
        HashSet<Vertex> toRemove = new HashSet<Vertex>();
        Vertex v;
        Iterator<Vertex> it = vertices.iterator();
        for (int i=0; i<vertices.size(); i++) {
            v = it.next();
            if (v.degree()==2) toRemove.add(v);
        }

        // Loop thru' deg 2 vertices
        it = toRemove.iterator();
        try {
            for (int i=0; i<toRemove.size(); i++) {
                v = it.next();
                this.removeDegreeTwoVertex(v);
            }
        }
        catch (AlgorithmException anErr) {
            javax.swing.JOptionPane.showMessageDialog(null,"Serious problem cropping degree two vertices: "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
            // System.out.println("Serious problem cropping degree two vertices: "+anErr.getMessage());
        }
    }

    /** Contract an edge down */
    public Vertex contractEdge(Edge e) throws AlgorithmError {
        if (!edges.contains(e)) {
            throw new AlgorithmError("Cannot contract an edge not included in graph.");
        }
        return e.contract();
    }

    /** Merge two vertices */
    private Vertex mergeVertices(Vertex v, Vertex w) {
        Vertex res = addNewVertex("["+v.label+","+w.label+"]"); // replacement vertex

        Vertex u;
        Edge e;
        double x = 0.0;

        try {

            HashSet<Vertex> nbh = v.getNeighbours();
            Iterator<Vertex> it = nbh.iterator();
            for (int i=0; i<nbh.size(); i++) {
                u = it.next();
                e = v.getEdge(u);
                x = e.getLength();
                e.disconnect();
                if ((!(u.equals(v)))&&(!(u.equals(w)))) {
                    this.connect(res, u, x);
                }
                if (u.equals(v)) {
                    this.connect(res, res, x);
                }
            }

            nbh = w.getNeighbours();
            it = nbh.iterator();
            for (int i=0; i<nbh.size(); i++) {
                u = it.next();
                e = w.getEdge(u);
                x = e.getLength();
                e.disconnect();
                if ((!(u.equals(v)))&&(!(u.equals(w)))) {
                    this.connect(res, u, x);
                }
                if (u.equals(w)) {
                    this.connect(res, res, x);
                }
            }

            vertices.remove(v);
            vertices.remove(w);

        }
        catch (AlgorithmException anErr) {
             javax.swing.JOptionPane.showMessageDialog(null,"Serious problem merging two vertices: "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }

        return res;
    }

    /** Add a degree 2 vertex in the middle of an edge */
    public Vertex divide(Edge e, String s) {
         return e.divide(s);
    }

    /* -------------------------------------------------------------------------- */

    /* Start vertex inner class */

    public class Vertex implements java.io.Serializable {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        protected HashMap<Vertex,Edge> connections; // Map from vertex neighbours to corresponding edges NB: implies at most one edge between every pair of vertices!
        public String label; // Taxon label at leaf

        /** Return a hash set containing a copy of the neighbours */
        public HashSet<Vertex> getNeighbours() {
            HashSet res = new HashSet();
            res.addAll(connections.keySet());
            return res;
        }

        /** Given an adjacent edge get the corresponding vertex */
        public Vertex getNeighbour(Edge e) {
            if (!connections.containsValue(e)) return null;
            Iterator<Vertex> it = connections.keySet().iterator();
            for (int i=0; i<connections.size(); i++) {
                Vertex v = it.next();
                Edge f = connections.get(v);
                if (e==f) return v;
            }
            return null;
        }

        /** Return an iterator through the neighbours */
        public Iterator<Vertex> neighboursIterator() {
            return connections.keySet().iterator();
        }

        /** Return edge corresponding to a neighbour */
        public Edge getEdge(Vertex v) throws AlgorithmException {
            if (!connections.containsKey(v)) {
                throw new AlgorithmException("Request for edge object for non-neighbouring vertices");
            }
            return connections.get(v);
        }

        /** Return a hashset of all attached edges */
        public HashSet<Edge> getEdges() {
            HashSet<Edge> h = new HashSet<Edge>();
            h.addAll(connections.values());
            return h;
        }

        /** Return an iterator thru attached edges */
        public Iterator<Edge> edgeIterator() {
            return connections.values().iterator();
        }

        /** Test for connection to another vertex */
        public boolean isConnectedTo(Vertex v) {
            return connections.keySet().contains(v);
        }

        /** Connect to another vertex via edge.
            Avoid using this! */
        protected void addConnection(Vertex w, Edge e) {
            connections.put(w, e);
        }

        /** Remove a neighbour: called by Edge.disconnect called by Graph.cut */
        private void removeNeighbour(Vertex v) {
            edges.remove(connections.get(v));
            connections.remove(v);
        }

        /** Get degree */
        public int degree() {
            return connections.size();
        }

        /** Trivial constructor */
        protected Vertex(String s) {
            connections = new HashMap();
            label = s;
            vertices.add(this);
        }

        /** Return the distance of a neighbour or -1 if not connected */
        public double getNeighbourDistance(Vertex v) {
            Object e = connections.get(v);
            if (e==null) return -1.0;
            return ((Edge)e).getLength();
        }

        /** toString currently used for debugging */
        public String toString() {
            String s = new String();
            if (label.isEmpty()) {
                s += "Unlabelled vertex ";
            }
            else {
                s = s + label + " ";
            }
            s += "degree "+connections.size();
            return s;
        }

    }


    /* -------------------------------------------------------------------------- */

    /* Start edge inner class */

    public class Edge implements java.io.Serializable{
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;

        protected Vertex vertex1;
        protected Vertex vertex2;
        private double length;

       /** Main constructor */
        protected Edge(Vertex v, Vertex w, double x) {
            vertex1 = v;
            vertex2 = w;
            length = x;
        }

        /** toString currently used for debugging */
        public String toString() {
            String s = new String();
            s += "'"+vertex1.label+"' connected to '"+vertex2.label+"' length="+length;
            return s;
        }

        /** Provide access to (unordered) vertices */
        public Vertex[] getVertices() {
            Vertex[] h = new Vertex[2];
            h[0] = vertex1;
            h[1] = vertex2;
            return h;
        }

        /** Return length */
        public double getLength() {
            return length;
        }

        /** Set length */
        public void setLength(double x) {
            length = x;
        }

        /** Remove edge from corresponding vertices */
        private void disconnect() {
            edges.remove(this);
            vertex1.removeNeighbour(vertex2);
            vertex2.removeNeighbour(vertex1);
        }

        /** Contract down an edge an replace the two vertices with just one */
        private Vertex contract() {
            edges.remove(this);
            vertex1.removeNeighbour(vertex2);
            vertex2.removeNeighbour(vertex1);
            return mergeVertices(vertex1, vertex2);
        }

        /** Clone an edge
        Corresp is a hashmap from old vertices to new.  */
        private void clone(HashMap<Vertex,Vertex> correspondence, Graph g) throws AlgorithmError {
            Vertex w1 = correspondence.get(vertex1);
            Vertex w2 = correspondence.get(vertex2);
            g.connect(w1, w2, this.length);
        }

        /** Add a degree 2 vertex to an edge */
        private Vertex divide(String s) {
            this.disconnect();
            Vertex v = addNewVertex(s);
            double x = this.length*0.5;
            try {
                connect(vertex1, v, x);
                connect(v, vertex2, x);
            }
            catch (AlgorithmError anErr) {
                javax.swing.JOptionPane.showMessageDialog(null,"Error subdividing an edge. This should not be possible. "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
            }
            return v;
        }

        /** Get neighbouring edges */
        public HashSet<Edge> getNeighbouringEdges() {
            HashSet<Edge> h = new HashSet();
            h.addAll(vertex1.getEdges());
            h.addAll(vertex2.getEdges());
            h.remove(this);
            return h;
        }

        /** Change one end of an edge */
        public void changeVertex(Vertex v, Vertex w) throws AlgorithmError {
            if (vertex1==v) {
                this.disconnect();
                vertex1 = w;
                addEdge(this);
                return;
            }
            if (vertex2==v) {
                this.disconnect();
                vertex2 = w;
                addEdge(this);
                return;
            }
            throw new AlgorithmError("Bad change to one end of an edge");
        }

    }

    /* -------------------------------------------------------------------------- */

}
