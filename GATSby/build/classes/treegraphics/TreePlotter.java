/**
    TreePlotter

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

package treegraphics;

import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Set;
import treebase.Graph;
import treebase.Graph.Vertex;
import treebase.Tree;

/**
 * Plot unrooted trees
 */



public class TreePlotter implements DrawTree {

    protected HashMap<Vertex, Point> vertexLocations; // Map from vertices to points
    protected Tree theTree;
    private HashMap vertexOrderTags; // Map from vertex to string used for ordering
    private CompareVerticesByTag vertexComparator; // used for sorting vertices
    private Graph.Vertex centralVertex;
    private double maxLen;

    private java.text.DecimalFormat df = new java.text.DecimalFormat("####.######");

    static final int NODE_RADIUS = 3;


    /** Constructor */
    public TreePlotter(Tree t) {
        theTree = t;
        vertexOrderTags = new HashMap();
        centralVertex = theTree.pickCentralVertex(false);
        maxLen = theTree.getMaxLeafDistance(centralVertex, false);
    }

    /** Draw the tree */
    public void draw(Graphics g, TreeViewer theViewer, Dimension sizeOfTree) {
        Graph.Vertex centralVertex = theTree.pickCentralVertex(false);
        maxLen = theTree.getMaxLeafDistance(centralVertex, false);
        vertexLocations = new HashMap();

        // Get storage for the vertex ordering
        vertexOrderTags.clear();
        // Get tags for vertex ordering
        getVertexOrderTags(centralVertex, null);
        vertexComparator = new CompareVerticesByTag();
 
        Dimension s = calculateSize(theViewer, g);
        sizeOfTree.width = s.width;
        sizeOfTree.height = s.height;
        Point where = new Point(s.width/2, s.height/2);
        double direction = 0.0;
        int n = theTree.numTaxa();
        draw(centralVertex, null, g, theViewer, where, n, direction);
    }

    /** Recursive draw */
    private void draw(Graph.Vertex v, Graph.Vertex fromV, Graphics g, TreeViewer theViewer, Point where, int totalNumLeaves, double direction) {

        Point location = new Point(where.x, where.y);
        vertexLocations.put(v, location);

        // Deal with the case of no children
        if ((v.degree() == 1)||(v.degree() == 0)) {
            // Draw the node
            drawLeafNode(g, v, theViewer, location, direction);
            return;
        }
        else {
            // Order the children so that they are drawn the same way each time
            ArrayList orderedChildren = new ArrayList();
            orderedChildren.addAll(v.getNeighbours());
            if (fromV != null) orderedChildren.remove(fromV);
            Collections.sort(orderedChildren, vertexComparator);

            // Set up variables prior to loop
            double angularExtent = (6.2832*countLeaves(v,fromV))/(1.0*totalNumLeaves);
            double currentAngle = direction - 0.5*angularExtent;

            Graph.Vertex child;

             // Loop thru' the children
             ListIterator it = orderedChildren.listIterator();
             for (int i=0; i<orderedChildren.size(); i++) {
                child = (Graph.Vertex) it.next();
                double childAngularExtent = (6.2832*countLeaves(child,v))/(1.0*totalNumLeaves);
                double childDirection = currentAngle + 0.5*childAngularExtent;
                currentAngle += childAngularExtent;
                Point childPosition = new Point();
                double dt = v.getNeighbourDistance(child);
                childPosition.x = location.x + theViewer.polarToX(childDirection, dt);
                childPosition.y = location.y + theViewer.polarToY(childDirection, dt);
                this.draw(child, v, g, theViewer, childPosition, totalNumLeaves, childDirection);
                // Draw line to parent
                g.drawLine(location.x, location.y, childPosition.x, childPosition.y);
                //child.drawNode(g);
            } // End loop thru' children
            // Draw the node
            drawNode(g, v, location);
        }

    }

    /** Draw just the node */
    protected void drawNode(Graphics g, Graph.Vertex v, Point location) {
        // Do nothing! Override later if necessary
    }


    /** Draw node with a taxon label */
    private void drawLeafNode(Graphics g, Graph.Vertex v, TreeViewer theViewer, Point location, double direction) {
        drawNode(g, v, location);
        Graphics2D g2d = (Graphics2D)g;
        int p=0;
        int q=0;
        if (Math.cos(direction)>0.0) {
            p = location.x + (int) Math.round(theViewer.BORDER_WIDTH*Math.cos(direction));
            q = location.y + (int) Math.round(theViewer.BORDER_WIDTH*Math.sin(direction));
        }
        else {
            p = location.x + (int) Math.round((this.getLabelLength(g, v.label)+theViewer.BORDER_WIDTH)*Math.cos(direction));
            q = location.y + (int) Math.round((this.getLabelLength(g, v.label)+theViewer.BORDER_WIDTH)*Math.sin(direction));
            direction += 3.14159;
        }
        g2d.translate(p, q);
        g2d.rotate(direction);
        g2d.translate(-p, -q);
        g2d.drawString(v.label, p, q);
        g2d.translate(p, q);
        g2d.rotate(-direction);
        g2d.translate(-p, -q);

    }

    /** Get the length of a string */
    private static int getLabelLength(Graphics g, String s) {
        if (s==null) return 0;
        if (s.length()==0) return 0;
        java.awt.FontMetrics fm = g.getFontMetrics();
        java.awt.geom.Rectangle2D r = fm.getStringBounds(s, g);
        int w = (int) Math.round(r.getWidth());
        return w;
    }

    /** Get the longest label length */
    public int getLongestLabelLength(Graphics g) {
        if (g==null) return 0;
        FontMetrics fm = g.getFontMetrics();
        HashSet theLeaves = theTree.getTaxa();
        int longest = 0;
        Iterator it = theLeaves.iterator();
        for (int i=0; i<theLeaves.size(); i++) {
            String name = (String) it.next();
            java.awt.geom.Rectangle2D r = fm.getStringBounds(name, g);
            int w = (int) Math.round(r.getWidth());
            longest = (w>longest ? w : longest);
        }
        return longest+20 ;
    }

    /** Get the leaves descended from a vertex */
    private void getLeaves(HashSet theLeaves, Graph.Vertex v, Graph.Vertex w) {
        if (v.degree()<=1) {
            theLeaves.add(v);
            return;
        }

        HashSet h = v.getNeighbours();
        if (w != null) h.remove(w);
        Iterator it = h.iterator();
        Graph.Vertex u;
        for (int i=0; i<h.size(); i++) {
            u = (Graph.Vertex) it.next();
            getLeaves(theLeaves, u, v);
        }
    }
    private int countLeaves(Graph.Vertex v, Graph.Vertex w) {
        HashSet h = new HashSet();
        getLeaves(h,v,w);
        return h.size();
    }

    /** Calculate the size of the tree without drawing it */
    public java.awt.Dimension calculateSize(TreeViewer theViewer, Graphics g){
        double dt = maxLen;
        int dx = theViewer.polarToX(0.0, dt);
        int dy = theViewer.polarToY(1.571, dt);
        int l = getLongestLabelLength(g);

        dx += l + theViewer.BORDER_WIDTH;
        dy += l + theViewer.BORDER_WIDTH;
        return new Dimension(2*dx, 2*dy);
    }

    /** Measure the max path length */
    public double getMaxPathLength() {
        return maxLen;
    }


    /** Return the Vertex at a given coordinate */
    public Graph.Vertex getVertexAtPoint(java.awt.Point where){
        Set vertices = vertexLocations.keySet();
        Iterator it = vertices.iterator();
        for (int i=0; i<vertices.size(); i++) {
            Graph.Vertex v = (Graph.Vertex) it.next();
            Point location = (Point) vertexLocations.get(v);
            if ((Math.abs(where.x-location.x)<=NODE_RADIUS)&&(Math.abs(where.y-location.y)<=NODE_RADIUS)) return v;
        }
        return null;
    }

    /** Return the Edge at a given coordinate */
    public Graph.Edge getEdgeAtPoint(java.awt.Point where){
        Graph.Edge e;
        Graph.Vertex[] hv;
        Point loc1, loc2;
        Iterator ith;
        Iterator it = theTree.getEdgeIterator();
        for (int i=0; i<theTree.numEdges(); i++) {
            e = (Graph.Edge) it.next();
            hv = e.getVertices();
            loc1 = (Point) vertexLocations.get(hv[0]);
            loc2 = (Point) vertexLocations.get(hv[1]);

            int dx = loc2.x - loc1.x;
            int dy = loc2.y - loc1.y;
            double l = Math.sqrt(dx*dx + dy*dy);
            double x1,x2,y1,y2;
            x1 = dx/l;
            y1 = dy/l;
            x2 = -y1;
            y2 = x1;
            dx = where.x - loc1.x;
            dy = where.y - loc1.y;
            double comp1 = dx*x1 + dy*y1;
            double comp2 = dx*x2 + dy*y2;

            double d1, d2;
            d1 = l;
            d2 = 1.5;

            if ((comp1<=d1)&&(comp1>=0.0)&&(comp2<=d2)&&(comp2>=-d2)) return e;
        }
        return null;
    }

    /** Get tooltip labels */
    public String getEdgeLabel(treebase.Graph.Edge e) {
        return df.format(e.getLength());
    }
    public String getVertexLabel(treebase.Graph.Vertex v){
        return v.label;
    }

    /* Order the set of vertices: recursive call down tree  */
    private String getVertexOrderTags(Graph.Vertex v, Graph.Vertex fromV) {
        if (v.degree()<2) {
            vertexOrderTags.put(v, v.label);
            return v.label;
        }

        HashSet h = v.getNeighbours();
        if (fromV != null) h.remove(fromV);
        Iterator it = h.iterator();
        Graph.Vertex w;
        String maxTag=""; 
        String childTag = "";
        for (int i=0; i<h.size(); i++) {
            w = (Graph.Vertex) it.next();
            childTag = getVertexOrderTags(w, v);
            if (i==0) maxTag = childTag;
            if (childTag.compareTo(maxTag)>0) maxTag = childTag;
        }
        vertexOrderTags.put(v, maxTag);
        return maxTag;
    }
    
    private class CompareVerticesByTag implements Comparator {
        public CompareVerticesByTag() {
        }
        public int compare(Object o1, Object o2) {
            String s1 = (String) vertexOrderTags.get(o1);
            String s2 = (String) vertexOrderTags.get(o2);
            if ((s1==null)&&(s2==null)) return 0;
            if (s1==null) return -1;
            if (s2==null) return 1;
            return s1.compareTo(s2);
        }
    }


}
