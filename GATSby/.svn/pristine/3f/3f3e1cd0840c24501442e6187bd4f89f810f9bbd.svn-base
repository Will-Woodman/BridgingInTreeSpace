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
import treebase.RootedTree;

/**
 * Plot rooted trees
 */


public class RootedTreePlotter implements DrawTree {

    private HashMap vertexLocations; // Map from vertices to points
    private RootedTree theTree;
    private HashMap vertexOrderTags; // Map from vertex to string used for ordering
    private CompareVerticesByTag vertexComparator; // used for sorting vertices
    private double maxLen;

    private java.text.DecimalFormat df = new java.text.DecimalFormat("####.######");

    static final int NODE_RADIUS = 3;

    /** Constructor */
    public RootedTreePlotter(RootedTree t) {
        theTree = t;
        vertexOrderTags = new HashMap();
        maxLen = theTree.getMaxLeafDistance(t.getRoot(), false);
    }


    /** Draw the tree */
    public void draw(Graphics g, TreeViewer theViewer, Dimension sizeOfTree) {
        // Get storage for the vertex ordering
        vertexOrderTags.clear();
        // Get tags for vertex ordering
        getVertexOrderTags(theTree.getRoot(), null);
        vertexComparator = new CompareVerticesByTag();
        maxLen = theTree.getMaxLeafDistance(theTree.getRoot(), false);
        vertexLocations = new HashMap();

        // Recursive draw
        Point where = new Point(theViewer.BORDER_WIDTH, theViewer.BORDER_WIDTH);
        Point location = new Point();
        draw(theTree.getRoot(), g, theViewer, where, sizeOfTree, location);
    }


    /** Draw the tree -- cartesian method */
    public void draw(RootedTree.VertexWithAncestor v, Graphics g, TreeViewer theViewer, Point where, Dimension sizeOfTree, Point vertexLocation) {

        Point location = new Point();

        // Calculate the x-coord of the node
        int dx = 0;
        if (v != theTree.getRoot()) {
            Graph.Edge e = theTree.getAncestralEdge(v);
            dx = theViewer.convertBranchLenToPix(e.getLength());
        }

        int x =  where.x + dx;
        sizeOfTree.width = dx;
        int y;

        // Deal with the case of no children
        if ((v.degree() == 1)||(v.degree() == 0)) {

            // Get the coords: node is draw just below of where.y
            y = where.y + theViewer.BORDER_WIDTH;
            int label_length = this.getLongestLabelLength(g);
            sizeOfTree.width += 2*theViewer.BORDER_WIDTH + label_length;
            sizeOfTree.height = 2*theViewer.BORDER_WIDTH;
            vertexLocation.x = x;
            vertexLocation.y = y;
            location.x = x;
            location.y = y;
            vertexLocations.put(v, location);
            // Draw the line to parent
            g.drawLine(x, y, where.x, y);
            // Draw the node
            drawLeafNode(g, location, v.label);
            return;

        }
        else {
            // Order the children so that they are drawn the same way each time
            ArrayList orderedChildren = new ArrayList();
            orderedChildren.addAll(v.getNeighbours());
            if (v!=theTree.getRoot())
                orderedChildren.remove(v.getAncestor());
            Collections.sort(orderedChildren, vertexComparator);

            // Set up variables prior to loop
            Point childPosition = new Point(x,where.y);
            Dimension childDimension = new Dimension();
            Point childVertexPos = new Point();
            int ymin=0;
            int ymax=0;
            int xmax=0;
            sizeOfTree.height = 0;

             // Loop thru' the children
             ListIterator it = orderedChildren.listIterator();
             for (int i=0; i<orderedChildren.size(); i++) {
                RootedTree.VertexWithAncestor child = (RootedTree.VertexWithAncestor) it.next();
                draw(child, g, theViewer, childPosition, childDimension, childVertexPos);
                // Draw the next child below the current child
                childPosition.y += childDimension.height;
                sizeOfTree.height += childDimension.height;
                if (i==0) {
                    ymin = childVertexPos.y;
                    ymax = childVertexPos.y;
                    xmax = childDimension.width;
                }
                else {
                    if (childVertexPos.y<ymin) ymin = childVertexPos.y;
                    if (childVertexPos.y>ymax) ymax = childVertexPos.y;
                    if (childDimension.width>xmax) xmax = childDimension.width;
                }
            } // End loop thru' children
            sizeOfTree.width += xmax;
            if (v==theTree.getRoot()) sizeOfTree.height += 2*theViewer.BORDER_WIDTH;

            // get y position for the node:
            // half way between the top and bottom child nodes
            y = (ymax+ymin)/2;
            location.x = x;
            location.y = y;
            vertexLocations.put(v, location);
            vertexLocation.x = x;
            vertexLocation.y = y;

            // Draw line to parent
            g.drawLine(x, y, where.x, y);

            // Draw horizontal line
            g.drawLine(x,ymin,x,ymax);

            // Draw the node
            drawNode(g, location);

        }

    }

    /** Draw just the node -- horizontal mode*/
    protected void drawLeafNode(Graphics g, Point location, String theLabel) {
        drawNode(g, location);
        // Put a label next to the node
        g.drawString(theLabel, location.x+NODE_RADIUS+3, location.y+NODE_RADIUS+3);

    }

    /** Draw just the node */
    protected void drawNode(Graphics g, Point location) {
        // Do nothing! Override later if necessary
    }


    /** Calculate the size of the tree without drawing it */
    public java.awt.Dimension calculateSize(TreeViewer theViewer, Graphics g){
        double dt = maxLen;
        int dx = theViewer.convertBranchLenToPix(dt);
        int numDegOne = 0;

        Iterator it = theTree.getVertexIterator();
        Graph.Vertex v;
        for (int i=0; i<theTree.numVertices(); i++) {
            v = (Graph.Vertex) it.next();
            if (v.degree()==1) numDegOne++;
        }
        int dy = 2*theViewer.BORDER_WIDTH*numDegOne;
        int l = getLongestLabelLength(g);

        dx += l + 2*theViewer.BORDER_WIDTH;
        dy += 2*theViewer.BORDER_WIDTH;
        return new Dimension(dx, dy);
    }



    /** Measure the max path length */
    public double getMaxPathLength() {
        return maxLen;
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

    /** Get label length */
    public int getLabelLength(Graphics g, String s) {
        if (s==null) return 0;
        if (s.length()==0) return 0;
        java.awt.FontMetrics fm = g.getFontMetrics();
        java.awt.geom.Rectangle2D r = fm.getStringBounds(s, g);
        int w = (int) Math.round(r.getWidth());
        return w;
    }

    /** Get tooltip labels */
    public String getEdgeLabel(treebase.Graph.Edge e) {
        return df.format(e.getLength());
    }
    public String getVertexLabel(treebase.Graph.Vertex v){
        return v.label;
    }

    /** Return the Vertex at a given coordinate */
    public Graph.Vertex getVertexAtPoint(Point where) {
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
    public Graph.Edge getEdgeAtPoint(Point where) {
        RootedTree.VertexWithAncestor v;
        Point loc, upLoc;
        Iterator it = theTree.getVertexIterator();
        for (int i=0; i<theTree.numVertices(); i++) {
            v = (RootedTree.VertexWithAncestor) it.next();
            if (v != theTree.getRoot()) {
                loc = (Point) vertexLocations.get(v);
                upLoc = (Point) vertexLocations.get(v.getAncestor());

                if ((where.x<=loc.x)&&(where.x>=upLoc.x)&&(where.y>=loc.y-1)&&(where.y<=loc.y+1)) {
                    return theTree.getAncestralEdge(v);
                }
            }
        }

        return null;

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
