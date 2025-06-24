/* PathPlotter.java

    Copyright (C) 2013  Tom M. W. Nye

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

package geodesicgraphics;

/**
 * Implementation of drawTree interface for plotting a collection of trees along a path.
 * Takes care of rotational ordering and scaling.
 * Adapted from previous version.
 */

import treebase.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;


import treebase.TreeAsSplits;

public class PathPlotter implements treegraphics.DrawTree {

    private RootedTree[] theTrees; // Set of all trees
    private int numFrames; // Number of trees along the path
    private TreeSpacePath thePath; // The path object
    private double[] paramVals; // Param values along the path
    private double maxLen; // Max tree length
    private RootedTree currentTree;
    private int indexCurrentTree;

    private HashMap<Graph.Vertex, java.awt.Point> vertexLocations; // Map from vertices to points

    static final int NODE_RADIUS = 3;

    /* Stuff related to ordering */
    private ArrayList<String> rotationalOrder;
    private String minTaxon;
    private HashMap<Graph.Vertex, String> vertexOrderTags; // Map from vertex to string used for ordering vertices
    private CompareVerticesByTag vertexComparator; // used for sorting vertices
    private CompareTaxonSetsByMinTaxon taxonSetComparator;


    /** Constructor from a TreeSpacePath object */
    public PathPlotter(TreeSpacePath p, double[] x) {

        thePath = p;
        paramVals = x;
        numFrames = x.length;

        /* Initialize the set of trees */
        theTrees = new RootedTree[numFrames];
        Split rootSplit = thePath.getRootSplit();
        double l;
        maxLen = 0.0;
        try {
            for (int i=0; i<numFrames; i++) {
                TreeAsSplits theMap = thePath.getTreeOnPath(paramVals[i]);
                RootedTree theTree = new RootedTree(theMap.getMap(), rootSplit);
                theTrees[i] = theTree;
                l = theTree.getMaxLeafDistance(theTree.getRoot(), false); // Get the maximum distance of a leaf from the root
                if (l>maxLen) maxLen = l;
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error making a PathPlotter - bad root. "+anErr.getMessage());
        }

        indexCurrentTree = 0;
        currentTree = theTrees[indexCurrentTree];

        /* Related to ordering */
        HashSet<String> taxa = theTrees[0].getTaxa();
        ArrayList<String> orderedTaxa = new ArrayList();
        orderedTaxa.addAll(taxa);
        Collections.sort(orderedTaxa);
        minTaxon = orderedTaxa.get(0);
        taxonSetComparator = new CompareTaxonSetsByMinTaxon();
        rotationalOrder = thePath.getRotationalOrder();
        vertexOrderTags = new HashMap();
        try {
            for (int i=0; i<numFrames; i++) {
                setUpVertexOrdering(theTrees[i]);
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error obtaining vertex ordering for geodesic graphics. "+anErr.getMessage());
        }
        vertexComparator = new CompareVerticesByTag();

    }

    /** Compare pairs of vertices using the tags set by setUpVertexOrdering */
    private class CompareVerticesByTag implements Comparator {
        public CompareVerticesByTag() {
        }
        public int compare(Object o1, Object o2) {
            String s1 = vertexOrderTags.get((Graph.Vertex)o1);
            String s2 = vertexOrderTags.get((Graph.Vertex)o2);
            if ((s1==null)||(s2==null)) System.out.println("AlgorithmError with  setUpVertexOrdering: "+
                    "apparent request to compare non-adjacent vertices.");
            return s1.compareTo(s2);
        }
    }


    /** Calculate the size of the tree without drawing it */
    public java.awt.Dimension calculateSize(treegraphics.TreeViewer theViewer, java.awt.Graphics g){
        double dt = maxLen;
        int dx = theViewer.polarToX(0.0, dt);
        int dy = theViewer.polarToY(1.571, dt);
        int l = getLongestLabelLength(g);

        dx += l + theViewer.BORDER_WIDTH;
        dy += l + theViewer.BORDER_WIDTH;
        return new java.awt.Dimension(2*dx, 2*dy);
    }

    /** Measure the max path length */
    public double getMaxPathLength() {
        return maxLen;
    }

     /** Get the length of a string */
    private static int getLabelLength(java.awt.Graphics g, String s) {
        if (s==null) return 0;
        if (s.length()==0) return 0;
        java.awt.FontMetrics fm = g.getFontMetrics();
        java.awt.geom.Rectangle2D r = fm.getStringBounds(s, g);
        int w = (int) Math.round(r.getWidth());
        return w;
    }

    /** Get the longest label length */
    private int getLongestLabelLength(java.awt.Graphics g) {
        if (g==null) return 0;
        java.awt.FontMetrics fm = g.getFontMetrics();
        HashSet theLeaves = theTrees[0].getTaxa();
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

    /** Get tooltip labels */
    public String getEdgeLabel(treebase.Graph.Edge e) {
        return String.format("%7.7f", e.getLength());
    }
    public String getVertexLabel(treebase.Graph.Vertex v){
        return v.label;
    }

    /** Return the Vertex at a given coordinate */
    public Graph.Vertex getVertexAtPoint(java.awt.Point where){
        Set<Graph.Vertex> vertices = vertexLocations.keySet();
        Iterator<Graph.Vertex> it = vertices.iterator();
        for (int i=0; i<vertices.size(); i++) {
            Graph.Vertex v = it.next();
            java.awt.Point location = vertexLocations.get(v);
            if ((Math.abs(where.x-location.x)<=NODE_RADIUS)&&(Math.abs(where.y-location.y)<=NODE_RADIUS)) return v;
        }
        return null;
    }

    /** Return the Edge at a given coordinate */
    public Graph.Edge getEdgeAtPoint(java.awt.Point where){
        HashSet<Graph.Edge> edges = currentTree.getEdges();
        Graph.Edge e;
        Graph.Vertex[] hv;
        Graph.Vertex v1, v2;
        java.awt.Point loc1, loc2;
        Iterator<Graph.Edge> it = edges.iterator();
        for (int i=0; i<edges.size(); i++) {
            e = it.next();
            hv = e.getVertices();
            loc1 = vertexLocations.get(hv[0]);
            loc2 = vertexLocations.get(hv[1]);
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
            d2 = 1.0;

            if ((comp1<=d1)&&(comp1>=0.0)&&(comp2<=d2)&&(comp2>=-d2)) return e;
        }
        return null;
    }

    /** Change to a particular tree */
    public void changeTree(int k) {
        currentTree = theTrees[k];
        indexCurrentTree = k;
    }


    /** Draw the tree */
    public void draw(java.awt.Graphics g, treegraphics.TreeViewer theViewer, java.awt.Dimension sizeOfTree) {

        /* Clear the store of where the vertices are */
        vertexLocations = new HashMap();

        /* Dimension calcs */
        java.awt.Dimension s = calculateSize(theViewer, g);
        sizeOfTree.width = s.width;
        sizeOfTree.height = s.height;
        java.awt.Point where = new java.awt.Point(s.width/2, s.height/2);
        double direction = 0.0;

        /* Recursive draw */
        draw(currentTree.getRoot(), g, theViewer, where, direction);
    }

    /** Recursive draw */
    private void draw(RootedTree.VertexWithAncestor v, java.awt.Graphics g, treegraphics.TreeViewer theViewer, java.awt.Point where, double direction) {

        java.awt.Point location = new java.awt.Point(where.x, where.y);
        vertexLocations.put(v, location);

        // Deal with the case of no children
        if ((v.degree() == 1)||(v.degree() == 0)) {
            // Draw the node
            drawLeafNode(g, theViewer, location, direction, v.label);
            return;
        }
        else {
            // Order the children so that they are drawn the same way each time
            ArrayList<Graph.Vertex> orderedChildren = new ArrayList();
            orderedChildren.addAll(currentTree.getChildren(v));
            Collections.sort(orderedChildren, vertexComparator);

            // Set up variables prior to loop
            double angularExtent = (6.2832*countLeaves(v, currentTree))/(1.0*currentTree.numTaxa());
            double currentAngle = direction - 0.5*angularExtent;

            RootedTree.VertexWithAncestor child;

             // Loop thru' the children
             for (int i=0; i<orderedChildren.size(); i++) {
                child = (RootedTree.VertexWithAncestor) orderedChildren.get(i);
                double childAngularExtent = (6.2832*countLeaves(child, currentTree))/(1.0*currentTree.numTaxa());
                double childDirection = currentAngle + 0.5*childAngularExtent;
                if (v==currentTree.getRoot()) {
                    childDirection = i*3.1416;
                }
                currentAngle += childAngularExtent;
                java.awt.Point childPosition = new java.awt.Point();
                double dt = v.getNeighbourDistance(child);
                childPosition.x = location.x + theViewer.polarToX(childDirection, dt);
                childPosition.y = location.y + theViewer.polarToY(childDirection, dt);
                this.draw(child, g, theViewer, childPosition, childDirection);
                // Draw line to parent
                g.drawLine(location.x, location.y, childPosition.x, childPosition.y);
                //child.drawNode(g);
            } // End loop thru' children
            // Draw the node
            drawNode(g, location);
        }

    }

    /** Draw just the node */
    protected void drawNode(java.awt.Graphics g, java.awt.Point location) {
        // Do nothing! Override later if necessary
    }


    /** Draw node with a taxon label */
    private void drawLeafNode(java.awt.Graphics g, treegraphics.TreeViewer theViewer, java.awt.Point location, double direction, String theLabel) {
        drawNode(g, location);
        double theta = direction;
        java.awt.Graphics2D g2d = (java.awt.Graphics2D)g;
        int p=0;
        int q=0;
        if (Math.cos(direction)>0.0) {
            p = location.x + (int) Math.round(theViewer.BORDER_WIDTH*Math.cos(direction));
            q = location.y + (int) Math.round(theViewer.BORDER_WIDTH*Math.sin(direction));
        }
        else {
            p = location.x + (int) Math.round((this.getLabelLength(g, theLabel)+theViewer.BORDER_WIDTH)*Math.cos(direction));
            q = location.y + (int) Math.round((this.getLabelLength(g, theLabel)+theViewer.BORDER_WIDTH)*Math.sin(direction));
            direction += 3.14159;
        }
        g2d.translate(p, q);
        g2d.rotate(direction);
        g2d.translate(-p, -q);
        g2d.drawString(theLabel, p, q);
        g2d.translate(p, q);
        g2d.rotate(-direction);
        g2d.translate(-p, -q);
    }

    private static int countLeaves(RootedTree.VertexWithAncestor v, RootedTree theTree) {
        HashSet<String> h = theTree.getDescendantTaxa(v);
        return h.size();
    }


    /** Given a tree label its vertices in the best way in order to match the circular ordering. */
    private void setUpVertexOrdering(RootedTree theTree) throws AlgorithmError {
        RootedTree.VertexWithAncestor root = theTree.getRoot();
        if (root.degree()!=2) throw new AlgorithmError("Root did not have degree 2 in PathPlotter.");

        // Order the initial two vertices
        Iterator it = root.neighboursIterator();
        RootedTree.VertexWithAncestor vA = (RootedTree.VertexWithAncestor) it.next();
        RootedTree.VertexWithAncestor vB = (RootedTree.VertexWithAncestor) it.next();
        HashSet<String> vAtaxa = theTree.getDescendantTaxa(vA);
        if (vAtaxa.contains(minTaxon)) {
            // Order vA then vB
            vertexOrderTags.put(vA, "A");
            vertexOrderTags.put(vB, "B");
        }
        else {
            // Order vB then vA
            vertexOrderTags.put(vA, "B");
            vertexOrderTags.put(vB, "A");
        }

        /* Recursively order remaining vertices */
        ArrayList<Graph.Vertex> orderedVertices = theTree.getOrderedListOfVertices();
        recursiveSetUpVertexOrdering(vA, theTree, orderedVertices);
        recursiveSetUpVertexOrdering(vB, theTree, orderedVertices);

    }

    /** Set up labels to order vertices */
    private void recursiveSetUpVertexOrdering(RootedTree.VertexWithAncestor v, RootedTree theTree, ArrayList<Graph.Vertex> orderedVertices) {
        if (v.degree()<3) {
            // No ordering necessary!
            return;
        }

        if (v.degree()==3) {
            /* Deg 3 case: just two options */
            HashSet<Graph.Vertex> children = theTree.getChildren(v);
            RootedTree.VertexWithAncestor vA, vB;
            Iterator it = children.iterator();
            vA = (RootedTree.VertexWithAncestor) it.next();
            vB = (RootedTree.VertexWithAncestor) it.next();

            HashSet<String> tX = new HashSet();
            tX.addAll(rotationalOrder);
            HashSet<String> tA = theTree.getDescendantTaxa(vA);
            HashSet<String> tB = theTree.getDescendantTaxa(vB);
            tX.removeAll(tA);
            tX.removeAll(tB);

            int s1 = scoreVertexPairAgainstRotationalOrder(tA, tB, tX);
            int s2 = scoreVertexPairAgainstRotationalOrder(tB, tA, tX);

            if (s1<s2) {
                vertexOrderTags.put(vA, "A");
                vertexOrderTags.put(vB, "B");
            }
            else if (s2<s1) {
                vertexOrderTags.put(vA, "B");
                vertexOrderTags.put(vB, "A");
            }
            else {
                // Scores equal.
                if (orderedVertices.indexOf(vA)<orderedVertices.indexOf(vB)) {
                    vertexOrderTags.put(vA, "A");
                    vertexOrderTags.put(vB, "B");
                }
                else {
                    vertexOrderTags.put(vA, "B");
                    vertexOrderTags.put(vB, "A");
                }
            }

            recursiveSetUpVertexOrdering(vA, theTree, orderedVertices);
            recursiveSetUpVertexOrdering(vB, theTree, orderedVertices);

            return;
        }
        else {
            /* Degree > 3 */
            ArrayList<HashSet<String>> setsOfTaxa = new ArrayList();
            HashSet<Graph.Vertex> children = theTree.getChildren(v);
            HashMap<HashSet<String>, RootedTree.VertexWithAncestor> vertexFromTaxonSet = new HashMap();
            Iterator it = children.iterator();
            while (it.hasNext()) {
                RootedTree.VertexWithAncestor w = (RootedTree.VertexWithAncestor)it.next();
                HashSet<String> taxonSet = currentTree.getDescendantTaxa(w);
                setsOfTaxa.add(taxonSet);
                vertexFromTaxonSet.put(taxonSet, w);
            }
            Collections.sort(setsOfTaxa, taxonSetComparator);
            findOptimalTaxonSetOrder(setsOfTaxa);
            for (int i=0; i<children.size(); i++) {
                RootedTree.VertexWithAncestor w = vertexFromTaxonSet.get( setsOfTaxa.get(i) );
                vertexOrderTags.put( w , String.valueOf(i));
                recursiveSetUpVertexOrdering(w, theTree, orderedVertices);
            }

        }

    }

    /** Scoring function for two sets of descendant taxa */
    private int scoreVertexPairAgainstRotationalOrder(HashSet<String> v1, HashSet<String> v2, HashSet<String> other) {
        String tA, tB;
        int score = 0;
        for (int i=0; i<rotationalOrder.size(); i++) {

            tA = rotationalOrder.get(i);
            if ((i+1)<rotationalOrder.size()) {
                tB = rotationalOrder.get(i+1);
            }
            else {
                tB = rotationalOrder.get(0);
            }

            if ((v1.contains(tB)&&v2.contains(tA))||(v2.contains(tB)&&other.contains(tA))||(other.contains(tB)&&v1.contains(tA))) {
                // The groups of descendant taxa violate the desired order for tA and tB
                score++;
            }
        }
        return score;
    }

    /** As above, but for a list of taxon sets */
    private int scoreTaxonSetsAgainstRotationalOrder(ArrayList<HashSet<String>> setsOfTaxa) {

        HashSet<String> otherTaxa = new HashSet();
        otherTaxa.addAll(currentTree.getTaxa());
        for (int i=0; i<setsOfTaxa.size(); i++) {
            otherTaxa.removeAll(setsOfTaxa.get(i));
        }

        String tA, tB;
        int score = 0, indA, indB;
        for (int i=0; i<rotationalOrder.size(); i++) {

            tA = rotationalOrder.get(i);
            if ((i+1)<rotationalOrder.size()) {
                tB = rotationalOrder.get(i+1);
            }
            else {
                tB = rotationalOrder.get(0);
            }

            indA = setsOfTaxa.size();
            indB = setsOfTaxa.size();
            for (int j=0; j<setsOfTaxa.size(); j++) {
                if (setsOfTaxa.get(j).contains(tA)) indA = j;
                if (setsOfTaxa.get(j).contains(tB)) indB = j;
            }

            if (!((indA==indB)||(indB==indA+1)||((indA==setsOfTaxa.size())&&(indB==0)))) score++;

        }
        return score;
    }


    /** Compare pairs of taxon sets */
    private class CompareTaxonSetsByMinTaxon implements Comparator<HashSet<String>> {
        public CompareTaxonSetsByMinTaxon() {
        }
        public int compare(HashSet<String> h1, HashSet<String> h2) {
            String s1 = minTaxon(h1);
            String s2 = minTaxon(h2);
            return s1.compareTo(s2);
        }
        private String minTaxon(HashSet<String> h) {
            Iterator<String> it = h.iterator();
            String s = it.next();
            String t;
            while (it.hasNext()) {
                t = it.next();
                if (t.compareTo(s)<0) s=t;
            }
            return s;
        }
    }

    /** Sort an array of taxon sets to match rotational order */
    private void findOptimalTaxonSetOrder(ArrayList<HashSet<String>> setsOfTaxa) {

        int numSets = setsOfTaxa.size();
        ArrayList<HashSet<String>> optimalOrder = new ArrayList();
        optimalOrder.addAll(setsOfTaxa);
        //Collections.shuffle(optimalOrder); // A way to avoid local minima?

        int score = scoreTaxonSetsAgainstRotationalOrder(optimalOrder);
        int oldScore = score+1;
        int testScore, optScore, optPos;
        HashSet<String> taxonSet;

        while (score<oldScore) {

            oldScore = score;

            for (int k=0; k<numSets; k++) {
                taxonSet = setsOfTaxa.get(k);
                optScore = score;
                optPos = optimalOrder.indexOf(taxonSet);
                // Try re-inserting set k in each position
                for (int i=0; i<numSets; i++) {
                    optimalOrder.remove(taxonSet);
                    optimalOrder.add(i, taxonSet);
                    testScore = scoreTaxonSetsAgainstRotationalOrder(optimalOrder);
                    if (testScore<optScore) {
                        optScore = testScore;
                        optPos = i;
                    }
                }
                score = optScore;
                optimalOrder.remove(taxonSet);
                optimalOrder.add(optPos, taxonSet);
                // Completes search for position for taxon k
            }

        }

        setsOfTaxa.clear();
        setsOfTaxa.addAll(optimalOrder);
    }


}