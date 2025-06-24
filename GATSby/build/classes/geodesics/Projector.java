/*
 * Projector.java

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

package geodesics;

/**
 * Methods dealing with projection onto a geodesic segment
 */

import java.util.ArrayList;
import treebase.TreeAsSplits;
import treebase.AlgorithmException;

public class Projector {

    private double TOLERANCE_ON_PROJ = 0.0005; // Tolerance on projection, as a proportion 0 to 1.
    private boolean IGNORE_PENDANT_EDGES = true;


    /** Constructor: set tolerance */
    public Projector(double tol, boolean leavesIgnored) {
        TOLERANCE_ON_PROJ = tol;
        IGNORE_PENDANT_EDGES = leavesIgnored;
    }

    /** Set convergence control for projection */
    public void setTolerance(double tol) {
        TOLERANCE_ON_PROJ = tol;
    }
    
    /** Get tolerance */
    public double getTolerance() {
        return TOLERANCE_ON_PROJ;
    }

    public void ignorePendantEdges() {
        IGNORE_PENDANT_EDGES = true;
    }
    public void includePendantEdges() {
        IGNORE_PENDANT_EDGES = false;
    }
    public boolean arePendantEdgesIgnored() {
        return IGNORE_PENDANT_EDGES;
    }

    /** Inner class: objective function to minimise orthogonal distance */
    private class GeodesicObjective implements Optimisation1D.Objective {

        private Geodesic theLine;
        private TreeAsSplits projTree;
        public Geodesic perpLine;

        public GeodesicObjective(Geodesic l, TreeAsSplits t) {
            theLine = l;
            projTree = t;
        }

        /** Evaluate at x\in[0,1] */
        public double evaluate(double x) throws AlgorithmException {
            TreeAsSplits treeOnLine = theLine.getTree(x);
            if (perpLine!=null) {
                 Geodesic h = perpLine.shiftEndPoints(treeOnLine, projTree);
                 perpLine = h;
            }
            else {
                perpLine = new Geodesic(treeOnLine, projTree);
            }
            if (IGNORE_PENDANT_EDGES) return perpLine.getInternalLength();
            else return perpLine.getLength();
        }

        public void setTree(TreeAsSplits t) {
            projTree = t;
            perpLine = null;
        }
    }

    /** Project a tree onto a geodesic segment.
     Output: double[index]
     index = 0 is perpendicular distance
     index = 1 is geodesic distance from end of theLine
     index = 2 is an indicator. Value 0 means in interior of line, -1 or +1 correspond to the ends
    */
    public double[] projectOntoLineSegment(Geodesic theLine, TreeAsSplits theTree) {

        GeodesicObjective theObjective = new GeodesicObjective(theLine, theTree);

        double len;
        if (IGNORE_PENDANT_EDGES) len = theLine.getInternalLength();
        else len = theLine.getLength();

        double[] dist;
        try
        {
            dist = Optimisation1D.optimise(theObjective, new double[]{0.0, 1.0}, TOLERANCE_ON_PROJ);
        }
        catch (treebase.AlgorithmError anErr) {
            System.out.print("WARNING: error occurred while computing a geodesic for projection: "+anErr.getMessage());
            System.out.println(" This should not be possible.");
            return new double[]{0.0, 0.0, 0.0};
        }
        catch (treebase.AlgorithmException anEx) {
            System.out.print("WARNING: error occurred evaluating geodesic objective: "+anEx.getMessage());
            System.out.println(" This should not be possible.");
            return new double[]{0.0, 0.0, 0.0};
        }


        double[] res = new double[3];

        // Convert [0,1] coordinate on segment to geodesic distance
        res[1] = dist[0]*len;

        // Perp distance
        res[0] = dist[1];

        // Detect projection onto end points
        if (Math.abs(dist[0])<2.0*TOLERANCE_ON_PROJ)
            res[2]=-1.0;
        else if (Math.abs(dist[0]-1.0)<2.0*TOLERANCE_ON_PROJ)
            res[2]=1.0;
        else
            res[2]=0.0;

        return res;
    }

    /** Brute force project */
    public double[] projectionCheck(Geodesic theLine, TreeAsSplits theTree, int n) {
        double[] res = new double[3];
        double y;

        try {
            for (int i=0; i<=n; i++) {
                double s = i/(1.0*n);
                TreeAsSplits treeOnLine = theLine.getTree(s);
                Geodesic h = new Geodesic(theTree, treeOnLine);
                if (IGNORE_PENDANT_EDGES) y = h.getInternalLength();
                else y = h.getLength();
                if ((i==0)||(y<res[0])) {
                    res[1] = s*theLine.getLength();
                    res[0] = y;
                }
            }
        }
        catch (AlgorithmException anErr) {
            System.out.println("Error building geodesic during projection debug.");
        }
        return res;
    }


    /** Project many trees onto a geodesic segment.
     Output double[index][tree number+1]

     numTrees+1 is a summary: sums of squares for index 0,1 and an indicator for index 2
     (zero means all trees projected onto the interior)

     index = 0 is perpendicular distance
     index = 1 is geodesic distance from end of theLine
     index = 2 is an indicator. Value 0 means in interior of line, -1 or +1 correspond to the ends.
    */
    public double[][] projectOntoLineSegment(Geodesic theLine, ArrayList<TreeAsSplits> theTrees, boolean adjust) {

        GeodesicObjective theObjective = new GeodesicObjective(theLine, null);

        // Get length of line segment
        double len;
        if (IGNORE_PENDANT_EDGES) len = theLine.getInternalLength();
        else len = theLine.getLength();

        // Set up storage for output
        double[][] res = new double[3][theTrees.size()+1];
        double interiorCount = 0.0;

        // Loop thru' trees
        double perp = 0.0;
        double min = len, max = 0.0;
        for (int i=0; i<theTrees.size(); i++) {
            TreeAsSplits tree = theTrees.get(i);
            theObjective.setTree(tree);

            double[] dist;
            try
            {
                dist = Optimisation1D.optimise(theObjective, new double[]{0.0, 1.0}, TOLERANCE_ON_PROJ);
            }
            catch (treebase.AlgorithmError anErr) {
                System.out.print("WARNING: error occurred while computing a geodesic for projection: "+anErr.getMessage());
                System.out.println(" This should not be possible.");
                dist = new double[]{0.0, 0.0, 0.0};
            }
            catch (treebase.AlgorithmException anEx) {
                System.out.print("WARNING: error occurred evaluating geodesic objective: "+anEx.getMessage());
                System.out.println(" This should not be possible.");
                dist =  new double[]{0.0, 0.0, 0.0};
            }

            // Convert [0,1] coordinate on segment to geodesic distance
            res[1][i] = dist[0]*len;
            if (res[1][i]<min) min=res[1][i];
            if (res[1][i]>max) max=res[1][i];

            // Perp distance
            res[0][i] = dist[1];

            // Detect projection onto end points
            if (Math.abs(dist[0])<2.0*TOLERANCE_ON_PROJ) {
                res[2][i]=-1.0;
                interiorCount+=1.0;
            }
            else if (Math.abs(dist[0]-1.0)<2.0*TOLERANCE_ON_PROJ) {
                res[2][i]=1.0;
                interiorCount=+1.0;
            }
            else {
                res[2][i]=0.0;
            }

            perp += res[0][i]*res[0][i];

        } // End loop thru' trees

        if (adjust) {
            perp += min*min;
            perp += (max-len)*(max-len);
        }
        res[0][theTrees.size()] = perp;
        res[2][theTrees.size()] = interiorCount;
        /* Compute parallel sum squares */
        double mean = 0.0;
        for (int i=0; i<theTrees.size(); i++) {
            mean += res[1][i];
        }
        mean = mean/theTrees.size();
        double para = 0.0;
        for (int i=0; i<theTrees.size(); i++) {
            double z = res[1][i]-mean;
            para += z*z;
        }
        res[1][theTrees.size()] = para;

        return res;

    }

    /** Project but just return the perp sum of squares */
    public double getPerpSumSquares(Geodesic theLine, ArrayList<TreeAsSplits> theTrees, boolean adjust, int numProcessors) {
        if (numProcessors<2) {
            return getPerpSumSquares(theLine, theTrees, adjust);
        }
        int chunk = Math.round(((float)theTrees.size())/((float)numProcessors));
        int start=0, end=chunk;
        double[][] res = new double[3][theTrees.size()+1];

        ProjectionChunk[] processes = new ProjectionChunk[numProcessors];
        for (int i=0; i<numProcessors; i++) {
            // Create processes
            if (i==(numProcessors-1)) end=theTrees.size();
            ProjectionChunk p = new ProjectionChunk(theLine, theTrees, start, end, res);
            /* Question: would it be better to use a different arraylist of trees for each process? */
            p.start();
            processes[i] = p;
            if (i<(numProcessors-1)) {
                chunk = Math.round(((float)theTrees.size()-end)/((float)numProcessors-i-1));
                start = end;
                end = start+chunk;
            }
        }
        // Wait
        double perp=0.0;
        double len = (IGNORE_PENDANT_EDGES) ? theLine.getInternalLength() : theLine.getLength();
        double min=len, max=0.0;
        try {
             for (int i=0; i<numProcessors; i++) {
                processes[i].join();
                perp += processes[i].getPerpDistance();
                if (processes[i].min<min) min = processes[i].min;
                if (processes[i].max>max) max = processes[i].max;
            }
        }
        catch (InterruptedException ex) {
           System.out.println("Error: something funny happened with a thread interrupt while projecting trees.");
        }
        if (adjust) {
            perp += min*min + (max-len)*(max-len);
        }

        return perp;

    }
    /** Project but just return the perp sum of squares.
     Adjust=true means add on end segment lengths. */
    public double getPerpSumSquares(Geodesic theLine, ArrayList<TreeAsSplits> theTrees, boolean adjust) {

        GeodesicObjective theObjective = new GeodesicObjective(theLine, null);

        // Get length of line segment
        double len = (IGNORE_PENDANT_EDGES) ? theLine.getInternalLength() : theLine.getLength();

        // Loop thru' trees
        double perp = 0.0;
        double min=len, max=0.0;
        for (int i=0; i<theTrees.size(); i++) {
            TreeAsSplits tree = theTrees.get(i);
            theObjective.setTree(tree);

            double[] dist;
            try
            {
                dist = Optimisation1D.optimise(theObjective, new double[]{0.0, 1.0}, TOLERANCE_ON_PROJ);
            }
            catch (treebase.AlgorithmError anErr) {
                System.out.print("WARNING: error occurred while computing a geodesic for projection: "+anErr.getMessage());
                System.out.println(" This should not be possible.");
                dist = new double[]{0.0, 0.0, 0.0};
            }
            catch (treebase.AlgorithmException anEx) {
                System.out.print("WARNING: error occurred evaluating geodesic objective: "+anEx.getMessage());
                System.out.println(" This should not be possible.");
                dist =  new double[]{0.0, 0.0, 0.0};
            }

            perp += dist[1]*dist[1];

            if (adjust) {
                if (dist[0]*len>max) max = dist[0]*len;
                if (dist[0]*len<min) min = dist[0]*len;
            }

        } // End loop thru' trees

        if (adjust) {
            perp += min*min + (max-len)*(max-len);
        }

        return perp;

    }

    /** Project many trees onto a geodesic segment: parallel processing! */

    public double[][] projectOntoLineSegmentParallel(Geodesic theLine, ArrayList<TreeAsSplits> theTrees, boolean adjust, int numProcessors)  {
        int chunk = Math.round(((float)theTrees.size())/((float)numProcessors));
        int start=0, end=chunk;
        double[][] res = new double[3][theTrees.size()+1];

        ProjectionChunk[] processes = new ProjectionChunk[numProcessors];
        for (int i=0; i<numProcessors; i++) {
            // Create processes
            if (i==(numProcessors-1)) end=theTrees.size();
            ProjectionChunk p = new ProjectionChunk(theLine, theTrees, start, end, res);
            /* Question: would it be better to use a different arraylist of trees for each process? */
            p.start();
            processes[i] = p;
            if (i<(numProcessors-1)) {
                chunk = Math.round(((float)theTrees.size()-end)/((float)numProcessors-i-1));
                start = end;
                end = start+chunk;
            }
        }
        // Wait
        double min=0.0, max=0.0, perp=0.0;
        int numOnEnds = 0;
        try {
             for (int i=0; i<numProcessors; i++) {
                processes[i].join();
                if (i==0) {
                    min=processes[0].min; max=processes[0].max;
                }
                else {
                    if (processes[i].min<min) min=processes[i].min;
                    if (processes[i].max>max) max=processes[i].max;
                }
                perp += processes[i].getPerpDistance();
                numOnEnds += processes[i].getNumOnEnds();
             }
             res[0][theTrees.size()] = perp;
             res[2][theTrees.size()] = 1.0*numOnEnds;
        }
        catch (InterruptedException ex) {
           System.out.println("Error: something funny happened with a thread interrupt while projecting trees.");
        }

        if (adjust) {
            double len = processes[0].len;
            res[0][theTrees.size()] += min*min;
            res[0][theTrees.size()] += (max-len)*(max-len);
        }

        /* Compute parallel sum squares */
        double mean = 0.0;
        for (int i=0; i<theTrees.size(); i++) {
            mean += res[1][i];
        }
        mean = mean/theTrees.size();
        double para = 0.0;
        for (int i=0; i<theTrees.size(); i++) {
            double z = res[1][i]-mean;
            para += z*z;
        }
        res[1][theTrees.size()] = para;

        return res;

    }

    /** Inner class for calculation */
    public class ProjectionChunk extends Thread {

        private int start, end;
        private Geodesic theLine;
        private ArrayList<TreeAsSplits> theTrees;
        private double[][] res; // Results: null if all you want is the sum of square perp distances
        private double perp;
        public double min, max, len;
        private int numOnEnds;

        public ProjectionChunk(Geodesic g, ArrayList<TreeAsSplits> a, int s, int e, double[][] r) {
            theLine = g;
            theTrees = a;
            start=s;
            end=e;
            res = r;
            len = (IGNORE_PENDANT_EDGES) ? theLine.getInternalLength() : theLine.getLength();
            max = 0.0;
            min = len;
            numOnEnds = 0;
        }

        public ProjectionChunk(Geodesic g, ArrayList<TreeAsSplits> a, int s, int e) {
            theLine = g;
            theTrees = a;
            start=s;
            end=e;
            res = null;
            len = (IGNORE_PENDANT_EDGES) ? theLine.getInternalLength() : theLine.getLength();
            max = 0.0;
            min = len;
            numOnEnds = 0;
        }

        public void run() {
            GeodesicObjective theObjective = new GeodesicObjective(theLine, null);


            // Loop thru' trees
            for (int i=start; i<end; i++) {
                TreeAsSplits tree = theTrees.get(i);
                theObjective.setTree(tree);

                double[] dist;
                try
                {
                    dist = Optimisation1D.optimise(theObjective, new double[]{0.0, 1.0}, TOLERANCE_ON_PROJ);
                }
                catch (treebase.AlgorithmError anErr) {
                    System.out.print("WARNING: error occurred while computing a geodesic for projection: "+anErr.getMessage());
                    System.out.println(" This should not be possible.");
                    dist = new double[]{0.0, 0.0, 0.0};
                }
                catch (treebase.AlgorithmException anEx) {
                    System.out.print("WARNING: error occurred evaluating geodesic objective: "+anEx.getMessage());
                    System.out.println(" This should not be possible.");
                    dist =  new double[]{0.0, 0.0, 0.0};
                }

                if (res!=null) {
                    // Convert [0,1] coordinate on segment to geodesic distance
                    res[1][i] = dist[0]*len;

                    // Perp distance
                    res[0][i] = dist[1];

                    // Detect projection onto end points
                    if (Math.abs(dist[0])<2.0*TOLERANCE_ON_PROJ) {
                        res[2][i]=-1.0;
                        numOnEnds++;
                    }
                    else if (Math.abs(dist[0]-1.0)<2.0*TOLERANCE_ON_PROJ) {
                        res[2][i]=1.0;
                        numOnEnds++;
                    }
                    else {
                        res[2][i]=0.0;
                    }

                    perp += res[0][i]*res[0][i];
                }
                else { // res==null, so we just want a perp distance summary
                    perp += dist[1]*dist[1];
                }

                if (dist[0]*len>max) max = dist[0]*len;
                if (dist[0]*len<min) min = dist[0]*len;

            } // End loop thru' trees
        }

        public double getPerpDistance() {
            return perp;
        }

        public int getNumOnEnds() {
            return numOnEnds;
        }

    } // End inner class
    
}