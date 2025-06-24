/**
 * EKMaxFlow
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

package geodesics;

/**
 * Implementation of Edmonds-Karp max flow algorithm, specifically on a 
 * bipartite graph with weighted vertices where vertex is a <Split> object and 
 * the bipartite graph represents incompatibility between splits. 
 * The max flow corresponds to a minimum weight vertex cover. 
 * See Owen 2010 paper on polynomial time algorithm for calculation of geodesic 
 * between two trees. 
 * 
 * Optimization code taken directly from Wikipedia.
 *
 * Usage:
 * - Call the constructor, which also computes the max flow / min cover
 * - Call getVertexPartitions to extract the sets required for the refinement
 *   of a partition in Owen's algorithm. 
 */

import treebase.*;
import java.util.*;

public class EKMaxFlow {

    public static final boolean DEBUG_ON = true;

    public static final int WHITE = 0, GRAY = 1, BLACK = 2;
    private double[][] flow, capacity, res_capacity;
    private int[] parent, color, queue;
    private double[] min_capacity;
    private int size, source, sink, first, last;
    public double max_flow;
    private boolean[] S;
    private double[] weightsA, weightsB;

    private HashSet coverVerticesA, coverVerticesB;
    private ArrayList splitsA, splitsB;


    /** Constructor: get up the bipartite graph from two sets of weighted splits */
    public EKMaxFlow(ArrayList partA, ArrayList partB, double[] wA, double[] wB) {
        
        // Store splits
        splitsA = new ArrayList(partA);
        splitsB = new ArrayList(partB);
        weightsA = wA;
        weightsB = wB;
        
        /* Initialize incompatibility graph*/
        size = splitsA.size()+splitsB.size()+2;
        source = 0; // source is vertex 0
        sink = size-1; // sink is last vertex
        capacity = new double[size][size];
        int numEdges = 0;

        // Set up capacities
        // 1. From source
        for (int i=1; i<=splitsA.size(); i++) {
            capacity[source][i] = weightsA[(i-1)];
        }
        // 2. To sink
        for (int j=1; j<=splitsB.size(); j++) {
            capacity[splitsA.size()+j][sink] = weightsB[(j-1)];
        }
        // 3. From A to B
        Split sA, sB;
        ListIterator itA, itB;
        itA = splitsA.listIterator();
        for (int i=1; i<=splitsA.size(); i++) {
            sA = (Split) itA.next();
            itB = splitsB.listIterator();
            for (int j=1; j<=splitsB.size(); j++) {
                sB = (Split) itB.next();
                if (!(sA.isCompatible(sB))) {
                    capacity[i][splitsA.size()+j] = 1.0; 
                    numEdges++;
                }
            }
        }

        if (numEdges==0) {
            max_flow = 1.0;
            return;
        }

        // Run maximization
        maxFlow();
        
    } 

    /** Get HashSets corresponding to the partition of vertices induced by the flow */
    public HashSet getACover() {
        return coverVerticesA;
    }
    public HashSet getBCover() {
        return coverVerticesB;
    }


    /* ------------------------------------------------------------------ */

    /* Algorithm copied from wikipedia! */

    private void maxFlow()  // Edmonds-Karp algorithm with O(V³E) complexity
    {
            flow = new double[size][size];
            res_capacity = new double[size][size];
            parent = new int[size];
            min_capacity = new double[size];
            color = new int[size];
            queue = new int[size];
            S = new boolean[size];

            for (int i = 0; i < size; i++)
                    for (int j = 0; j < size; j++)
                            res_capacity[i][j] = capacity[i][j];


            while (BFS(source))
            {
                    max_flow += min_capacity[sink];
                    int v = sink, u;
                    while (v != source)
                    {
                            u = parent[v];
                            flow[u][v] += min_capacity[sink];
                            flow[v][u] -= min_capacity[sink];
                            res_capacity[u][v] -= min_capacity[sink];
                            res_capacity[v][u] += min_capacity[sink];
                            v = u;
                    }
            }

            // Prepare the cover!!!

            coverVerticesA = new HashSet();
            coverVerticesB = new HashSet();

            for (int i=1; i<=splitsA.size(); i++) {
                if (!S[i]) coverVerticesA.add(splitsA.get(i-1));
            }

            for (int j=1; j<=splitsB.size(); j++) {
                if (S[splitsA.size()+j]) coverVerticesB.add(splitsB.get(j-1));
            }

            if (coverVerticesA.size()==splitsA.size()) {
                max_flow = 1.0;
            }
            if (coverVerticesB.size()==splitsB.size()) {
                max_flow = 1.0;
            }


            if (DEBUG_ON) {
                double sAbar = 0.0;
                for (int i=1; i<=splitsA.size(); i++) {
                    if (!S[i]) sAbar += weightsA[i-1];
                }
                // System.out.println("Weight bar{S_A} = "+sAbar);

                double sB = 0.0;
                for (int j=1; j<=splitsB.size(); j++) {
                    if (S[splitsA.size()+j]) sB += weightsB[j-1];
                }
                //System.out.println("Weight S_B = "+sB);
                
                // Sanity check
                if (Math.abs(sAbar+sB-max_flow)>1.0E-6) {
                    System.out.println("Error: max flow != min weight cover");
                }
            }

    }

    private boolean BFS(int source)  // Breadth First Search in O(V<sup>2</sup>)
    {
            for (int i = 0; i < size; i++)
            {
                    S[i] = false;
                    color[i] = WHITE;
                    min_capacity[i] = Double.MAX_VALUE;
            }
            S[source] = true;

            first = last = 0;
            queue[last++] = source;
            color[source] = GRAY;

            while (first != last)  // While "queue" not empty..
            {
                    int v = queue[first++];
                    for (int u = 0; u < size; u++)
                            if (color[u] == WHITE && res_capacity[v][u] > 0)
                            {
                                    min_capacity[u] = Math.min(min_capacity[v], res_capacity[v][u]);
                                    parent[u] = v;
                                    color[u] = GRAY;
                                    if (u == sink) return true;
                                    queue[last++] = u;
                                    S[u] = true;
                            }
            }
            return false;
    }

}