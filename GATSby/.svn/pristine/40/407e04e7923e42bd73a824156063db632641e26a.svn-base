/*
 * TopologicalOperations
   Copyright (C) 2013  Sarah E. Heaps

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

    Contact the author at:  <sarah.heaps@ncl.ac.uk>
 */

package treebase;

import treebase.Graph.Vertex;
import treebase.Graph.Edge;

/**
 * Interface representing the core topological operations
 */

public interface TopologicalOperations {
 
    public Edge performNNI(Edge e, Vertex[][] v, double newLength, boolean choice) throws AlgorithmException;
    public Edge reverseNNI(Edge oldEdge, Edge newEdge, Vertex[][] v, boolean choice) throws AlgorithmException;
    public Edge[][] performSPR(Edge pruneEdge, Edge graftEdge, double lengthMergedEdge, double lengthSplitEdge1, double lengthSplitEdge2, double lengthPrunedEdge) throws AlgorithmError;
    public Edge[] reverseSPR(Edge[][] edgeArray) throws AlgorithmError;
    public Edge[][] performRootMove(Edge newRootEdge, double lengthMergedEdge, double lengthSplitEdge1, double lengthSplitEdge2) throws AlgorithmError;
    public Edge[] reverseRootMove(Edge[][] edgeArray) throws AlgorithmError;

}
