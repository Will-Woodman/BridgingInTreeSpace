/**
    DrawTree
    Interface enabling tree to be drawn in a TreeViewer

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

/**
 * Interface enabling tree to be drawn in a TreeViewer
 */


public interface DrawTree {

    /** Draw the tree */
    public void draw(java.awt.Graphics g, TreeViewer theViewer, java.awt.Dimension treeArea);

    /** Calculate the size of the tree without drawing it */
    public java.awt.Dimension calculateSize(TreeViewer theViewer, java.awt.Graphics g);

    /** Return the Vertex at a given coordinate */
    public treebase.Graph.Vertex getVertexAtPoint(java.awt.Point where);

    /** Return the Edge at a given coordinate */
    public treebase.Graph.Edge getEdgeAtPoint(java.awt.Point where);

    /** Measure the max path length */
    public double getMaxPathLength();

    /** Get tooltip labels */
    public String getEdgeLabel(treebase.Graph.Edge e);
    public String getVertexLabel(treebase.Graph.Vertex v);

}
