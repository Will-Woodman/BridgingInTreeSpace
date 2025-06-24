/*
    TreeSpacePath
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

package geodesicgraphics;

/**
 * Interface for obtaining trees on paths through tree space
 */

public interface TreeSpacePath {

    /** Return the tree at a particular point */
    public treebase.TreeAsSplits getTreeOnPath(double s);

    /** Get bounds on the parameter */
    public double[] getPathBounds();

    /** Get a rotational order */
    public java.util.ArrayList<String> getRotationalOrder();

    /** Get a shared root split */
    public treebase.Split getRootSplit();

}