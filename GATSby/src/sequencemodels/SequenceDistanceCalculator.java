/** 
    DistanceCalculator
    Methods for calculating distances between sequences.

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

package sequencemodels;

/**
 * Methods for calculating distances between sequences.
 * Not all sequence models will implement these -- hence inclusion in a
 * separate interface. 
 */

public interface SequenceDistanceCalculator {

    /** Distance between two sequences */
    public double calcDistance(int[] a, int[] b, double infiniteDist) throws treebase.AlgorithmError;

    /** Calc array of distances from an alignment */
    public double[][] calcDistanceMatrix(Alignment a, double infiniteDist);

}
