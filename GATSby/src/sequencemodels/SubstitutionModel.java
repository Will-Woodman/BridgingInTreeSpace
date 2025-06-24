/**
    SubstitutionModel
    Interface representing a Markov substitution model for sequence evolution

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
 * Interface representing a Markov substitution model for sequence evolution
 */


public interface SubstitutionModel extends java.io.Serializable {
    
    /* Utilities ------------------------------------------------------------ */

    public Alphabet getAlphabet();
    
    
    /* Rate and transition calculations ------------------------------------- */
    
    public double[] getStationaryDistrib();

    public double stationaryProb(int i);
    
    /** Get the transition matrix for a certain branch length */
    public double[][] getTransitionMatrix(double t);
    /** Same method as above, but without creating a new mtx. Saves memory! */
    public void getTransitionMatrix(double t, double[][] theMtx);

}

