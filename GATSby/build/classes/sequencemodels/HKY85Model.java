/*
    HKY85Model
    Implementation of the HKY85 model for DNA substitution

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

package sequencemodels;

/**
 * Implementation of the HKY85 model for DNA substitution
 * See Ziheng Yang, pages 15-18
 */

public class HKY85Model extends TN93Model {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /** Constructor: fix rate matrix etc.
     p[] must be an array length 4 containing the stationary distribution in order A, G, C, T
     kappa is the transition / transversion ratio*/
    public HKY85Model(double[] p, double kappa) {
        super(p, kappa, kappa);
    }

    public void changeParameters(double[] p, double kappa) {
        setParameters(p, kappa, kappa);
    }

}