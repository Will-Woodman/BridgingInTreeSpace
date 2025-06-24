/**
    PoissonModel
    Implementation of the Poisson model for amino acids

    Copyright (C) 2012  Sarah E. Heaps

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

package sequencemodels;

/**
 * Implementation of the Poisson model for amino acids, ie the empirical matrix for
 * amino acids exchangeabilities with 1's everywhere off the diagonal.
 */

public class PoissonModel extends EmpiricalExchangeabilityMatrixModel {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    /** Constructor: fix rate matrix etc.
     p[] must be an array length 20 containing the stationary distribution in order A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V */
    public PoissonModel(double[] p) {
        // Fix alphabet
        super();
        // Fix params
        setParameters(p);
    }
    
    protected void setParameters(double[] p) {
        int numAminoAcids = theAlphabet.size();
        
         // Fix params etc
        System.arraycopy(p, 0, stationaryDistrib, 0, numAminoAcids);

        // Fix Q
        for(int i=0; i<numAminoAcids; i++) {
            for(int j=0; j<numAminoAcids; j++) {
                if(i!=j) Q[i][j] = p[j];
                else Q[i][j] = p[j] - 1.0;
            }
        }

        // Normalize to overall rate 1.0
        StaticLikelihoodCalculator.normalizeQ(Q, stationaryDistrib);

    }
    
    /* For likelihood calculations ------------------------------------------ */

    public double[][] getTransitionMatrix(double t) {

        double sum = 0.0;
        for(int i=0; i<20; i++) sum += (1.0-stationaryDistrib[i])*stationaryDistrib[i];
        double lambda = 1.0/sum;
        
        double[][] P = new double[20][20];
        for(int i=0; i<20; i++) {
            for(int j=0; j<20; j++) {
                if(i==j) P[i][j] = Math.min(Math.max(0.0, stationaryDistrib[i]+(1-stationaryDistrib[i])*Math.exp(-lambda*t)),1.0);
                else P[i][j] = Math.min(Math.max(0.0, stationaryDistrib[j]*(1-Math.exp(-lambda*t))),1.0);
            }
        }
        return P;
    }
    public void getTransitionMatrix(double t, double[][] P) {
        
        double sum = 0.0;
        for(int i=0; i<20; i++) sum += (1.0-stationaryDistrib[i])*stationaryDistrib[i];
        double lambda = 1.0/sum;
        
        for(int i=0; i<20; i++) {
            for(int j=0; j<20; j++) {
                if(i==j) P[i][j] = Math.min(Math.max(0.0, stationaryDistrib[i]+(1-stationaryDistrib[i])*Math.exp(-lambda*t)),1.0);
                else P[i][j] = Math.min(Math.max(0.0, stationaryDistrib[j]*(1-Math.exp(-lambda*t))),1.0);
            }
        }
    }
    

    
}
