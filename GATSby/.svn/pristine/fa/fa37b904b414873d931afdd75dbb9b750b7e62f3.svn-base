/**
    GTRModel
    Implementation of the GTR model

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

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;

/**
 * Implementation of the GTR model for a general alphabet
 * See Ziheng Yang, Computational Molecular Evolution, Section 1.5.2
 * Note that for an alphabet of size K rho[K-1,K] is fixed at 1 for identifiability
 */

public class GTRModel implements SubstitutionModel {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /* Model parameters */
    protected double rho[]; // rho is the upper triangle of the exchangeability matrix
                            // ie it contains rho[1,2], ..., rho[1,K], rho[2,1], ..., rho[K-1,K]
    protected final int alphabetSize;
    private Alphabet theAlphabet;
    private double[][] Q;
    private double[] stationaryDistrib;
    
    /* For computation of diagonal form */
    protected double[][] U, Uinv; // mtxs of evecs
    protected double[] evals; // row vec of evals
    
    /** Constructor: fix rate matrix etc.
     p[] must be an array length alphabetSize containing the stationary distribution
     rho[] must be an array of length alphabetSize(alphabetSize-1)/2 - 1 containing the parameters from the upper triangle of the exchangeability
     matrix ignoring the final term (which is fixed at 1), ie rho[1,2], ..., rho[1,K], rho[2,1], ..., rho[K-2,K] */
    public GTRModel(double[] p, double[] r, Alphabet a) {
        // Fix alphabet
        theAlphabet = a;
        alphabetSize = this.theAlphabet.size();
        
        // Fix params
        setParameters(p, r);
    }
    
    @Override
    public Alphabet getAlphabet() {
        return theAlphabet;
    }
    
    protected void setParameters(double[] p, double[] r) {
        // Fix params etc
        stationaryDistrib = new double[alphabetSize];
        rho = new double[(alphabetSize*(alphabetSize-1))/2];
        System.arraycopy(p, 0, stationaryDistrib, 0, p.length);
        System.arraycopy(r, 0, rho, 0, r.length);
        rho[r.length] = 1.0;

        // Fix Q
        int cnt = 0;
        Q = new double[alphabetSize][alphabetSize];
        for(int i=0; i<alphabetSize-1; i++) {
            for(int j=i+1; j<alphabetSize; j++) {
                Q[i][j] = rho[cnt]*stationaryDistrib[j];
                Q[j][i] = rho[cnt]*stationaryDistrib[i];
                Q[i][i] -= Q[i][j];
                Q[j][j] -= Q[j][i];
                cnt++;
            }
        }

        // Normalize to overall rate 1.0
        double fac = StaticLikelihoodCalculator.normalizeQ(Q, stationaryDistrib);
        for(int i=0; i<rho.length; i++) {
            rho[i] /= fac;
        }
 
        // Get diag form
        computeDiagForm();

    }
    
    protected void computeDiagForm() {
        DoubleMatrix2D Qmtx = DoubleFactory2D.dense.make(Q);
        DenseDoubleEigenvalueDecomposition theDiagForm = new DenseDoubleEigenvalueDecomposition(Qmtx);
        DoubleMatrix2D Umtx = theDiagForm.getV();
        U = Umtx.toArray();
        Uinv = cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra.DEFAULT.inverse(Umtx).toArray();
        evals = theDiagForm.getRealEigenvalues().toArray();
    }

    /* For likelihood calculations ------------------------------------------ */
    
    /** Get the transition matrix for a certain branch length */
    public double[][] getTransitionMatrix(double t) {
        int n = theAlphabet.size();
        double[][] res = new double[n][n];
        // Create diagonal mtx
        double[] diag = new double[n];
        for (int k=0; k<n; k++) diag[k] = Math.exp(t*evals[k]);
        // Compute ij-th element of U * diag
        double x;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                // sum U_{ik} U^{-1}_{kj}*D_k
                x = 0.0;
                for (int k=0; k<n; k++) {
                    x += U[i][k]*Uinv[k][j]*diag[k];
                }
                res[i][j] = Math.min(1.0,Math.max(0.0,x));
            }
        }
        return res;
    }
    /** Same method as above, but without creating a new mtx. Saves memory! */
    public void getTransitionMatrix(double t, double[][] theMtx) {
        int n = theAlphabet.size();
        // Create diagonal mtx
        double[] diag = new double[n];
        for (int k=0; k<n; k++) diag[k] = Math.exp(t*evals[k]);
        // Compute ij-th element of U * diag
        double x;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                // sum U_{ik} U^{-1}_{kj}*D_j
                x = 0.0;
                for (int k=0; k<n; k++) {
                    x += U[i][k]*Uinv[k][j]*diag[k];
                }
                theMtx[i][j] = Math.min(1.0,Math.max(0.0,x));
            }
        }
    }
    @Override
    public double[] getStationaryDistrib() {
        double[] x = new double[stationaryDistrib.length];
        System.arraycopy(stationaryDistrib, 0, x, 0, stationaryDistrib.length);
        return x;
    }
    @Override
    public double stationaryProb(int i) {
        return stationaryDistrib[i];
    }
    /*@Override
    public double getQEntry(int i, int j) {
        return Q[i][j];
    }
    @Override
    public double getRhoEntry(int i, int j) throws treebase.AlgorithmException {
        if(i==j) {
            throw new treebase.AlgorithmException("AlgorithmExcepion: diagonal elements in rho matrix not defined.");
        }
        int max = (i > j) ? i : j;
        int min = (i < j) ? i : j;
        
        int cnt = -1;
        for(int k=0; k<min; k++) {
            for(int l=k+1; l<alphabetSize; l++) {
                cnt++;
            }
        }
        for(int l=min+1; l<=max; l++) {
            cnt++;
        }
        return rho[cnt];
    }*/
    
    
}