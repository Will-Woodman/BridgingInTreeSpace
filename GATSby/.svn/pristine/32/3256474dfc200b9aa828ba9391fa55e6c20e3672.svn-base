/**
    EmpiricalExchangeabilityMatrix
    Abstract class representing empirical matrix models for amino acids

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
 * Abstract class representing empirical matrix models for amino acids
 **/

public abstract class EmpiricalExchangeabilityMatrixModel implements SubstitutionModel {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Model parameters */
    protected Alphabet theAlphabet = AminoAcidAlphabet.getInstance();
    protected double[][] Q;
    protected double[] stationaryDistrib;
    
    /* For computation of diagonal form */
    protected double[][] U, Uinv; // mtxs of evecs
    protected double[] evals; // row vec of evals    
    
    /** Constructor: fix rate matrix etc.
     p[] must be an array length 20 containing the stationary distribution in order A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V */
    public EmpiricalExchangeabilityMatrixModel() {
        // Fix params
        stationaryDistrib = new double[theAlphabet.size()];
        Q = new double[theAlphabet.size()][theAlphabet.size()];
    }
    
    protected void computeDiagForm() {
        DoubleMatrix2D Qmtx = DoubleFactory2D.dense.make(Q);
        DenseDoubleEigenvalueDecomposition theDiagForm = new DenseDoubleEigenvalueDecomposition(Qmtx);
        DoubleMatrix2D Umtx = theDiagForm.getV();
        U = Umtx.toArray();
        Uinv = cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra.DEFAULT.inverse(Umtx).toArray();
        evals = theDiagForm.getRealEigenvalues().toArray();
    }
    
    @Override
    public Alphabet getAlphabet() {
        return theAlphabet;
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
        return Q[i][j]/stationaryDistrib[j];
    }*/

}
