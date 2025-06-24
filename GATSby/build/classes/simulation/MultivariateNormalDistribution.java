/**
 *  MultivariateNormalDistribution
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

package simulation;

/**
 * Multivariate normal distribution.
 *
 */

import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import java.util.Arrays;
import treebase.AlgorithmException;

public class MultivariateNormalDistribution implements java.io.Serializable, RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;
    
    /** Instance variables */
    private final int dimension;
    private double[] mu;
    private SymmetricPositiveDefiniteMatrix sigmaOrPrecision;/* Variance can be summarised by
                                            variance OR precision matrix. There
                                            is a time-saving in density evaluation
                                            if a precision matrix is used.*/
    private boolean isPrecision;
    public boolean performChecks;
    private Normal coltObject;
    private double[] workingVector;
    private double[] workingVector2;
    
    /** Default constructor */
    public MultivariateNormalDistribution(int n, boolean performChecks) throws AlgorithmException {
        dimension = n;
        mu = new double[n];
        Arrays.fill(mu, 0.0);
        sigmaOrPrecision = new SymmetricPositiveDefiniteMatrix(n, performChecks);
        isPrecision = true;
        coltObject = new Normal(0.0, 1.0, Random.getEngine());
        workingVector = new double[dimension]; 
        workingVector2 = new double[dimension];
        this.performChecks = performChecks;
    }
    
    public MultivariateNormalDistribution(double[] mean, double[][] varOrPrec, boolean isPrec, boolean performChecks) throws AlgorithmException {
        this(mean, varOrPrec, isPrec, performChecks, Random.getEngine());
    }
    
    /** Constructors which allow use of random number generators which are not
     the global generator*/
    public MultivariateNormalDistribution(double[] mean, double[][] varOrPrec, boolean isPrec, boolean performChecks, int seed) throws AlgorithmException {
        this(mean, varOrPrec, isPrec, performChecks, new DoubleMersenneTwister(seed));
    } 
    public MultivariateNormalDistribution(double[] mean, double[][] varOrPrec, boolean isPrec, boolean performChecks, DoubleMersenneTwister randomEngine) throws AlgorithmException {
        dimension = mean.length;
        mu = new double[dimension];
        System.arraycopy(mean, 0, mu, 0, dimension);
        sigmaOrPrecision = new SymmetricPositiveDefiniteMatrix(varOrPrec, performChecks);
        isPrecision = isPrec;
        coltObject = new Normal(0.0, 1.0, randomEngine);
        workingVector = new double[dimension]; 
        workingVector2 = new double[dimension];
        this.performChecks = performChecks;
    }
    
    public void resetRandomEngineSeed() {
        coltObject = new Normal(0.0, 1.0, Random.getEngine());
    }
    
    public void setMean(double[] mean) throws AlgorithmException {
        if(performChecks && mean.length!=dimension) throw new AlgorithmException("Dimension of new mean is not correct.");
        System.arraycopy(mean, 0, mu, 0, dimension);
    }
    
    public int getDimension() {
        return dimension;
    }
    
    public double[] getMean() {
        double[] mean = new double[dimension];
        System.arraycopy(mu, 0, mean, 0, dimension);
        return mean;
    }
    
    public void setSigma(double[][] variance) throws AlgorithmException {
        sigmaOrPrecision.setMatrix(variance);
        isPrecision = false;
    }
    public void setSigma(SymmetricPositiveDefiniteMatrix variance) throws AlgorithmException {
        sigmaOrPrecision.setMatrix(variance);
        isPrecision = false;
    }
    
    public void setPrecision(double[][] precision) throws AlgorithmException {
        sigmaOrPrecision.setMatrix(precision);
        isPrecision = true;
    }
    public void setPrecision(SymmetricPositiveDefiniteMatrix precision) throws AlgorithmException {
        sigmaOrPrecision.setMatrix(precision);
        isPrecision = true;
    }
    
    public double[][] getSigmaValue() throws AlgorithmException {
        //AVOID using this if sigmaOrPrecision is a precision matrix
        if(isPrecision){
            return sigmaOrPrecision.getInverseMatrix();
        } else {
            return sigmaOrPrecision.getMatrix();
        }
    }
    
    public double[][] getPrecisionValue() throws AlgorithmException {
        //AVOID using this if sigmaOrPrecision is a variance matrix
        if(!isPrecision){
            return sigmaOrPrecision.getInverseMatrix();
        } else {
            return sigmaOrPrecision.getMatrix();
        }
    }
    
    /** Use with caution */
    public SymmetricPositiveDefiniteMatrix getSigmaOrPrecision() throws AlgorithmException {
        return sigmaOrPrecision;
    }
    
    public boolean getIsPrecision() {
        return isPrecision;
    }
    
    public double logpdf(double[] xArr) throws AlgorithmException {
        int i;
        double ret=-(dimension/2.0)*Math.log(2.0*Math.PI),logDet=sigmaOrPrecision.getLogDeterminant();
        for(i=0; i<dimension; i++) {
            workingVector[i] = xArr[i] - mu[i];
        }
        if(isPrecision) {
            logDet *= -1.0;
            sigmaOrPrecision.upperTriangleMult(workingVector, workingVector2);
        } else {
            sigmaOrPrecision.lowerTriangleSolve(workingVector, workingVector2);
        }
        ret-=0.5*logDet;
        for(i=0; i<dimension; i++) {
            ret -= 0.5*workingVector2[i]*workingVector2[i]; // Actually could lose workingVector here
        }
        return ret;
    }
    
    public double pdf(double[] xArr) throws AlgorithmException {
        return Math.exp(this.logpdf(xArr));
    }
    
    public double[] sample() throws AlgorithmException {
        int i;
        for(i=0;i<dimension;i++)
            workingVector[i] = coltObject.nextDouble();
        if(isPrecision){
            sigmaOrPrecision.upperTriangleSolve(workingVector, workingVector2);
        } else {
            sigmaOrPrecision.lowerTriangleMult(workingVector, workingVector2);
        }
        double[] xArr = new double[dimension];
        for(i=0;i<dimension;i++)
           xArr[i]=workingVector2[i]+mu[i];
        return xArr;
    }
    
    public static class SymmetricPositiveDefiniteMatrix implements java.io.Serializable {
        
        private int dimension;
        protected double[][] theMatrix;
        private CholeskyDecomposition decompObject;
        private boolean performChecks;
        
        public SymmetricPositiveDefiniteMatrix(double[][] mat, boolean performChecks) throws AlgorithmException {
            dimension = mat.length;
            theMatrix = new double[dimension][dimension];
            checkSizeSymmetry(mat);
            for (int i = 0; i < dimension; i++) {
                System.arraycopy(mat[i], 0, theMatrix[i], 0, dimension);
            }
            decompObject = new CholeskyDecomposition(dimension, performChecks);
            decompObject.factor(theMatrix);
            this.performChecks = performChecks;
        }
        
        public SymmetricPositiveDefiniteMatrix(int size, boolean performChecks) throws AlgorithmException {
            dimension = size;
            theMatrix = new double[dimension][dimension];
            for(int i=0; i<dimension; i++) {
                Arrays.fill(theMatrix[i], 0.0);
                theMatrix[i][i] = 1.0;
            }
            decompObject = new CholeskyDecomposition(dimension, performChecks);
            decompObject.factor(theMatrix);
            this.performChecks = performChecks;
        }
        
        public SymmetricPositiveDefiniteMatrix(DenseDoubleMatrix2D mat, boolean performChecks) throws AlgorithmException {
            dimension = mat.rows();
            theMatrix = new double[dimension][dimension];
            checkSizeSymmetry(mat);
            theMatrix = mat.toArray();
            decompObject = new CholeskyDecomposition(dimension, performChecks);
            decompObject.factor(theMatrix);
            this.performChecks = performChecks;
        }
        
        public SymmetricPositiveDefiniteMatrix(SymmetricPositiveDefiniteMatrix theMat) throws AlgorithmException {
            dimension = theMat.getDimension();
            theMatrix = theMat.getMatrix();
            decompObject = new CholeskyDecomposition(dimension, performChecks);
            decompObject.setFromDecomposition(theMat.decompObject);
            performChecks = theMat.performChecks;
        }
        
        private void checkSizeSymmetry(double[][] mat) throws AlgorithmException {
            if (mat.length != dimension) throw new AlgorithmException("No. rows in the matrix is not correct.");
            for (int i = 0; i < dimension; i++) {
                if (mat[i].length != dimension) throw new AlgorithmException("The matrix is not square.");
                for(int j=0; j<dimension; j++) {
                    if (Math.abs(mat[i][j] - mat[j][i]) > 1e-5) throw new AlgorithmException("The matrix is not symmetric.");
                }
            }
        }
        
        private void checkSizeSymmetry(DenseDoubleMatrix2D mat) throws AlgorithmException {
            if (mat.rows() != dimension) throw new AlgorithmException("No. rows in the matrix is not correct.");
            if (mat.columns() != dimension) throw new AlgorithmException("The matrix is not square.");
            for (int i = 0; i < dimension; i++) {
                for(int j=0; j<dimension; j++) {
                    if (Math.abs(mat.getQuick(i,j) - mat.getQuick(j,i)) > 1e-5) throw new AlgorithmException("The matrix is not symmetric.");
                }
            }
        }
        
        public void setMatrix(double[][] mat) throws AlgorithmException {
            if(performChecks) checkSizeSymmetry(mat);
            for (int i = 0; i < dimension; i++) {
                System.arraycopy(mat[i], 0, theMatrix[i], 0, dimension);
            }
            decompObject.factor(theMatrix);
        }
        
        /* This is essentially a copy method */
        public void setMatrix(SymmetricPositiveDefiniteMatrix mat) throws AlgorithmException {
            if(performChecks) if (mat.getDimension() != dimension) throw new AlgorithmException("No. rows in the matrix is not correct.");
            for (int i = 0; i < dimension; i++) {
                System.arraycopy(mat.theMatrix[i], 0, theMatrix[i], 0, dimension);
            }
            decompObject.setFromDecomposition(mat.decompObject);
        }
        
        public int getDimension() {
            return dimension;
        }
        
        public double getEntry(int i, int j) {
            return theMatrix[i][j];
        }
        
        public double[][] getMatrix() {
            double[][] mat = new double[dimension][dimension];
            for (int i = 0; i < dimension; i++) {
                System.arraycopy(theMatrix[i], 0, mat[i], 0, dimension);
            }
            return mat;
        }
        
        public void getMatrix(double[][] mat) {
            for (int i = 0; i < dimension; i++) {
                System.arraycopy(theMatrix[i], 0, mat[i], 0, dimension);
            }
        }
        
        public void getMatrix(DenseDoubleMatrix2D mat) {
            mat.assign(theMatrix);
        }
        
        public double[][] getInverseMatrix() throws AlgorithmException {
            double[][] mat = new double[dimension][dimension];
            decompObject.inverse(mat);
            return mat;
        }
        
        public double getLogDeterminant() throws AlgorithmException {
            return decompObject.logDeterminant();
        }
        
        /**
         * Computes the transpose of a square matrix in place.
         */
        public static void transposeSqMatrix(double[][] matrix) throws AlgorithmException {
            int m = matrix.length;
            if (matrix[0].length != m) throw new AlgorithmException("Matrix is not square.");
            for (int i = 0; i < m; i++) {
                for (int j = i + 1; j < m; j++) {
                    double temp = matrix[i][j];
                    matrix[i][j] = matrix[j][i];
                    matrix[j][i] = temp;
                }
            }
        }
        
        /** 
         * Returns lower triangle of Cholesky decomposition in DenseDoubleMatrix2D object 
         */
        public void getLowerTriangle(DenseDoubleMatrix2D mat) {
            mat.assign(decompObject.lowerTriangle);
        }
        
        /**
         * Wrappers for Cholesky functions
         */
        public void solve(double[] b, double[] x) throws AlgorithmException {
            decompObject.solve(b, x);
        }
        
        public void solve(double[][] b, double[][] x) throws AlgorithmException {
            decompObject.solve(b, x);
        }

        public void lowerTriangleMult(double[] y, double[] b) throws AlgorithmException {
            decompObject.lowerTriangleMult(y, b);
        }

        public void upperTriangleMult(double[] y, double[] b) throws AlgorithmException {
            decompObject.upperTriangleMult(y, b);
        }

        public void lowerTriangleSolve(double[] b, double[] y) throws AlgorithmException {
            decompObject.lowerTriangleSolve(b, y);
        }

        public void upperTriangleSolve(double[] b, double[] y) throws AlgorithmException {
            decompObject.upperTriangleSolve(b, y);
        }

        
        private class CholeskyDecomposition implements java.io.Serializable {

            private final int dimension;
            private double[][] lowerTriangle;
            private boolean factorisationCalled = false;
            private double[] workingVector;
            private boolean performChecks;

            public CholeskyDecomposition(int n, boolean performChecks) {
                dimension = n;
                lowerTriangle = new double[n][n];
                workingVector = new double[n];
                this.performChecks = performChecks;
            }
            
            private void setFromDecomposition(CholeskyDecomposition decomp) throws AlgorithmException {
                if(performChecks && decomp.dimension != dimension) throw new AlgorithmException("Matrices do not have matching dimensions.");
                factorisationCalled = decomp.factorisationCalled;
                for(int i=0; i<dimension; i++) System.arraycopy(decomp.lowerTriangle[i], 0, lowerTriangle[i], 0, dimension);
            }

            /**
             * Computes the Cholesky decomposition of the matrix: note method
             * assumes the theMatrix is square
             */
            public void factor(double[][] theMatrix) throws AlgorithmException {
                int i, j, k;
                double sum;
                if (performChecks && theMatrix.length != dimension) {
                    throw new AlgorithmException("Error in"
                            + " computing Cholesky decomposition: matrix does not have the"
                            + " correct number of rows.");
                }
                for (i = 0; i < dimension; i++) {
                    for (j = i; j < dimension; j++) {
                        sum = theMatrix[i][j];
                        for (k = i - 1; k >= 0; k--) {
                            sum -= lowerTriangle[i][k] * lowerTriangle[j][k];
                        }
                        if (i == j) {
                            if (sum <= 0.0) {
                                printMatrix(theMatrix);
                                throw new NotPositiveDefiniteException ("Error in computing"
                                        + " Cholesky decomposition: matrix is not positive"
                                        + " definite."); // theMatrix, given rounding errors, is
                            }                                                 // not positive definite
                            lowerTriangle[i][i] = Math.sqrt(sum);
                        } else {
                            lowerTriangle[j][i] = sum / lowerTriangle[i][i];
                        }
                    }
                }
                for (i = 0; i < dimension; i++) {
                    for (j = 0; j < i; j++) {
                        lowerTriangle[j][i] = 0.0;
                    }
                }
                factorisationCalled = true;
            }
            
            /**
             * For the matrix A on which factor() was most recently invoked,
             * solves A x = b, returning the solution in x
             */
            public void solve(double[] b, double[] x) throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                int i, k;
                double sum;
                if (performChecks && (b.length != dimension || x.length != dimension)) {
                    throw new AlgorithmException("Error "
                            + "in solving linear equation. Vectors do not have the correct length.");
                }
                // Solve L y = b where y = workingVector
                for (i = 0; i < dimension; i++) {
                    sum = b[i];
                    for (k = i - 1; k >= 0; k--) {
                        sum -= lowerTriangle[i][k] * workingVector[k];
                    }
                    workingVector[i] = sum / lowerTriangle[i][i];
                }

                // Solve L^T x = y
                for (i = dimension - 1; i >= 0; i--) {
                    sum = workingVector[i];
                    for (k = i + 1; k < dimension; k++) {
                        sum -= lowerTriangle[k][i] * x[k];
                    }
                    x[i] = sum / lowerTriangle[i][i];
                }
            }

            /**
             * For the matrix A on which factor() was most recently invoked,
             * solves A X = B, returning the solution in X
             */
            public void solve(double[][] b, double[][] x) throws AlgorithmException {
                transposeSqMatrix(b);
                for (int i = 0; i < b.length; i++) {
                    solve(b[i], x[i]);
                }
                transposeSqMatrix(x);
            }

            /**
             * For the matrix A = L L^T on which factor() was most recently
             * invoked, computes L y = b returning the answer in b
             */
            public void lowerTriangleMult(double[] y, double[] b) throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                int i, j;
                if (performChecks && (y.length != dimension || b.length != dimension)) {
                    throw new AlgorithmException("Error "
                            + "in performing multiplication. Vectors do not have the correct length.");
                }
                for (i = 0; i < dimension; i++) {
                    b[i] = 0.0;
                    for (j = 0; j <= i; j++) {
                        b[i] += lowerTriangle[i][j] * y[j];
                    }
                }
            }

            /**
             * For the matrix A = L L^T on which factor() was most recently
             * invoked, computes L^T y = b returning the answer in b
             */
            public void upperTriangleMult(double[] y, double[] b) throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                int i, j;
                if (performChecks && (y.length != dimension || b.length != dimension)) {
                    throw new AlgorithmException("Error "
                            + "in performing multiplication. Vectors do not have the correct length.");
                }
                for (i = 0; i < dimension; i++) {
                    b[i] = 0.0;
                    for (j = i; j < dimension; j++) {
                        b[i] += lowerTriangle[j][i] * y[j];
                    }
                }
            }

            /**
             * For the matrix A = L L^T on which factor() was most recently
             * invoked, solves L y = b returning the answer in y
             */
            public void lowerTriangleSolve(double[] b, double[] y) throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                int i, j;
                double sum;
                if (performChecks && (b.length != dimension || y.length != dimension)) {
                    throw new AlgorithmException("Error "
                            + "in solving linear equation. Vectors do not have the correct length.");
                }
                for (i = 0; i < dimension; i++) {
                    sum = b[i];
                    for (j = 0; j < i; j++) {
                        sum -= lowerTriangle[i][j] * y[j];
                    }
                    y[i] = sum / lowerTriangle[i][i];
                }
            }

            /**
             * For the matrix A = L L^T on which factor() was most recently
             * invoked, solves L^T y = b returning the answer in y
             */
            public void upperTriangleSolve(double[] b, double[] y) throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                int i, j;
                double sum;
                if (performChecks && (b.length != dimension || y.length != dimension)) {
                    throw new AlgorithmException("Error "
                            + "in solving linear equation. Vectors do not have the correct length.");
                }
                for (i = dimension - 1; i >= 0; i--) {
                    sum = b[i];
                    for (j = i + 1; j < dimension; j++) {
                        sum -= lowerTriangle[j][i] * y[j];
                    }
                    y[i] = sum / lowerTriangle[i][i];
                }
            }

            /**
             * For the matrix A on which factor() was most recently invoked,
             * inverts A and returns the answer in invMatrix Note method assumes
             * the invMatrix is square
             */
            public void inverse(double[][] invMatrix) throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                int i, j, k;
                double sum;
                if (performChecks && (invMatrix.length != dimension)) {
                    throw new AlgorithmException("Error "
                            + "in inverting matrix. Return matrix does not have the correct "
                            + "number of rows.");
                }
                // Solve L y = b storing y in x
                for (i = 0; i < dimension; i++) {
                    for (j = 0; j <= i; j++) {
                        sum = (i == j ? 1.0 : 0.0);
                        for (k = i - 1; k >= j; k--) {
                            sum -= lowerTriangle[i][k] * invMatrix[j][k];
                        }
                        invMatrix[j][i] = sum / lowerTriangle[i][i];
                    }
                }
                // Solve L^T x = y
                for (i = dimension - 1; i >= 0; i--) {
                    for (j = 0; j <= i; j++) {
                        sum = (i < j ? 0.0 : invMatrix[j][i]);
                        for (k = i + 1; k < dimension; k++) {
                            sum -= lowerTriangle[k][i] * invMatrix[j][k];
                        }
                        invMatrix[j][i] = sum / lowerTriangle[i][i];
                        invMatrix[i][j] = invMatrix[j][i];
                    }
                }
            }

            /**
             * For the matrix A on which factor() was most recently invoked,
             * returns the log determinant
             */
            public double logDeterminant() throws AlgorithmException {
                if (!factorisationCalled) {
                    throw new AlgorithmException("Error: a Cholesky "
                            + "factorisation has not yet been invoked.");
                }
                double sum = 0.0;
                for (int i = 0; i < dimension; i++) {
                    sum += Math.log(lowerTriangle[i][i]);
                }
                return 2.0 * sum;
            }
        }
        
    }
    
    public static class NotPositiveDefiniteException extends AlgorithmException {
        
        public NotPositiveDefiniteException () {
            super();
        }
        
        public NotPositiveDefiniteException (String msg) {
            super(msg);
        }
        
    }
    
    
    
    /* ---------------------------------------------------------------------- */
    /* Test area: 
    public static void main(String[] args) throws AlgorithmException {
        
        int i, j;
        double[] temp;
        double[] x1 = {0.5,-0.4,-2.3,0.0};
        double[] x2 = {2.9,-1.8,1.0,-0.2};
        double[] x3 = {9.4,10.7,9.3,11.1};
        double[] mean = {-1.0,0.5,0.0,3.0};
        double[] mean2 = {-2.8, 6.9,1.3,0.5};
        double[][] variance = {{4.7,-0.3,1.0,1.9},
                               {-0.3,4.1,-0.5,-1.5},
                               {1.0,-0.5,0.9,0.6},
                               {1.9,-1.5,0.6,2.2}};
        double[][] matrix = {{49.4,-4.9,16.3,36.1},
                             {-4.9,13.0,-6.5,-10.7},
                             {16.3,-6.5,8.4,14.3},
                             {36.1,-10.7,14.3,32.0}};
        
        
        boolean isPrec = false;
        MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(mean,variance,isPrec, true);
        
        for(i=0;i<100000;i++){
            temp = mvn.sample();
            for(j=0;j<4;j++){
                System.out.printf("%7.7f ",temp[j]);
            }
            System.out.println("");
        }
        //System.out.printf("%7.7f %7.7f %7.7f%n",mvn.logpdf(x1),mvn.logpdf(x2),mvn.logpdf(x3));
        mvn.setMean(mean2);
        mvn.setSigma(matrix);
        for(i=0;i<100000;i++){
            temp = mvn.sample();
            for(j=0;j<4;j++){
                System.out.printf("%7.7f ",temp[j]);
            }
            System.out.println("");
        }
        //System.out.printf("%7.7f %7.7f %7.7f%n",mvn.logpdf(x1),mvn.logpdf(x2),mvn.logpdf(x3));
        mvn.setPrecision(matrix);
        for(i=0;i<100000;i++){
            temp = mvn.sample();
            for(j=0;j<4;j++){
                System.out.printf("%7.7f ",temp[j]);
            }
            System.out.println("");
        }
        //System.out.printf("%7.7f %7.7f %7.7f%n",mvn.logpdf(x1),mvn.logpdf(x2),mvn.logpdf(x3));
        
    }
    * 
    * */
    
    public static void printMatrix(double[][] m) {
        for(int i=0; i<m.length; i++) {
            for(int j=0; j<m[i].length; j++) {
                System.out.print(String.format("%7.7f ",m[i][j]));
            }
            System.out.println();
        }
    }
    
    public static void printVector(double[] v) {
        for(int i=0; i<v.length; i++) {
             System.out.print(String.format("%7.7f ",v[i]));
        }
        System.out.println();
    }
    
}