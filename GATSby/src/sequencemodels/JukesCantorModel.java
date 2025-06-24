/**
    JukesCantorModel
    Implementation of the JC69 substitution model

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

import cern.colt.matrix.tdouble.DoubleFactory2D;

/**
 * Implementation of the JC69 substitution model
 */

public class JukesCantorModel implements SubstitutionModel, SequenceDistanceCalculator {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;
    
    private Alphabet theAlphabet = DNAAlphabet.getInstance();
    private double[] stationaryDistrib;
    private double[][] Q;

    /** Constructor: fix rate matrix etc. No parameters! */
    public JukesCantorModel() {

        // Fix stationary distrib
        stationaryDistrib = new double[4];
        stationaryDistrib[0] = 0.25; stationaryDistrib[1] = 0.25;
        stationaryDistrib[2] = 0.25; stationaryDistrib[3] = 0.25;
        // Fix Q
        Q = new double[4][4];
        for (int i=0; i<4; i++) {
            Q[i][i] = -1.0;
            for (int j=0; j<i; j++) {
                Q[i][j] = 0.33333333;
                Q[j][i] = 0.33333333;
            }
        }
        
    }
    
    @Override
    public Alphabet getAlphabet() {
        return theAlphabet;
    }
    
    /* For likelihood calculations ------------------------------------------ */

    /** Get the transition matrix for a certain branch length
     Override default calculation! */
    @Override
    public double[][] getTransitionMatrix(double t) {
        double[][] P = new double[4][4];
        double x = 0.25+0.75*Math.exp(-1.33333333*t);
        double y = 0.25-0.25*Math.exp(-1.33333333*t);
        for (int i=0; i<4; i++) {
            P[i][i] = x;
            for (int j=0; j<i; j++) {
                P[i][j]=y;
                P[j][i]=y;
            }
        }
        return P;
    }
    @Override
    public void getTransitionMatrix(double t, double[][] P) {
        double x = 0.25+0.75*Math.exp(-1.33333333*t);
        double y = 0.25-0.25*Math.exp(-1.33333333*t);
        for (int i=0; i<4; i++) {
            P[i][i] = x;
            for (int j=0; j<i; j++) {
                P[i][j]=y;
                P[j][i]=y;
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
        else return 4.0/3.0;
    }*/
   
    /* Distance calculations ----------------------------------------------- */

    /** Distance between two sequences */
    @Override
    public double calcDistance(int[] a, int[] b, double infiniteDist) throws treebase.AlgorithmError {
        if (a.length!=b.length) throw new treebase.AlgorithmError("Sequence length mismatch in JC69 distance calculation");
        // Count up where a and b match
        double p = 0.0;
        for (int i=0; i<a.length; i++) {
            if ((a[i]>=0)&&(b[i]>=0)&&(a[i]!=b[i])) p+=1.0 ;
            else if ((a[i]<0)||(b[i]<0)) p+=0.75;
        }
        p = p/a.length; // p is proportion mismatch
        if (p>=0.75) return infiniteDist;
        return -0.75*Math.log(1.0-1.3333333*p);
    }

    /** Distance calculation from alignment */
    @Override
    public double[][] calcDistanceMatrix(Alignment a, double infiniteDist) {
        double[][] d = new double[a.getNumTaxa()][a.getNumTaxa()];
        double x;

        try {
            for (int i=0; i<a.getNumTaxa(); i++) {
                d[i][i] = 0.0;
                for (int j=0; j<i; j++) {
                    x = calcDistance(a.getRow(i), a.getRow(j),infiniteDist);
                    d[i][j] = x;
                    d[j][i] = x;
                }
            }
        }
        catch (treebase.AlgorithmError e) {
            javax.swing.JOptionPane.showMessageDialog(null,"Error with alignment -- difference numbers of characters in rows. ","Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }

        return d;
    }
    

    /* --------------------------------------------------------------------- */

 


}
