/**
    StaticLikelihoodCalculator
    Class containing static methods for likelihood calculation

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
 * Class containing static methods for likelihood calculation
 */

import treebase.*;
import java.util.*;

public class StaticLikelihoodCalculator {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    /* Static methods for use by classes which implement SubstitutionModel ----- */
    /** Normalize (and return normalisation constant) */
    public static double normalizeQ(double[][] Q, double[] stationaryDistrib) {
        double r = getOverallRate(Q, stationaryDistrib);
        double x;
        for (int i=0; i<stationaryDistrib.length; i++) {
            for (int j=0; j<stationaryDistrib.length; j++) {
                x = Q[i][j];
                Q[i][j] = x/r;
            }
        }
        return r;
    }
    /** Get overall rate */
    public static double getOverallRate(double[][] Q, double[] stationaryDistrib) {
        double r = 0.0;
        for (int i=0; i<stationaryDistrib.length; i++) {
            r += -stationaryDistrib[i]*Q[i][i];
        }
        return r;
    }
    
    /* Likelihood calcs ----------------------------------------------------- */

    protected static Graph.Vertex findBaseVertex(Tree t) {
        Iterator<Graph.Vertex> it = t.getVertexIterator();
        Graph.Vertex v=null;
        for (int i=0; i<t.numVertices(); i++) {
            Graph.Vertex w = it.next();
            if (w.degree()>1) {
                v = w;
                break;
            }
        }
        return v;
    }

    /** Compute likelihood without specifying which vertex to use as root. */
    public static double logLikelihood(Tree theTree, Alignment theAlignment, SubstitutionModel substModel, boolean saveCPU) throws AlgorithmError {
        Graph.Vertex v = findBaseVertex(theTree);
        return logLikelihood(theTree, theAlignment, substModel, v, saveCPU);
    }
    /** Compute likelihood specifying a particular vertex to use as root. */
    public static double logLikelihood(Tree theTree, Alignment theAlignment, SubstitutionModel substModel, Graph.Vertex v, boolean saveCPU) throws AlgorithmError {
        if (theAlignment.getAlphabet()!=substModel.getAlphabet()) throw new AlgorithmError("Wrong alphabet in sequence model likelihood calculation. ");

        int numStates = theAlignment.getAlphabet().size();
        double[][] transMatrix = new double[numStates][numStates];
        
        if (saveCPU) {
            // Compute likelihood for all sites at once
            double[][] x;
            x = recursiveCalcLikelihood(v, null, theAlignment, substModel, 1.0, transMatrix);
            double l = 0.0, s;
            for (int site=0; site<theAlignment.getNumSitePatterns(); site++) {
                s = 0.0;
                for (int j=0; j<numStates; j++) {
                    s += x[site][j]*substModel.stationaryProb(j); // Take product of stationary distrib with likelihood given that letter at root
                }
                l += theAlignment.getNumSitesPerPattern(site)*Math.log(s);
            }
            return l;
        }
        else {
            // Compute likelihood one site at a time
            double l = 0.0;
            double[] x;
            for (int site=0; site<theAlignment.getNumSitePatterns(); site++) {
                x = recursiveCalcLikelihood(v, null, theAlignment, substModel, site, 1.0, transMatrix);
                double s = 0.0;
                for (int j=0; j<numStates; j++) {
                    s += x[j]*substModel.stationaryProb(j); // Take product of stationary distrib with likelihood given that letter at root
                }
                l += theAlignment.getNumSitesPerPattern(site)*Math.log(s);
            }
            return l;
        }

    }

    /* Rate mixture model likelihood */

    /** Compute likelihood for rate mixture model */
    public static double logLikelihood(Tree theTree, Alignment theAlignment, SubstitutionModel substModel, double[] relativeRates, double[] rateProbabilities, boolean saveCPU) throws AlgorithmError {
        Graph.Vertex v = findBaseVertex(theTree);
        return logLikelihood(theTree, theAlignment, substModel, relativeRates, rateProbabilities, v, saveCPU);
    }
    /** Compute likelihood specifying a particular vertex to use as root. */
    public static double logLikelihood(Tree theTree, Alignment theAlignment, SubstitutionModel substModel, double[] relativeRates, double[] rateProbabilities, Graph.Vertex v, boolean saveCPU) throws AlgorithmError {
        if (theAlignment.getAlphabet()!=substModel.getAlphabet()) throw new AlgorithmError("Wrong alphabet in sequence model likelihood calculation. ");

        int numStates = theAlignment.getAlphabet().size();
        double[][] transMatrix = new double[numStates][numStates];
        
        if(saveCPU) {
            // Compute likelihood for all sites at once
            double[][] x=new double[theAlignment.getNumSitePatterns()][numStates];
            for (int k=0; k<relativeRates.length; k++) {
                double[][] y = recursiveCalcLikelihood(v, null, theAlignment, substModel, relativeRates[k], transMatrix);
                arrayAdd(y,x,rateProbabilities[k]);
            }
            
            double l = 0.0, s;
            for (int site=0; site<theAlignment.getNumSitePatterns(); site++) {
                s = 0.0;
                for (int j=0; j<numStates; j++) {
                    s += x[site][j]*substModel.stationaryProb(j); // Take product of stationary distrib with likelihood given that letter at root
                }
                l += theAlignment.getNumSitesPerPattern(site)*Math.log(s);
            }
            return l;
        }
        else {
            // Compute likelihood one site at a time
            double l = 0.0;
            double[] x;

            for (int site=0; site<theAlignment.getNumSitePatterns(); site++) {
                x = new double[numStates];
                for (int k=0; k<relativeRates.length; k++) {
                    double[] y = recursiveCalcLikelihood(v, null, theAlignment, substModel, site, relativeRates[k], transMatrix);
                    arrayAdd(y,x,rateProbabilities[k]);
                }

                double s = 0.0;
                for (int j=0; j<numStates; j++) {
                    s += x[j]*substModel.stationaryProb(j); // Take product of stationary distrib with likelihood given that letter at root
                }
                l += theAlignment.getNumSitesPerPattern(site)*Math.log(s);
            } // End loop thru sites

            return l;
        }

    
    }
    
    
    /* Likelihood utils ----------------------------------------------------- */

    /** Calculate vector of likelihoods at single vertex given a single site in the alignment. */
    public static double[] recursiveCalcLikelihood(Graph.Vertex v, Graph.Vertex fromV, Alignment theAlignment, SubstitutionModel substModel, int siteNum, double relRate, double[][] transMatrix) throws AlgorithmError {
        int n = theAlignment.getAlphabet().size();
        double[] x = new double[n];
        if (v.degree()==1) {
            // At a leaf
            int letter = theAlignment.getLetter(v.label, siteNum);
            if (letter>=0) {
                x[letter] = 1.0;
                return x;
            }
            else {
                Arrays.fill(x,1.0);
                return x;
            }
        }

        // Set all values to 1.0
        for (int j=0; j<n; j++) x[j] = 1.0;

        // Loop thru' children
        HashSet children = v.getNeighbours();
        double s;
        if (fromV!=null) children.remove(fromV);
        Iterator it = children.iterator();
        for (int i=0; i<children.size(); i++) {
            Graph.Vertex w = (Graph.Vertex) it.next();
            // Get likelihood vector for child
            double[] y = recursiveCalcLikelihood(w, v, theAlignment, substModel, siteNum, relRate, transMatrix);
            // Get the transition mtx
            substModel.getTransitionMatrix(w.getNeighbourDistance(v)*relRate,transMatrix);
            // Now multiply y by transMatrix
            for (int j=0; j<n; j++) {
                s = 0.0;
                for (int k=0; k<n; k++) {
                    s += transMatrix[j][k]*y[k];
                }
                x[j] = x[j]*s;
            }
        }
        
        return x;

    }

    /** Calculate vector of likelihoods at single vertex for every site in the alignment. */
    public static double[][] recursiveCalcLikelihood(Graph.Vertex v, Graph.Vertex fromV, Alignment theAlignment, SubstitutionModel substModel, double relRate, double[][] transMatrix) throws AlgorithmError {
        int n = theAlignment.getAlphabet().size();
        if (v.degree()==1) {
            double[][] x = new double[theAlignment.getNumSitePatterns()][n];
            int[] seq = theAlignment.getRow(v.label);
            for (int i=0; i<theAlignment.getNumSitePatterns(); i++) {
                if (seq[i]>=0) {
                    x[i][seq[i]] = 1.0;
                }
                else {
                    Arrays.fill(x[i],1.0);
                }
            }
            return x;
        }

        double[][] x = new double[theAlignment.getNumSitePatterns()][n];
        double s;

        // Loop thru' children
        HashSet children = v.getNeighbours();
        if (fromV!=null) children.remove(fromV);
        Iterator it = children.iterator();
        for (int i=0; i<children.size(); i++) {
            Graph.Vertex w = (Graph.Vertex) it.next();
            // Get likelihood vectors for child
            double[][] y = recursiveCalcLikelihood(w, v, theAlignment, substModel, relRate, transMatrix);
             // Get the transition mtx
            substModel.getTransitionMatrix(w.getNeighbourDistance(v)*relRate,transMatrix);
            // For each site multiply the likelihood vector by the transition matrix
            for (int site=0; site<theAlignment.getNumSitePatterns(); site++) {
                double[] likeVec = y[site]; // Get the likelihood vector
                for (int j=0; j<n; j++) {
                    s = 0.0;
                    for (int k=0; k<n; k++) {
                        s += transMatrix[j][k]*likeVec[k];
                    }
                    if (i==0) x[site][j] = s;
                    else x[site][j] = s*x[site][j];
                }
            }

        }
        return x;

    }

    
    /* Utilities ------------------------------------------------------------ */

    /** Add one array to another */
    public static void arrayAdd(double[] y, double[] x, double w) {
        for (int i=0; i<y.length; i++) {
            x[i] += w*y[i];
        }
    }

    /** Add one array to another */
    public static void arrayAdd(double[][] y, double[][] x, double w) {
        for (int i=0; i<y.length; i++) {
            for (int j=0; j<y[0].length; j++) {
                x[i][j] += w*y[i][j];
            }
        }
    }
     
   
}