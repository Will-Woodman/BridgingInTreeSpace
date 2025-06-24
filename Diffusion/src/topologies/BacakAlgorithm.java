/*
 * BacakAlgorithm.java

    Copyright (C) 2016  Tom M. W. Nye

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

package topologies;

import geodesics.Geodesic;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.Math.min;
import java.util.ArrayList;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
 * Implementation of the unweighted Bacak algorithm for the Frechet mean 
 */
public class BacakAlgorithm {
    
    public static TreeAsSplits unweightedFM(ArrayList<TreeAsSplits> theTrees, int nits) {
        
        TreeAsSplits currentTree = theTrees.get(0);
        try {
            for (int i=0; i<nits; i++) {
                for (int j=0; j<theTrees.size(); j++) {
                    Geodesic g = new Geodesic(currentTree, theTrees.get(j));
                    double u = 2.0/(i+1);
                    TreeAsSplits nextTree = g.getTree(u/(1.0+u));
                    currentTree = nextTree;
                }
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error constructing geodesic in deterministic Bakac algorithm.");
        }
        return currentTree;

        
    }
    
    public static TreeAsSplits unweightedGM(ArrayList<TreeAsSplits> theTrees, int nits) {
        
        TreeAsSplits currentTree = theTrees.get(0);
        try {
            for (int i=0; i<nits; i++) {
                for (int j=0; j<theTrees.size(); j++) {
                    Geodesic g = new Geodesic(currentTree, theTrees.get(j));
                    double lambda = 1.0/(i+1);
                    double t=min(1, lambda/g.getInternalLength());
                    TreeAsSplits nextTree = g.getTree(1-t);
                    currentTree = nextTree;
                }
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error constructing geodesic in deterministic Bakac algorithm.");
        }
        return currentTree;

        
    }
    
    
    public static void main(String args[]) throws java.io.IOException, AlgorithmError {
        String inputFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/4 Taxon Gene Trees/coal_trees_0.01_0.1_AGT3.txt";
        File inFile = new File(inputFilename);
        BufferedReader br = new BufferedReader(new FileReader(inFile));
        
        ArrayList<TreeAsSplits> theTrees = new ArrayList();
       
        String s;
        while ((s = br.readLine())!=null) {
            Tree t = null;
            TreeAsSplits tSplits = null;
            try {
                t = new Tree(s);
                tSplits= new TreeAsSplits(t);
            } catch (AlgorithmException ex) {
                throw new IOException("Bad tree: "+s);
            }
            theTrees.add(tSplits);
        }
        
        TreeAsSplits theMean = unweightedFM(theTrees,1000);
        System.out.println(theMean);
        System.out.println(theMean.getTree().toTopologyString());
        TreeAsSplits theMedian = unweightedGM(theTrees,1000);
        System.out.println(theMedian);
        
        System.out.println(theMedian.getTree().toTopologyString());
        
        
        
    }
    
    
    
    
    /*public static void main(String[] args) throws java.io.IOException, AlgorithmError {

        if ((args.length<2)||(args.length>3)) {
            System.out.println("geophytterplus.BacakAlgorithm computes Frechet mean and variance. ");
            System.out.println("\nSyntax:-");
            System.out.println("BacakAlgorithm [-i] nits filename");
            System.out.println("where ");
            System.out.println("*-i option: use internal edge lengths only for variance");
            System.out.println("* filename is a file containing Newick strings");
            System.out.println("* nits is number of iterations");
            System.exit(1);
        }

        int maxIts = 0;
        TreeAsSplitsDataSet theData=null;
        
        int c = 0;
        if (args[0].equals("-i")) {
            c=1;
            GeodesicWithLengthOption.setPendantsOff();
        }
        
        try {
            File theFile = new File(args[c+1]);  
            theData = new TreeAsSplitsDataSet(theFile, false, false);
        }
        catch (java.io.IOException ioErr) {
            System.out.println("Fatal error reading tree file. Please check the format. ");
            System.exit(1);
        }

        maxIts = Integer.parseInt(args[c]);


        TreeAsSplits mean = unweightedFM(theData.theTrees, maxIts);
        System.out.println(mean.toString());
        
        double dsq = 0.0;
        for (int i=0; i<theData.numTrees; i++) {
            GeodesicWithLengthOption g = new GeodesicWithLengthOption(mean, theData.getTree(i));
            double l = g.getLength();
            dsq += l*l;
        }

        double var = dsq/theData.numTrees;
        System.out.println(var);
    }
*/
    
}
