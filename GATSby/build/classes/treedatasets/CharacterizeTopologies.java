/*
 * CharacterizeTopologies.java

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

package treedatasets;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import treebase.Split;
import treebase.TreeAsSplits;

/**
 * Take an array of trees and find all the topologies
 */

public class CharacterizeTopologies {
    
    public static int[] identifyDifferentTopologies(ArrayList<TreeAsSplits> theTrees, ArrayList<String> outputTopologies) {
        
        int[] topologyIndex = new int[theTrees.size()];
        ArrayList<TreeAsSplits> topologies = new ArrayList();
        
        topologies.add(theTrees.get(0));
        outputTopologies.add(theTrees.get(0).toString(false));
        topologyIndex[0] = 0;
        
        for (int i=1; i<theTrees.size(); i++) {
            
            int ind = -1;
            TreeAsSplits t = theTrees.get(i);
            for (int j=0; j<topologies.size(); j++) {
                if (t.matchesTopology(topologies.get(j))) {
                    ind = j;
                }
            }
            if (ind<0) {
                // New topology
                topologyIndex[i] = topologies.size();
                topologies.add(t);
                outputTopologies.add(t.toString(false));
            }
            else {
                // Existing topology
                topologyIndex[i] = ind;
            }
        }
        
        return topologyIndex;
    }
    
    
    public static void outputTopologiesInTrees(TreeAsSplitsDataSet theData, File outputFile, boolean minEdgeLen) throws IOException {
        int[] topologyIndex = new int[theData.numTrees];
        ArrayList<TreeAsSplits> topologies = new ArrayList();
        
        topologies.add(theData.getTree(0));
        ArrayList<String> outputTopologies = new ArrayList();
        ArrayList<String> resolved = new ArrayList();
        outputTopologies.add(theData.getTree(0).toString(false));
        if (theData.getTree(0).fullyResolved()) resolved.add("");
        else resolved.add("unresolved");
        topologyIndex[0] = 0;
        
        for (int i=1; i<theData.numTrees; i++) {
            
            int ind = -1;
            TreeAsSplits t = theData.getTree(i);
            for (int j=0; j<topologies.size(); j++) {
                if (t.matchesTopology(topologies.get(j))) {
                    ind = j;
                }
            }
            if (ind<0) {
                // New topology
                topologyIndex[i] = topologies.size();
                topologies.add(t);
                outputTopologies.add(t.toString(false));
                if (t.fullyResolved()) resolved.add("");
                else resolved.add("unresolved");
            }
            else {
                // Existing topology
                topologyIndex[i] = ind;
            }
        }
        
        /* Do output */
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        int[] countTopologies = new int[outputTopologies.size()];
        for (int i=0; i<theData.numTrees; i++) {

            TreeAsSplits t = theData.getTree(i);
            out.print(topologyIndex[i]+1);
            countTopologies[topologyIndex[i]]++;
            if (minEdgeLen) {
                // Find min edgelength
                double l=-1;
                HashSet<Split> internalSplits = t.getNonTrivialSplits();
                Iterator<Split> it = internalSplits.iterator();
                while (it.hasNext()) {
                    Split p = it.next();
                    double x = t.getSplitLength(p);
                    if ((l<0)||(x<l)) {
                        l=x;
                    }
                }
               out.print(String.format(" %7.7f", l));
            }
            out.println();
            
        }
        /* Now output topos */
        System.out.println("");
        for (int i=0; i<outputTopologies.size(); i++) {
            out.println("# topology "+(i+1)+" = "+outputTopologies.get(i)+" "+resolved.get(i));
            System.out.println(outputTopologies.get(i)+" "+countTopologies[i]);
        }
        out.flush();
        out.close();
    }
    
       
    public static void main(String[] args) throws IOException {

//        String dataFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/apicomplexa/apicomplexa_non_out.txt"; 
//        String outputFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/apicomplexa_topologies.txt"; 
        
//        String dataFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/api_run1_proj_trees.txt"; 
//        String outputFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/apicomplexa_projected_topologies.txt"; 
        
//        String dataFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/lungfish/lungfish_filtered.txt"; 
//        String outputFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/lungfish_topologies.txt"; 
        
//        String dataFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/lung_filtered_run1_projected_trees.txt"; 
//        String outputFilename = "/home/ntmwn/research/PCA/GeoPhytter+/results/lungfish_projected_topologies.txt"; 

//        String dataFilename = "/data/pca/data/battistuzzi.txt"; 
//        String outputFilename = "/home/ntmwn/tmp/topologies.txt"; 
        
//        String dataFilename = "/data/pca/data/redelings_cleaned.txt"; 
//        String outputFilename = "/home/ntmwn/tmp/topologies.txt"; 
        
//        String dataFilename = "/data/pca/data/bacteria/keepOne_DTroll.nwk"; 
//        String outputFilename = "/home/ntmwn/tmp/topologies.txt"; 
        
        String dataFilename = "/home/tom/research/Diffusion/bridge/sims/sampled_trees_5_0.05.txt"; 
        String outputFilename = "/home/tom/tmp/topologies.txt"; 
        
        /* Read in the data */
        File inFile = new File(dataFilename);
        TreeAsSplitsDataSet theData = null;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        File outFile = new File(outputFilename);

        outputTopologiesInTrees(theData, outFile, true);
    }
    
}
