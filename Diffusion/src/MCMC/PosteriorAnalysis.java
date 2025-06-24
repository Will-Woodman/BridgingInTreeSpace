/*
 * PosteriorAnalysis.java

    Copyright (C) 2018  Tom M. W. Nye

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

package MCMC;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;

/**
 * Analyse posterior samples of trees
 */

public class PosteriorAnalysis {
    
    public static ArrayList<TreeAsSplits> extractTreesFromOutputFile(File theFile) throws IOException {
        return extractTreesFromOutputFile(theFile, 0);
    }
    
    public static ArrayList<TreeAsSplits> extractTreesFromOutputFile(File theFile, int burnIn) throws IOException {
        ArrayList<TreeAsSplits> theTrees = new ArrayList();
        
        // Read in lines and discard. Always do at least one for the header.
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        for (int l=0; l<=burnIn; l++) {
            String firstLine = br.readLine();
        }

        /* Keep reading lines */
        String s=br.readLine();
        while (s!=null) {
            if (!(s.startsWith("#"))) {
                int ind = s.indexOf(";");
                String treeString = s.substring(0, ind+1);
                try {
                    Tree t = new Tree(treeString);
                    TreeAsSplits u = new TreeAsSplits(t);
                    theTrees.add(u);
                } catch (AlgorithmException ex) {
                    System.out.println("Error making a tree in PosteriorAnalysis.");
                }
                
            }
            try {
                s=br.readLine();
            }
            catch (IOException anErr) {
                break;
            }
        }
        return theTrees;
    }
    
    
    private static class MyComparator implements Comparator<Entry<String, Integer>> {
        public int compare(Entry<String, Integer> o1, Entry<String, Integer> o2) {
            return -o1.getValue().compareTo(o2.getValue()); // Want big integers first!
        }
    }
    
    
    private static String getTopologyString(TreeAsSplits t) {
        try {
            Tree u = t.getTree();
            return u.toTopologyString();
        } catch (AlgorithmError ex) {
            System.out.println("AlgorithmError: Unable to get topology string in PosteriorAnalysis.");
        }
        return null;
    }

    
    public static LinkedHashMap<String,Integer> countTopologies(ArrayList<TreeAsSplits> theTrees) {
        HashMap<String,Integer> results = new HashMap();
        
        ArrayList<TreeAsSplits> topologies = new ArrayList();
        ArrayList<Integer> topologyCounts = new ArrayList();
        ArrayList<String> outputTopologies = new ArrayList();
        topologies.add(theTrees.get(0));
        topologyCounts.add(1);
        outputTopologies.add(getTopologyString(theTrees.get(0)));
        
        for (int i=0; i<theTrees.size(); i++) {
            
            int ind = -1;
            for (int j=0; j<topologies.size(); j++) {
                if (theTrees.get(i).matchesTopology(topologies.get(j))) {
                    ind = j;
                }
            }
            if (ind<0) {
                // New topology
                topologies.add(theTrees.get(i));
                outputTopologies.add(getTopologyString(theTrees.get(i)));
                topologyCounts.add(1);
            }
            else {
                // Existing topology
                int countsSoFar = topologyCounts.get(ind).intValue();
                topologyCounts.set(ind,countsSoFar+1);
            }
        }
        
        for (int i=0; i<outputTopologies.size(); i++) {
            results.put(outputTopologies.get(i), topologyCounts.get(i));
        }
        
        /* Now sort the hashmap */
        List<Entry<String, Integer>> list = new LinkedList<Entry<String, Integer>>(results.entrySet());
        Collections.sort(list, new MyComparator());
        LinkedHashMap<String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
        for (Entry<String, Integer> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        return sortedMap;

    }
    
   
    /*
    method to count topologies at various iterations in order to assess convergence of the MCMC
    Each outputted value is the count up to that point, rather than the count between points
    */
    public static LinkedHashMap<String,ArrayList<Integer>> countTopologiesAtPoints(ArrayList<TreeAsSplits> theTrees, int spacing) {
       int noOfPoints = theTrees.size()/spacing;
       HashMap<String,Integer> results = new HashMap();
       LinkedHashMap<String,ArrayList<Integer>> allResults = new LinkedHashMap(); 
       ArrayList<TreeAsSplits> topologies = new ArrayList();
       ArrayList<String> outputTopologies = new ArrayList();
       topologies.add(theTrees.get(0));

        outputTopologies.add(getTopologyString(theTrees.get(0)));
        ArrayList<Integer> topologyCounts = new ArrayList();
        topologyCounts.add(1);
        /*
        Count topologies up to the point k*spacing
        We work backwards so that all the topologies are added in the first run through
        */
        for (int k = noOfPoints; k>0;k--){    
            for (int l=0; l<topologyCounts.size();l++)
            {
                topologyCounts.set(l,0);
            }
            
        for (int i=0; i<k*spacing; i++){  
            //work out which topology the next sampled point has
            int ind = -1;
            for (int j=0; j<topologies.size(); j++) {
                if (theTrees.get(i).matchesTopology(topologies.get(j))) {
                    ind = j;
                }
            }
            if (ind<0) {
                // New topology
                topologies.add(theTrees.get(i));
                outputTopologies.add(getTopologyString(theTrees.get(i)));
                topologyCounts.add(1);
            }
            else {
                // Existing topology
            
                int countsSoFar = topologyCounts.get(ind).intValue();
                topologyCounts.set(ind,countsSoFar+1);
                
            }
        }
        
        for (int i=0; i<outputTopologies.size(); i++) {
            String theTop = outputTopologies.get(i);
            /* double[] vec =  allResults.get(theTop);
            if(vec==null){
                double[] vec2 = new double[noOfPoints];
                vec2[0]=topologyCounts.get(i);
                vec=vec2;
            }
            else{
            vec[k]=topologyCounts.get(i);
            allResults.put(theTop, vec);
            }
            */
            //match the topology index to the name in the HashMap
            if(allResults.get(theTop)==null){
                //if topology not seen before then create a new list of counts
                ArrayList<Integer> theList2 = new ArrayList<Integer>();
                theList2.add(topologyCounts.get(i));
                allResults.put(theTop,theList2);
            }
            else{
                //otherwise add the new count to the list for the corresponding topology 
            Integer theCount = topologyCounts.get(i);
            ArrayList theList =allResults.get(theTop);
            theList.add(theCount);
            allResults.put(theTop, theList);
            }
        }
        }
        /* Now sort the hashmap - not doing that for now at least*/
        return allResults;

    }
    
    public static void countTopologiesToScreen(File inFile) throws IOException {
        ArrayList<TreeAsSplits> theTrees = extractTreesFromOutputFile(inFile);
        countTopologiesToScreen(theTrees);
    }
   
    
    
    public static void countTopologiesToScreen(ArrayList<TreeAsSplits> theTrees) {
        LinkedHashMap<String,Integer> theCounts = countTopologies(theTrees);
        for (Map.Entry<String, Integer> entry : theCounts.entrySet()) {
            System.out.println(entry.getKey() + " " + entry.getValue());
        }
    }
    
    
     public static void countTopologiesToFile( File inFile, File outputFile) throws IOException{
        ArrayList<TreeAsSplits> theTrees = extractTreesFromOutputFile(inFile);
        countTopologiesToFile(theTrees, outputFile);
    }
     
     public static void countTopologiesToFile(ArrayList<TreeAsSplits> theTrees, File outputFile ) throws IOException{
         
    
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(outputFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outputFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
             }
                  LinkedHashMap<String,Integer> theCounts = countTopologies(theTrees);
             
                  
        for (Map.Entry<String, Integer> entry : theCounts.entrySet()) {
            out.println(entry.getKey() + " " + entry.getValue());
        }
        
               out.close();  
    
     }
    /*---- Output edge lengths -----------------------------------------------*/
    
    public static void outputEdgeLengthsForModalTopology(File outFile, File inFile) throws IOException {
        ArrayList<TreeAsSplits> t = extractTreesFromOutputFile(inFile);
        outputEdgeLengthsForModalTopology(outFile, t);
    }
    
    
    public static void outputEdgeLengthsForModalTopology(File outFile, File inFile, int burnIn) throws IOException {
        ArrayList<TreeAsSplits> t = extractTreesFromOutputFile(inFile, burnIn);
        outputEdgeLengthsForModalTopology(outFile, t);
    }
    
    
    public static void outputEdgeLengthsForModalTopology(File outFile, ArrayList<TreeAsSplits> theTrees) throws IOException {
    
        LinkedHashMap<String,Integer> myList = countTopologies(theTrees);
        String modalTop = myList.keySet().iterator().next();
        outputEdgeLengthsForSingleTopology(outFile, theTrees, modalTop);
    }
    
    
    public static void outputEdgeLengthsForSingleTopology(File outFile, File inFile, String topology) throws IOException {
        ArrayList<TreeAsSplits> t = extractTreesFromOutputFile(inFile);
        outputEdgeLengthsForSingleTopology(outFile, t, topology);
    }
    
    
    public static void outputEdgeLengthsForSingleTopology(File outFile, File inFile, int burnIn, String topology) throws IOException {
        ArrayList<TreeAsSplits> t = extractTreesFromOutputFile(inFile, burnIn);
        outputEdgeLengthsForSingleTopology(outFile, t, topology);
    }
    
    
    public static void outputEdgeLengthsForSingleTopology(File outFile, ArrayList<TreeAsSplits> theTrees, String topology) throws IOException {

        TreeAsSplits theTopo = null;
        for (int k=0; k<theTrees.size(); k++) {
            String s = getTopologyString(theTrees.get(k));
            if (s.equals(topology)) {
                theTopo = theTrees.get(k);
            }
        }
        if (theTopo==null) {
            System.out.println("Error finding tree with matching topology in PosteriorAnalysis. No match found.");
        }
                
        HashSet<Split> pendants = theTopo.getSplits();
        HashSet<Split> interiorSplits = theTopo.getNonTrivialSplits();
        pendants.removeAll(interiorSplits);
        
        ArrayList<Split> orderedInteriorSplits = new ArrayList();
        orderedInteriorSplits.addAll(interiorSplits);
        Collections.sort(orderedInteriorSplits);
        ArrayList<Split> orderedPendantSplits = new ArrayList();
        orderedPendantSplits.addAll(pendants);
        Collections.sort(orderedPendantSplits);
        
        FileWriter fileWriter = new FileWriter(outFile);
        fileWriter.write("# ");
	
        // Header
        for (int i=0; i<orderedInteriorSplits.size(); i++) {
            fileWriter.write(orderedInteriorSplits.get(i).toOrderedShortString());
            fileWriter.write(" ");
        }
        for (int i=0; i<orderedPendantSplits.size(); i++) {
            fileWriter.write(orderedPendantSplits.get(i).toOrderedShortString());
            fileWriter.write(" ");
        }
        fileWriter.write("\n");
       
        for (int j=0; j<theTrees.size(); j++) {
            for (int i=0; i<orderedInteriorSplits.size(); i++) {
                Split p = orderedInteriorSplits.get(i);
                fileWriter.write(String.format("%7.7f ", theTrees.get(j).getSplitLength(p)));
            }
            for (int i=0; i<orderedPendantSplits.size(); i++) {
                Split p = orderedPendantSplits.get(i);
                fileWriter.write(String.format("%7.7f ", theTrees.get(j).getSplitLength(p)));
             }
            fileWriter.write("\n");             
        }
        
        fileWriter.flush();
        fileWriter.close();
    }
    
    
    /*---- Ouput all numerical params ----------------------------------------*/
    
    public static void outputNumericalParamsForModalTopology(File outFile, File inFile) throws IOException {
        outputNumericalParamsForModalTopology(outFile, inFile, 0);
    }
    
    
    public static void outputNumericalParamsForModalTopology(File outFile, File inFile, int burnIn) throws IOException {
        outputNumericalParamsForSingleTopology(outFile, inFile, burnIn, "");
    }
    
    
    public static void outputNumericalParamsForSingleTopology(File outFile, File inFile, String topology) throws IOException {
        outputNumericalParamsForSingleTopology(outFile, inFile, 0, topology);
    }
    
       
    public static void outputNumericalParamsForSingleTopology(File outFile, File inFile, int burnIn, String topology) throws IOException {

        ArrayList<TreeAsSplits> theTrees = new ArrayList();
        ArrayList<Double> like = new ArrayList();
        ArrayList<Double> t0 = new ArrayList();
        
        // Read in lines and discard. Always do at least one for the header.
        BufferedReader br = new BufferedReader(new FileReader(inFile));
        for (int l=0; l<=burnIn; l++) {
            String firstLine = br.readLine();
        }

        /* Keep reading lines */
        String s=br.readLine();
        while (s!=null) {
            if (!(s.startsWith("#"))) {
                int ind = s.indexOf(";");
                String treeString = s.substring(0, ind+1);
                try {
                    Tree t = new Tree(treeString);
                    TreeAsSplits u = new TreeAsSplits(t);
                    theTrees.add(u);
                } catch (AlgorithmException ex) {
                    System.out.println("Error making a tree in PosteriorAnalysis.");
                }
                
                // Now parse the rest of the string to extract t0 and likelihood
                int space = s.indexOf(" ", ind+2);
                String t0Str = s.substring(ind+2, space);
                String likeStr = s.substring(space+1);
                t0.add(Double.valueOf(t0Str));
                like.add(Double.valueOf(likeStr));
            }
            try {
                s=br.readLine();
            }
            catch (IOException anErr) {
                break;
            }
        }
            
        TreeAsSplits theTopo = null;
        if (topology=="") {
            // Find the modal topology
            LinkedHashMap<String,Integer> myList = countTopologies(theTrees);
            topology = myList.keySet().iterator().next();
        }

        for (int k=0; k<theTrees.size(); k++) {
            String f = getTopologyString(theTrees.get(k));
            if (f.equals(topology)) {
                theTopo = theTrees.get(k);
            }
        }
        if (theTopo==null) {
            System.out.println("Error finding tree with matching topology in PosteriorAnalysis. No match found.");
        }
                
        HashSet<Split> pendants = theTopo.getSplits();
        HashSet<Split> interiorSplits = theTopo.getNonTrivialSplits();
        pendants.removeAll(interiorSplits);
        
        ArrayList<Split> orderedInteriorSplits = new ArrayList();
        orderedInteriorSplits.addAll(interiorSplits);
        Collections.sort(orderedInteriorSplits);
        ArrayList<Split> orderedPendantSplits = new ArrayList();
        orderedPendantSplits.addAll(pendants);
        Collections.sort(orderedPendantSplits);
        
        FileWriter fileWriter = new FileWriter(outFile);
        fileWriter.write("# ");
	
        // Header
        for (int i=0; i<orderedInteriorSplits.size(); i++) {
            fileWriter.write(orderedInteriorSplits.get(i).toOrderedShortString());
            fileWriter.write("|");
        }
        for (int i=0; i<orderedPendantSplits.size(); i++) {
            fileWriter.write(orderedPendantSplits.get(i).toOrderedShortString());
            fileWriter.write("|");
        }
        // Output t0
        fileWriter.write("t0|");
        // Output likelihood
        fileWriter.write("likelihood|");
        
        fileWriter.write("\n");
       
        for (int j=0; j<theTrees.size(); j++) {
            if (theTrees.get(j).matchesTopology(theTopo)) {
                for (int i=0; i<orderedInteriorSplits.size(); i++) {
                    Split p = orderedInteriorSplits.get(i);
                    fileWriter.write(String.format("%7.7f|", theTrees.get(j).getSplitLength(p)));
                }
                for (int i=0; i<orderedPendantSplits.size(); i++) {
                    Split p = orderedPendantSplits.get(i);
                    fileWriter.write(String.format("%7.7f|", theTrees.get(j).getSplitLength(p)));
                }

                // Output t0
                fileWriter.write(String.format("%7.7f|",t0.get(j).doubleValue()));
                // Output likelihood
                fileWriter.write(String.format("%7.7f",like.get(j).doubleValue()));

                fileWriter.write("\n");      
            }
        }
        
        fileWriter.flush();
        fileWriter.close();
    }
    

    /*------------------------------------------------------------------------*/
    
    public static void main(String[] args) throws IOException {
        File theFile = new File("/home/c1032934/Documents/Netbeans/TopInf20240503/x0t0Inference2024/5Taxa/Data_10pts_disp_0.1.txt");
        //ArrayList<TreeAsSplits> theTrees = extractTreesFromOutputFile(theFile);
        //countTopologiesToScreen(theTrees);
        File outFile = new File("/home/c1032934/Documents/Netbeans/TopInf20240503/x0t0Inference2024/5Taxa/Data_10pts_disp_0.1_Tops.txt");
        countTopologiesToFile(theFile,outFile);
    }
    
}
