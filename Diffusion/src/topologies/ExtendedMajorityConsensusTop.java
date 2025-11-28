/*
ExtendedMajorityConsensusTop
    Copyright (C) 2025  William M Woodman

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

    Contact the author at:  <w.m.woodman2@ncl.ac.uk>
                           
 */
package topologies;

import MCMC.PosteriorAnalysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import static java.lang.System.in;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import treebase.AlgorithmError;
import treebase.Split;
import treebase.TreeAsSplits;
import treebase.Tree;

/**
 *
 * @author will
 */
public class ExtendedMajorityConsensusTop {

    /**
      Main class for computing and extended Majority Rule Consensus Topology for a data set of trees
     */
    public static void main(String[] args) throws AlgorithmError {

        String inputFilename= "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/OldRokasData/yeast.txt";
        
        //read in the data and convert to TopologicalData
        File inFile = new File(inputFilename);
        TopologicalData theData = null;
        try {
            theData = new TopologicalData(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        
        //calculate the tree with the extended majority consensus top, where all edge lengths are one
        Tree theTree = getExtendedMajorityConsensusTop(theData);
        System.out.println(theTree.toTopologyString());
        
        String outputFileName = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/OldRokasData/yeast.txt";
 
        //save the tree
  File outputFile= new File(outputFileName);
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
             out.println(theTree.toString());
        
             out.close();
             
           
    }
    
    //returns a tree with the extended majority consensus top, where all edge lengths are one
     public static Tree getExtendedMajorityConsensusTop(TopologicalData data) throws AlgorithmError{
         HashSet<Split> theTrivialSplits = data.getTheTrivialSplits();
        
        //get the counts
        HashMap<Split, Integer> theSplitCounts = data.getTheSplitCounts();
        //sort the list of splits by most commonly occuring
        List<Map.Entry<Split, Integer>> list = new LinkedList<Map.Entry<Split, Integer>>(theSplitCounts.entrySet());
        Collections.sort(list, new MyComparator());
        LinkedHashMap<Split, Integer> sortedMap = new LinkedHashMap<Split, Integer>();
        for (Map.Entry<Split, Integer> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        
        Set<Split> theSplits = theSplitCounts.keySet();
        ArrayList<Split> splitsInConsensusTop = new ArrayList();
        HashMap<Split,Double>splitsInConsensusTopHashMap = new HashMap();
        int Compatible;
        int numFound =0;
        int nPrime = data.getNumTaxa()-3;
        
        for(Split key : sortedMap.keySet()){
            Compatible=1;
            for(Split split : splitsInConsensusTop){
                if(!key.isCompatible(split)) Compatible =0;
            }
            if(Compatible==1){
                splitsInConsensusTop.add(key);
                splitsInConsensusTopHashMap.put(key, 1.0);
                numFound = numFound+1;
            }
            if(numFound>=nPrime) break;
        }
        int theTotalCount=data.getTotalCount();
       
        //get the most common compatible splits
        //build the tree
        for(Split trivSplit : theTrivialSplits){
            splitsInConsensusTopHashMap.put(trivSplit, 1.0);
        }
        
        Tree theTree = new Tree(splitsInConsensusTopHashMap);
        return theTree;
     }
    
    
    private static class MyComparator implements Comparator<Map.Entry<Split, Integer>> {
        public int compare(Map.Entry<Split, Integer> o1, Map.Entry<Split, Integer> o2) {
            return -o1.getValue().compareTo(o2.getValue()); // Want big integers first!
        }
    }
    
}
