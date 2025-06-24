/*
 * TopologicalData.java

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

package topologies;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treebase.Split;
/**
 * Class representing counts of topologies.
 * Should have just one of these for any data set / MCMC run. 
 * Don't use this class to represented simulated collections of particles. 
 */

public class TopologicalData {
    
    /* Instance data */
    protected HashMap<String,Integer> counts;  // Hashmap from topology string to counts
    protected int totalCount;
    protected int numTaxa=0;
    protected HashMap<Split,Integer> splitCounts;
    protected HashSet<Split> trivialSplits;
    
    
    /* Utility methods ----------------------------------------------------- */
    
    //number of trees in the data set
    public int getTotalCount() {
        return totalCount;
    }
    
    public int getCount(String topString) {
        if (counts.containsKey(topString)) return counts.get(topString);
        else return 0;
    }
    
    public int getCount(Tree theTree) {
        String s = theTree.toTopologyString();
        return counts.get(s);
    }
    
    public int getNumTaxa() {
        return numTaxa;
    }
    
    public HashMap<String,Integer> getTheCounts(){
        return counts;
    }
    
    public HashMap<Split, Integer> getTheSplitCounts(){
        return splitCounts;
    }
    
    public HashSet<Split> getTheTrivialSplits(){
        return trivialSplits;
    }
            
    
    /* Constructors ----------------------------------------------------- */
    
    /**  */
    public TopologicalData(ArrayList<Tree> theTrees) {
        buildFromTreeData(theTrees);
        buildSplitsFromTreeData(theTrees);
    }
    
    
    /** Read in a list of Newick strings in a file*/
    public TopologicalData(File theFile) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(theFile));
        
        ArrayList<Tree> theTrees = new ArrayList();
        String s;
        while ((s = br.readLine())!=null) {
            Tree t = null;
            try {
                t = new Tree(s);
            } catch (AlgorithmException ex) {
                throw new IOException("Bad tree: "+s);
            }
            theTrees.add(t);
        }

        buildFromTreeData(theTrees);
        buildSplitsFromTreeData(theTrees);
        buildTrivialSplits(theTrees);
        
    }
    
    
    public void buildFromTreeData(ArrayList<Tree> theTrees) {
        totalCount = theTrees.size();
        counts = new HashMap();
        counts.put(theTrees.get(0).toTopologyString(), 1);
        
        for (int i=1; i<theTrees.size(); i++) {
            
            int newCount = -1;
            Tree t = theTrees.get(i);
            String s = t.toTopologyString();
            if (numTaxa==0) {
                numTaxa = t.numTaxa();
            }
            
            for (Map.Entry<String, Integer> entry : counts.entrySet()) {
                if (s.equals(entry.getKey())) {
                    // This topology has been seen before
                    newCount = entry.getValue().intValue()+1;
                    break;
                }
            }
        if (newCount<0) {
                // new topology
                counts.put(s, 1);
            }
            else {
                // Existing topology
                counts.remove(s);
                counts.put(s, newCount);
            }
        }
    }
    
      
       
            
       public void buildSplitsFromTreeData(ArrayList<Tree> theTrees) {
        // turn the array list of trees into an array list of trees as splits
        ArrayList<TreeAsSplits> theTreesAsSplits = new ArrayList();
        HashMap<Split,Integer> splitCountsNow = new HashMap();
        for(int i =0; i <theTrees.size(); i++){
            TreeAsSplits theTreeAsSplits = new TreeAsSplits(theTrees.get(i));
            theTreesAsSplits.add(theTreeAsSplits); 
        }
        //get the first lot of splits
        TreeAsSplits firstTreeAsSplits = theTreesAsSplits.get(0);
        HashSet<Split> firstSplits = firstTreeAsSplits.getNonTrivialSplits();
        int t=0;
        for(Split split : firstSplits){
        splitCountsNow.put(split,1);   
        }    
        
        for (int i=1; i<theTreesAsSplits.size(); i++) {
            TreeAsSplits currentTree = theTreesAsSplits.get(i);
            HashSet<Split> currentSplits = currentTree.getNonTrivialSplits();
            for(Split split : currentSplits){   
            int newCount = -1;
            /*Tree t = theTrees.get(i);
            String s = t.toTopologyString();
            if (numTaxa==0) {
                numTaxa = t.numTaxa();
            }
            */
            for (Map.Entry<Split, Integer> entry : splitCountsNow.entrySet()) {
                if (split.equals(entry.getKey())) {
                    // This topology has been seen before
                    newCount = entry.getValue().intValue()+1;
                    break;
                }
            }     
            if (newCount<0) {
                // new topology
                splitCountsNow.put(split, 1);
            }
            else {
                // Existing topology
                splitCountsNow.remove(split);
                splitCountsNow.put(split, newCount);
            }
        
            }
        } 
        splitCounts=splitCountsNow;
       }
       
     
       public void buildTrivialSplits(ArrayList<Tree> theTrees) {
        // turn the array list of trees into an array list of trees as splits
            TreeAsSplits theTreeAsSplits = new TreeAsSplits(theTrees.get(1)); 
            //HashSet<Split> trivSplits = new HashSet();
        HashSet<Split> theSplits = theTreeAsSplits.getSplits();
        HashSet<Split> trivSplits = new HashSet();
        int t=0;
        for(Split theSplit : theSplits){
        if(theSplit.isTerminal()!=null) {
            trivSplits.add(theSplit);
        }
        }    
        trivialSplits=trivSplits;
       }
  
       public static void main(String[] args) throws AlgorithmError {
        // TODO code application logic here
        String inputFilename= "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/4 Taxon Gene Trees/coal_trees_0.10_0.6_AGT0.txt";
        //ArrayList<TreeAsSplits> theTreesAsSplits = new ArrayList();
         File inFile = new File(inputFilename);
        TopologicalData theData = null;
        try {
            theData = new TopologicalData(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        
        //get the counts
        HashSet<Split> theSplitCounts = theData.getTheTrivialSplits();
        int theTotalCount=theData.getTotalCount();

}
       
       
}


