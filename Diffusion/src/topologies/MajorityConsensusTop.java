/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package topologies;

import java.io.File;
import static java.lang.System.in;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import treebase.AlgorithmError;
import treebase.Split;
import treebase.TreeAsSplits;
import treebase.Tree;

/**
 *
 * @author will
 */
public class MajorityConsensusTop {

    /**
     Main class for computing the Majority Rule Consensus Topology for a data set of trees
     */
    public static void main(String[] args) throws AlgorithmError {

        String inputFilename= "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/RokasData/1070_yeast_genetreesCleaned.txt";
      
        //read im the data and convert to TopologicalData
        File inFile = new File(inputFilename);
        TopologicalData theData = null;
        try {
            theData = new TopologicalData(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }
        
        //calculate the tree with the majority consensus top, where all edge lengths are one
        Tree theTree = getMajorityConsensusTop(theData);
         
        //print the MRCT
        System.out.println(theTree.toTopologyString());
        
    }
    
    //returns a tree with the majority consensus top, where all edge lengths are one
    public static Tree getMajorityConsensusTop(TopologicalData data) throws AlgorithmError{
        HashSet<Split> theTrivialSplits = data.getTheTrivialSplits();
        
        //get the counts
        HashMap<Split, Integer> theSplitCounts = data.getTheSplitCounts();
        //get total count returns the number of trees in the data set
        int theTotalCount=data.getTotalCount();
       
        //get the splits that appear in more than half of the gene trees
        Set<Split> theSplits = theSplitCounts.keySet();
        HashMap<Split,Double> splitsInConsensusTop = new HashMap();
        int theCount=0;
        
        //add in splits in more than 50% of the trees
        Double theLength=1.0;
        for(Split split :theSplits){
           theCount= 2* theSplitCounts.get(split);
        if( theCount>theTotalCount){
            splitsInConsensusTop.put(split, theLength);
        }
        }
        
        //add in the trivial splits with length 1
        for(Split trivSplit : theTrivialSplits){
            splitsInConsensusTop.put(trivSplit, theLength);
        }
        
        int temp=1;
        //build the tree from the splits
        Tree theTree = new Tree(splitsInConsensusTop);

        
        return theTree;
    }
    
}
