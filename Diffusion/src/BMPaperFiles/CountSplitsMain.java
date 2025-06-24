/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package BMPaperFiles;

import bridge.InferBrownianParamsMCMC;
import geodesics.Geodesic;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import static topologies.BacakAlgorithm.unweightedFM;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;
import static treebase.TreeAsSplits.rfDistance;
import treedatasets.TreeAsSplitsDataSet;

/**
 *
 * @author will
 */
public class CountSplitsMain {

    /**
     A method to check the number of count split in a data set (used to check how many times the long split
     * was included in the Rokas data)
     */
    public static void main(String[] args) throws IOException, AlgorithmException {

        String DataSetFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/RokasData/1070_yeast_genetreesCleaned.txt";
        File dataSetFile = new File(DataSetFilename);
        TreeAsSplitsDataSet dataSet = new TreeAsSplitsDataSet(dataSetFile);
        
        // get the counts of all the splits
        HashMap<Split,Integer> splitCounts = countSplits(dataSet);
        
        //print out the counts to the console
        for(Split split:splitCounts.keySet()){
            System.out.println(split.toShortString()+" "+splitCounts.get(split));
        }

    }
   
    
    
     public static HashMap<Split,Integer> countSplits(TreeAsSplitsDataSet data) {
        int numTrees = data.numTrees;
        ArrayList<TreeAsSplits> theTrees = data.theTrees;
        ArrayList<Split> theSplits1 = new ArrayList();
        HashMap<Split,Integer> theSplits = new HashMap();
        
        //get the list of all the splits from all the trees (with repeats)
        for(TreeAsSplits Tree:theTrees){
        for (Split split:Tree.getNonTrivialSplits()) {
            theSplits1.add(split);
        }
        }
        
        //get the counts for each distinct split in the list (without repeats)
        for(Split split:theSplits1){
            if(!theSplits.containsKey(split)){
                int count =0;
                for(Split split2:theSplits1){
                    if(split2.equals(split)) count=count+1;
                    theSplits.put(split,count);
            }
            }
                
            }
        
        
        return theSplits;
    }
     
  }
    
    
    
    

    
