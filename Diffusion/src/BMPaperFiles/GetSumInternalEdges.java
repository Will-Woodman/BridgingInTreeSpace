/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package BMPaperFiles;

import java.io.File;
import java.io.IOException;
import treebase.TreeAsSplits;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Tree;

/**
 *
 * @author will
 */
public class GetSumInternalEdges {

    /**
     Main class that reads in a tree and prints the sum of the internal edges
     */
    public static void main(String[] args) throws IOException, AlgorithmError, AlgorithmException {
    
        //String inputFilename = "/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_MCMCOutv11_splitModes_Tree.txt";
        String inputFilename = args[0];
        File inputFile=new File(inputFilename);
        
        TreeAsSplits theTree = new TreeAsSplits(new Tree(inputFile));
        System.out.println(theTree.sumLengths(true));
        
    
    }
    
}
