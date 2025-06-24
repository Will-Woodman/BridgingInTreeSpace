/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package simulateTops;

import static MCMC.PosteriorAnalysis.extractTreesFromOutputFile;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Map;
import treebase.TreeAsSplits;
import MCMC.PosteriorAnalysis;
import treebase.AlgorithmError;
import treebase.Tree;

/**
 *
 * @author will
 */
public class OutputTopologyStrings {

    /**
     Class to convert a file of trees (Newick string form) to a file of topologies for further analysis
     */
    public static void main(String[] args) throws IOException, AlgorithmError {
        
        String inputFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/x0Inference/4Taxa/20220511/Unbal_0.5_0.5_t0_8.txt";
        File inputFile=new File(inputFilename);
        String outputFilename="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/x0Inference/4Taxa/20220511/Unbal_0.5_0.5_t0_8_withTops.txt";
        
        //import the trees from the file
        ArrayList<TreeAsSplits> theTrees = extractTreesFromOutputFile(inputFile);
        
        //print the trees and their topologies to the new file
        File outputFile = new File(outputFilename);
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
             
             out.println("DataPoint Topology");
             
            
             
        for(int i=0;i<theTrees.size();i++){
            Tree theTree = theTrees.get(i).getTree(); 
            out.println(theTrees.get(i).toString() + " " + theTree.toTopologyString());
        }
        
        
               out.close();  
        
        
    
    }
    
}
