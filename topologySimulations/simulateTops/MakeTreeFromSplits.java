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
import treebase.AlgorithmException;
import treebase.Tree;

/**
 *
 * @author will
 */
public class MakeTreeFromSplits {

    /**
     Main class that takes a list of splits and their lengths and forms a tree, using a constructor from
     * TreeAsSplits: used to construct 'mode' trees from MCMC output on x0
     */
    public static void main(String[] args) throws IOException, AlgorithmError, AlgorithmException {
    
        String inputFilename = args[0];
        File inputFile=new File(inputFilename);
        String outputFilename=args[1];
        
        TreeAsSplits theTree = new TreeAsSplits(inputFile);
        //print to file
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
             
             
             

            out.println(theTree);
        
        
        
               out.close();  
        
        
    
    }
    
}
