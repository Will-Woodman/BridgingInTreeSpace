/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */


package topologies;
import MCMC.PosteriorAnalysis;
import java.io.File;
import static MCMC.PosteriorAnalysis.extractTreesFromOutputFile;
import static MCMC.PosteriorAnalysis.outputEdgeLengthsForModalTopology;
import static MCMC.PosteriorAnalysis.countTopologiesAtPoints;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import treebase.TreeAsSplits;
import java.util.LinkedHashMap;
import java.util.Map;
import treebase.AlgorithmError;

/**
 *
 * @author will
 */
public class convertToTopologies {

    /**
     * Main class to convert the outputted list of trees into a list of topologies and save to a new file
     */
    public static void main(String[] args) throws IOException, AlgorithmError {
        // TODO code application logic here
        String inputFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/RokasData/set0_MCMCoutv1.txt";
        File inputFile=new File(inputFilename);
        String outputFilename="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/Paper/RokasData/set0_MCMCoutv1/x0Topsv1.txt";
        File outputFile = new File(outputFilename);
        convertTreesToTopologies(inputFile,outputFile);
    }
    
    
    /**
     * Class to convert the outputted list of trees into a list of topologies and save to a new file
     */
    public static void convertTreesToTopologies(File inFile, File outFile) throws AlgorithmError, IOException{
        //import the trees
         ArrayList<TreeAsSplits> theTrees = extractTreesFromOutputFile(inFile);    
        int noOfTrees= theTrees.size();
        
//print the topologies to the output file
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(outFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
        
    }
             
             out.print("Topology ");
        
        for(TreeAsSplits tree: theTrees){
            out.println(tree.getTree().toTopologyString());
        }
        out.println();
        
               out.close();  
        
    }
       
    
}
