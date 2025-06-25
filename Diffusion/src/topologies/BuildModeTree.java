/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package topologies;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import treebase.AlgorithmError;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;

/*
Class to build a tree from a file of splits and lengths. 
*/
public class BuildModeTree {
        public static void main(String[] args) throws IOException, AlgorithmError {
        String inputFilename = args[0];
        String outputFilename= args[1];
        //String inputFilename = "/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_MCMCOutv11_splitModes.txt";
        //String outputFilename= "/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_MCMCOutv11_splitModes_Tree.txt";
        File inputFile=new File(inputFilename);      
        File outputFile=new File(outputFilename);  
        
       HashMap<Split, Double> splitsToLengths = new HashMap();
       String theSplitString = "";
       String notSplitString = "";
       Double ModalLength = 0.0;
       Split theSplit;
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                i=i+1;
               
            }
            else{
            String[] ss = s.split("]");
            theSplitString = ss[0].replaceAll("\\[", "");
            //System.out.println(theSplitString);
            theSplitString= theSplitString.replaceAll(" ", "");
            notSplitString = ss[1].replaceAll("\\[", "");
            notSplitString= notSplitString.replaceAll(" ", "");
            //System.out.println(notSplitString);
            ModalLength = Double.parseDouble(ss[2]);
           
            theSplit= new Split(theSplitString,notSplitString);
            splitsToLengths.put(theSplit,ModalLength);
            System.out.println(theSplit);
            }

        }
                    Tree theTree = new Tree(splitsToLengths);
                    TreeAsSplits theTreeAsSplits = new TreeAsSplits(theTree);
                    
                    System.out.println(theTreeAsSplits.toString());
                    
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
        
        out.println(theTreeAsSplits.toString());
        
               out.close();  
        }
        
}
