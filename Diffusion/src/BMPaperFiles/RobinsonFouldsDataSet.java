/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package BMPaperFiles;

import geodesics.Geodesic;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import static treebase.TreeAsSplits.rfDistance;
import treedatasets.TreeAsSplitsDataSet;

/**
 *
 * @author will
 */
public class RobinsonFouldsDataSet {

    /**
     Main class to calculate distances from s given tree to a data set
     */
    public static void main(String[] args) throws IOException, AlgorithmException {
        String DataSetFilename = args[0];
        String sourceTreeFilename = args[1];//file to output Frechet mean
        String outputFilename = args[2];//file to output distances from Frechet mean to the data
        
        File dataSetFile = new File(DataSetFilename);
        TreeAsSplitsDataSet dataSet = new TreeAsSplitsDataSet(dataSetFile);        
        //read in the source tree
        TreeAsSplits sourceTree = new TreeAsSplits (new Tree(new File(sourceTreeFilename)));

        File outputFile = new File(outputFilename);
        //print the Robinson Foulds and geodesic distances from the Frechet mean to the data points
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
             
             int numTrees = dataSet.numTrees;
             out.println("DataPoint" +" "+"RFD"+" "+"GD");
             for(int i=0;i<numTrees;i++){
                 Geodesic g = new Geodesic(dataSet.theTrees.get(i),sourceTree);
                 out.println(i+" "+rfDistance(dataSet.theTrees.get(i),sourceTree)+" "+g.getInternalLength());
             }
        out.close();
    
    }
    
    
    //compute the Frechet variance
    private static Double getFrechetVariance (TreeAsSplits mu, ArrayList<TreeAsSplits> theTrees) throws AlgorithmError{
         int numTrees= theTrees.size();
         double totalSquDist=0;
         for(int i=0;i<numTrees;i++){
         TreePair theTreePair = new TreePair(theTrees.get(i),mu);
         totalSquDist+= theTreePair.squDist;
         //System.out.println(theTreePair.squDist);
        }
         //System.out.println((double)((mu.getNumTaxa()-3)*numTrees));
         return(totalSquDist/((double)((mu.getNumTaxa()-3)*numTrees)));
         
     }
    
    //Auxilliary class for the Frechet variance
    private static class TreePair{
      TreeAsSplits treeA;
      TreeAsSplits treeB;
      double squDist;
      
      private TreePair(TreeAsSplits tA, TreeAsSplits tB) throws AlgorithmError{
          treeA=tA;
          treeB=tB;
          Geodesic g = new Geodesic(tA,tB);
          double dist = g.getInternalLength();
          squDist=dist*dist;

             
      }
     
  }
}
