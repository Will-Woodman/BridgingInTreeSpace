/*
FrechetMeanMainClass
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
package BMPaperFiles;

import geodesics.Geodesic;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import static topologies.BacakAlgorithm.unweightedFM;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;
import static treebase.TreeAsSplits.rfDistance;
import treedatasets.TreeAsSplitsDataSet;

/**
 Main class to estimate a Frechet sample mean and compute distances from it to the data
 */
public class FrechetMeanMainClass {
   
    public static void main(String[] args) throws IOException, AlgorithmException {
        
        String DataSetFilename = args[0];
        String FrechetParamsFilename = args[1];//file to output Frechet mean
        String outputFilename = args[2];//file to output distances from Frechet mean to the data
        
        File dataSetFile = new File(DataSetFilename);
        TreeAsSplitsDataSet dataSet = new TreeAsSplitsDataSet(dataSetFile);
        
        File FrechetParamsFile = new File(FrechetParamsFilename);
        //get the Frechet parameters
        TreeAsSplits FrechMean = unweightedFM(dataSet.theTrees,1000);
        double FrechVar = 0;
        FrechVar = getFrechetVariance(FrechMean, dataSet.theTrees);
        
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
                 Geodesic g = new Geodesic(dataSet.theTrees.get(i),FrechMean);
                 out.println(i+" "+rfDistance(dataSet.theTrees.get(i),FrechMean)+" "+g.getInternalLength());
             }
        out.close();
    
    //print out the Frechet estimators for the data set
     PrintWriter outt;
             try {
                 FileWriter outt1= new FileWriter(FrechetParamsFile, false);
                 BufferedWriter outt2 =new BufferedWriter(outt1);
                 outt = new PrintWriter(outt2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+FrechetParamsFile.getName()+" failed: "
                         +"writing output to console instead.");
                 outt = new PrintWriter(System.out);
            
}
             
        outt.println("FrechetMean FrechetVariance FrechetMeanTop");

        outt.println(FrechMean.toString() +" "+FrechVar+" "+FrechMean.getTree().toTopologyString());
        
        outt.close();
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
