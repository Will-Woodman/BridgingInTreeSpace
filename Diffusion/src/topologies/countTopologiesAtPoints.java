/*
countTopologiesAtPoints
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


package topologies;
import java.io.File;
import static MCMC.PosteriorAnalysis.extractTreesFromOutputFile;
import static MCMC.PosteriorAnalysis.countTopologiesAtPoints;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import treebase.TreeAsSplits;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author will
 */
public class countTopologiesAtPoints {

    /**
     * Main class to print the output from countTopologiesAtPoints. Use for the posterior sample on
     * x0 to assess convergence.
     */
    public static void main(String[] args) throws IOException {

        
        String inputFilename = args[0];
        String outputFilename=args[1];

        printTopologiesAtPoints(inputFilename,outputFilename);
        
        
    }
    
   
    public static void printTopologiesAtPoints(String inputFilename, String outputFilename) throws IOException{
        
        File inputFile=new File(inputFilename);

        ArrayList<TreeAsSplits> theTrees = extractTreesFromOutputFile(inputFile);
        int spacing =100;        
        LinkedHashMap<String,ArrayList<Integer>> output = countTopologiesAtPoints(theTrees, spacing);
        int noOfPoints = theTrees.size()/spacing;
        System.out.println(noOfPoints);
        System.out.println(theTrees.size());
        
        //print to screen:
        System.out.print("Topology ");
        
        for(int i=0; i<noOfPoints;i++){
            System.out.print((noOfPoints-i)*spacing+" ");
        }
        System.out.println();
            
        for (Map.Entry<String, ArrayList<Integer>> entry : output.entrySet()) {
            System.out.print(entry.getKey());
            for(int i=0; i<entry.getValue().size();i++){
            System.out.print(" "+entry.getValue().get(i));
            
                    }
            System.out.println();
        }
        
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
             
             out.print("Topology ");
        
        for(int i=0; i<noOfPoints;i++){
            out.print((noOfPoints-i)*spacing+" ");
        }
        out.println();
            
        for (Map.Entry<String, ArrayList<Integer>> entry : output.entrySet()) {
            out.print(entry.getKey());
            for(int i=0; i<entry.getValue().size();i++){
            out.print(" "+entry.getValue().get(i));
            
                    }
            out.println();
        }
        
               out.close();  
    }
       
    
}
