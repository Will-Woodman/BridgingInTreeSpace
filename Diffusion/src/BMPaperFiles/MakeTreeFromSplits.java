/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package BMPaperFiles;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import treebase.TreeAsSplits;
import treebase.AlgorithmError;
import treebase.AlgorithmException;

/**
 *
 * @author will
 */
public class MakeTreeFromSplits {
/*
MakeTreeFromSplits
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
