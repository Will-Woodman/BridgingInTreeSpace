/*
GetSumInternalEdges
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
