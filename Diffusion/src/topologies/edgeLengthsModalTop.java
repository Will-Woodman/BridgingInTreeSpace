/*
edgeLengthsModalTop
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
import MCMC.PosteriorAnalysis;
import static MCMC.PosteriorAnalysis.outputEdgeLengthsForModalTopology;
import java.io.IOException;
import static MCMC.PosteriorAnalysis.countTopologiesToFile;

/**
 *
 * @author will
 */
public class edgeLengthsModalTop {

    /**
     * main class to get some key information from the output of the MCMC (make sure trees are in the 
     * first column of the file
     */
    public static void main(String[] args) throws IOException {

        String inputFilename = args[0];
        String outputFilename= args[1];

        //edge lengths in the modal topology
        File inputFile=new File(inputFilename);      
        File outputFile=new File(outputFilename);  
        outputEdgeLengthsForModalTopology(outputFile, inputFile);
        
       
    }
    
}
