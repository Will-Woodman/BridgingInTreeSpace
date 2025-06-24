/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
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
