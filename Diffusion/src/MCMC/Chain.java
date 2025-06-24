/*
    Chain
    Copyright (C) 2012  Tom M. W. Nye

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

    Contact the author at:  <tom.nye@ncl.ac.uk>
                            <http://www.mas.ncl.ac.uk/~ntmwn/>
 */

package MCMC;

/**
 * Run a Metropolis Hastings MCMC chain
 */

import cern.jet.random.tdouble.DoubleUniform;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import treebase.AlgorithmException;

public class Chain implements java.io.Serializable, simulation.RandomEngineSeedSetter {
    
     /** Version number for serialization - increment when structural changes are made */
     private static final long serialVersionUID = 1L;

    /** Data:
     Consists of current state, prior, proposal and likelihood functions. */
     protected GlobalState theState, theBackupState;
     protected KernelWrapper theKernel;
     protected double realValuedLogLikelihood = Double.NaN; // Set initially to NaN
     protected DoubleMersenneTwister localCopyRandomEngine; // Provides memory of the thread
     protected boolean initialisedFromSavedChain;           // specific random engine in case
                                                            // of deserialization

     /** Constructor */
     public Chain(GlobalState initialState, GlobalState initialBackupState, KernelWrapper kernel) {
         theState = initialState;
         theBackupState = initialBackupState;
         theKernel = kernel;
         initialisedFromSavedChain = false;
         localCopyRandomEngine = null; // to be set during run method
         try {
               theKernel.checkCompatibility(theState);
               theKernel.checkCompatibility(theBackupState);
         }
         catch (AlgorithmException anErr) {
             // Display error message
             System.out.println("Fatal error: MCMC chain requested for incompatible types. "+anErr.getMessage());
             System.exit(1);
         }
     }
     
     /** Constructor from a serialized chain */
     public Chain(File savedChain) {
         
        Chain chain = null;
        // Deserialize chain
        try {
            java.io.FileInputStream fileIn = new java.io.FileInputStream(savedChain);
            java.io.ObjectInputStream in = new java.io.ObjectInputStream(fileIn);
            chain = (Chain) in.readObject();
            in.close();
            fileIn.close();
        } catch (java.io.IOException anEx) {
            anEx.printStackTrace();
        } catch (ClassNotFoundException anEx) {
            System.out.println("Chain class not found.");
            anEx.printStackTrace();
        }
        theState = chain.theState;
        theBackupState = chain.theBackupState;
        theKernel = chain.theKernel;
        realValuedLogLikelihood = chain.realValuedLogLikelihood;
        initialisedFromSavedChain = true;
        localCopyRandomEngine = chain.localCopyRandomEngine;
     }
     
     public void resetRandomEngineSeed() {
         theKernel.resetRandomEngineSeed();
         localCopyRandomEngine = simulation.Random.getEngine();
     }
     
     /** Run the chain
      Make outFile null for output to the console.
      Choose append = true if you want to append MCMC output to end of outFile.
      Note: passing in proposalsPerSweep ensures that if there are k states being sampled per complete MCMC draw, 
      then thin, nits and burnits are multiples of k. */
     
     public void run(int nits, int burnits, int thin, File outFile) {
         run(nits, burnits, thin, outFile, null, false, true);
     }
     
     public void run(int nits, int burnits, int thin, File outFile, File serFile, boolean append, boolean printAccRate) {
         if(initialisedFromSavedChain) simulation.Random.setEngine(localCopyRandomEngine);
         else localCopyRandomEngine = simulation.Random.getEngine();
         
         int proposalsPerSweep = theKernel.getNumMovesPerSweep();
         nits = nits*proposalsPerSweep;
         burnits = burnits*proposalsPerSweep;
         thin = thin*proposalsPerSweep;
         
         // Prepare output stream
         PrintWriter out;
         if (outFile==null) {
             out = new PrintWriter(System.out);
         }
         else {
             try {
                 out = new PrintWriter(new BufferedWriter(new FileWriter(outFile, append)));
             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
             }
         }
         //Print header:
         if(!append) printHeader(out);
         
         /* Main loop */
         for (int i=0; i<(nits+burnits); i++) {

             try {
                 generateSample();
             }
             catch (treebase.AlgorithmException anEx) {
                 System.out.println("Warning executing MCMC main loop: "+anEx.getMessage());
                 anEx.printStackTrace();
             }
             
             if(!Double.isNaN(theState.getLogLikelihood())) realValuedLogLikelihood = theState.getLogLikelihood();

             // Do any output here
             if (i>=burnits & ((i-burnits+1)%thin)==0) {
                 output(out);
             }
             
         } // End main loop

         // Print acceptance rates / numerical summaries as last lines:
         if(printAccRate) out.print(theKernel.printFinalSummary());

         // Flush output
         if (outFile==null) out.flush();
         else out.close();
         
         // Serialize Chain if required
         if (serFile != null) {
             try {
                 // Create byte-reader
                 java.io.FileOutputStream fileOut = new java.io.FileOutputStream(serFile);
                 // Wrap byte-reader to translate to object
                 java.io.ObjectOutputStream outSt = new java.io.ObjectOutputStream(fileOut);
                 outSt.writeObject(this);
                 outSt.close();
                 fileOut.close();
             } catch (java.io.IOException anEx) {
                 anEx.printStackTrace();
             }
         }

     }
     
     /** Print header */
     public void printHeader(PrintWriter out) {
         out.println(theState.getHeader().trim()+" "+"LogLikelihood");
     }
     
     /** Output current state */
     public void output(PrintWriter out) {
         out.println(theState.getValueString().trim()+" "+String.format("%7.7f",realValuedLogLikelihood));
     }
     
     /** Sample from kernel */
     public void generateSample() throws AlgorithmException{
         theKernel.sample(theState, theBackupState);
     }
     
     /** For the quasistatic stepping stone sampling - want to pass on the final bridges to the next MCMC run
     */
     public GlobalState getState(){
         return(theState);
     }
     
     public GlobalState getBackupState(){
         return(theBackupState);
     }
}
