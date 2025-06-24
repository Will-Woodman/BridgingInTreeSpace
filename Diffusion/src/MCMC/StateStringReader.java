/*
 * StateStringReader
    Copyright (C) 2013  Sarah E. Heaps

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

    Contact the author at:  <sarah.heaps@ncl.ac.uk>
    
 */

package MCMC;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * Abstract class which reads MCMC output from a file. For each line, the corresponding
 * state is constructed, and this can be post-processed in some way, printing the
 * output to screen or to another file.
 */

public abstract class StateStringReader {
    
    protected GlobalState theState;
    
    /** Constructor: argument is a state of the kind represented in the MCMC output.
        The names of the states in the initial global state must be the same as the 
        names used in the initial global state from which the MCMC output was produced 
        (or at least they must appear in the same order when arranged alphabetically). */
    public StateStringReader(GlobalState initialState) {
        theState = initialState;
    }
    
    /** Produce a header for the summaries obtained in getSummary */
    public abstract String getHeader();
    
    /** Post-process state, produce summary and represent in a string */
    public abstract String getSummary();
    
    public void run(File inFile, File outFile) {
        
        // Prepare output stream
        PrintWriter out;
        if (outFile==null) {
            out = new PrintWriter(System.out);
        }
        else {
            try {
                out = new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
            }
            catch (java.io.IOException anErr) {
                System.out.println("Warning: output to file "+outFile.getName()+" failed: "
                        +"writing output to console instead.");
                out = new PrintWriter(System.out);
            }
        }
        // Print header
        out.println(getHeader());
         
        try {
            BufferedReader br = new BufferedReader(new FileReader(inFile));
            String line = br.readLine(); //Ignore header-line
            line = br.readLine(); //Read first parameter draw
            while(line!=null) {
                if(!line.startsWith("#")) {
                    try {
                        theState.setValueFromString(line);
                    }
                    catch(MCMC.OutputStringProcessing.InsufficientInformationException anEx) {
                        for(int i=0; i<anEx.partString.size(); i++) {
                            // Recover states using the values of others if necessary
                            anEx.recover.get(i).setFromFullString(theState, anEx.partString.get(i), anEx.subStateNames.get(i));
                        }   
                    }
                    // Print output
                    out.println(getSummary());
                }
                line = br.readLine(); //Read next line
            }
            br.close();
        }
        catch(java.io.IOException anEx) {
            System.out.println("Problem reading input file: "+anEx.getMessage());
        }
        catch(treebase.AlgorithmException anEx) {
            System.out.println(anEx.getMessage());
        }
        
        // Flush output
        if (outFile==null) out.flush();
        else out.close();
        
    }
    
}
