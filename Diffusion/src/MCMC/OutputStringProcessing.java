/*
 * OutputStringProcessing
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

import java.util.ArrayList;

/**
 * Utility class for recreating States from printed MCMC output. If a state cannot
 * be rebuilt from its MCMC string representation but, instead, requires additional
 * information from another state, its setFromString() method can throw an
 * InsufficientInformationException. Provision of an instance or static field which
 * is a RecoverStringInformation object, with instructions for how to recreate the 
 * object given information on the global state, then allows reconstruction of 
 * the state.
 */
public class OutputStringProcessing {
        
    public static class InsufficientInformationException extends Exception {
    
        public ArrayList<RecoverStringInformation> recover;
        public ArrayList<String> partString;
        public ArrayList<String> subStateNames; // Names of states which throw the exceptions
    
        public InsufficientInformationException(RecoverStringInformation r, String string, String subStateName) {
            recover = new ArrayList<RecoverStringInformation>();
            recover.add(r);
            partString = new ArrayList<String>();
            partString.add(string);
            subStateNames = new ArrayList<String>();
            subStateNames.add(subStateName);
        }
        public InsufficientInformationException(ArrayList<RecoverStringInformation> r, ArrayList<String> string, ArrayList<String> subStateNames) {
            recover = new ArrayList<RecoverStringInformation>();
            recover.addAll(r);
            partString = new ArrayList<String>();
            partString.addAll(string);
            this.subStateNames = new ArrayList<String>();
            this.subStateNames.addAll(subStateNames);
        }
    }
    
    public static interface RecoverStringInformation {
        
        public void setFromFullString(GlobalState globalState, String partString, String subStateName) throws treebase.AlgorithmException;
        
    }
    
}