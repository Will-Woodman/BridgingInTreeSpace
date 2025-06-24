/*
    State
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
 * Abstract class representing the state in an MCMC chain i.e. a set of unknowns
 */

import java.util.TreeMap;

public abstract class State implements java.io.Serializable {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 2L;
    
    /* A state always has some data, but in this abstract superclass, the state data are
     completely unspecified. */

    /* Name = a name for this parameter e.g. "mean height" */
    protected final String name; // The name of the parameter
    protected TreeMap<String,State> subStates; // May be empty
    
    /* Constructor must set name */
    protected State(String name) {
        this.name = name;
        subStates = new TreeMap<String,State>();
    }
    
    public String getName() {return name;}
    
    /* For printing MCMC output*/
    public String getHeader() {return getName();} // Header for use in MCMC output file - override if not same as parameter name
    public abstract String getValueString(); // String output containing raw value, no parameter name
    
    /* For recreating states from printed MCMC output*/
    public abstract void setValueFromString(String str) throws treebase.AlgorithmException, OutputStringProcessing.InsufficientInformationException;
    public abstract int getLengthStringRepr(); // Return the number of unbroken blocks of characters required to represent the state through a string 
    
    /* Compute running mean */
    public abstract void updateRunningMean(State runningMean, int runningCount) throws treebase.AlgorithmException; // Update
                      // the state in runningMean using the runningCount and present state
    
    /* Accessing substates */
    public final State getSubState(String name) {
        State state = subStates.get(name);
        if(state!=null) return state;
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            state = pairs.getValue().getSubState(name);
            if(state!=null) return state;
        }
        return null;
    }
    
    //delete this in a minute
     public String getSubStates() {
        State state;
        String str="";
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            str+=pairs.getKey()+" ";
        }
        return str;
    }
    
    /* Accessing substates: this method will only work correctly if you only have 
       one state of class classname */
    public final State getSubState(Class className) {
        java.util.Iterator<State> it = subStates.values().iterator();
        while (it.hasNext()) {
            State topState = it.next();
            if(className.isInstance(topState)) return topState;
            else {
                State state = topState.getSubState(name);
                if(state!=null) return state;
            }
        }
        return null;
    }
    
}

    /* Legacy code -------------------------------------------------------------
     
    // Returns null if named state does not exist; otherwise returns the old value 
    public final State setSubState(String name, State x) {
        if(subStates.containsKey(name)) {
            State state = subStates.put(name, x);
            return state;
        }
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            State state = pairs.getValue().setSubState(name, x);
            if(state!=null) return state;
        }
        return null;
    }
    
    ------------------------------------------------------------------------- */
