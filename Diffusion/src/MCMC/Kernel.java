/*
    Kernel
    Copyright (C) 2013 Sarah E. Heaps

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

/**
 * Generates samples from a distribution
 */

import treebase.AlgorithmException;

public abstract class Kernel implements java.io.Serializable, simulation.RandomEngineSeedSetter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;
    
    protected String name;
    
    public String getName() {
        return name;
    }
    
    /** The sample methods return the log proposal density ratio. In the vast majority of cases, sample should
       simply modify the current state and so the backup state will not be required. The default behaviour is
       therefore to call the simple version of sample which does not require the backup state. However if sample 
       changes an object which is shared by the global and backup state, e.g. vertices or edges in a Tree object, 
       then the backup state will also required by sample. */
    public abstract double sample(GlobalState theState, String subStateName) throws treebase.AlgorithmException;
    public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
        return sample(theState, subStateName);
    }
    
    public abstract void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay); // whichWay = true if state->backup; otherwise false
    public abstract Kernel makeCopy();
    public abstract void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException; // Check this is a valid kernel for a particular state
    
}
