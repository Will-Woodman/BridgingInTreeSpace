/*
    KernelWrapper
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
 * Abstract class for performing MCMC updates
 */

import treebase.AlgorithmException;

public abstract class KernelWrapper implements java.io.Serializable, simulation.RandomEngineSeedSetter {
        
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    protected final String subStateName;
    
    protected KernelWrapper(String parameterName) {
        subStateName = parameterName;
    }
    
    public abstract void sample(GlobalState theState, GlobalState theBackupState) throws treebase.AlgorithmException;
    public abstract void sample(GlobalState theState, GlobalState theBackupState, double priorTemperature, double likelihoodTemperature) throws treebase.AlgorithmException;
    public abstract String printFinalSummary(); //For M-H kernels, this will include the acceptance rate. It might also include
                      // summaries of the posterior for latent variables. 
    public abstract void checkCompatibility(GlobalState x) throws AlgorithmException; // Check this is a valid kernel for a particular state
    public abstract void resetRandomEngineSeed();
    public abstract KernelWrapper makeCopy();
    public abstract int getNumMovesPerSweep();
    
}
