/*
    GibbsKernelWrapper
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
 * Uses given Kernel to perform Gibbs updates for the named state
 */

import treebase.AlgorithmException;

public class GibbsKernelWrapper extends KernelWrapper {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    protected GibbsKernel theKernel;
    
    public GibbsKernelWrapper(GibbsKernel kern, String subStateName) {
        super(subStateName);
        theKernel = kern;
    }
    
    @Override
    public int getNumMovesPerSweep() { return 1; }
    
    @Override
    public void sample(GlobalState theState, GlobalState theBackupState) throws treebase.AlgorithmException {
        theBackupState.setLogLikelihood(theState.getLogLikelihood());
        theState.setLogLikelihood(Double.NaN);     // Note: if the log-likelihood does not
                                                   // change during this move (i.e. because
                                                   // it involves prior parameters), the
                                                   // log-likelihood for theState can 
                                                   // be copied from theBackupState during the move
        theKernel.sample(theState, theBackupState, subStateName);
        theBackupState.setLogLikelihood(Double.NaN);
    }
    
    @Override
    public void sample(GlobalState theState, GlobalState theBackupState, double priorTemperature, double likelihoodTemperature) throws treebase.AlgorithmException {
        theBackupState.setLogLikelihood(theState.getLogLikelihood());
        theState.setLogLikelihood(Double.NaN);     // Note: if the log-likelihood does not
                                                   // change during this move (i.e. because
                                                   // it involves prior parameters), the
                                                   // log-likelihood for theState can 
                                                   // be copied from theBackupState during the move
        theKernel.sample(theState, subStateName, priorTemperature, likelihoodTemperature);
        theBackupState.setLogLikelihood(Double.NaN);
    }

    /* Override if you want to print something to do with latent variables etc */
    @Override
    public String printFinalSummary() {
        return "";
    }
    
    @Override
    public void checkCompatibility(GlobalState x) throws AlgorithmException {
        theKernel.checkCompatibility(x, subStateName);
    }

    @Override
    public void resetRandomEngineSeed() {
        theKernel.resetRandomEngineSeed();
    }

    @Override
    public KernelWrapper makeCopy() {
        return new GibbsKernelWrapper((GibbsKernel)theKernel.makeCopy(), subStateName);
    }
    
    public static abstract class GibbsKernel extends Kernel {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        public abstract void sample(GlobalState theState, String subStateName, double priorTemperature, double likelihoodTemperature) throws treebase.AlgorithmException;
          // Samples from a powered version of the kernel in which any likelihood factors are raised to the power likelihoodTemperature
          // and any prior factors are raised to the power priorTemperature. 
          // Note: this method is not required by Chain. It is needed if you want to use the power posterior method to compute marginal 
          // likelihoods or if you want to use MCMCMC. If you don't want to use either of these two methods, throw an UnsupportedOperationException.
    }
    
    public static abstract class GibbsKernelShortcut extends GibbsKernel {
        // Shortcut provided because copyToOrFromBackup should never be needed in a Gibbs move
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1L;
        
        public GibbsKernelShortcut() {
            name = "Gibbs proposal";
        }
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            throw new UnsupportedOperationException("The method copyToOrFromBackup should never be called for a GibbsKernel.");
        }

    }
    
}
