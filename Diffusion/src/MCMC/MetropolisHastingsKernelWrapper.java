/*
    MetropolisHastingsKernelWrapper
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
 * Uses its Kernel and DensityCalculator to perform Metropolis Hastings updates
 * for the named state
 */

import cern.jet.random.tdouble.DoubleUniform;
import treebase.AlgorithmException;

public class MetropolisHastingsKernelWrapper extends KernelWrapper {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    protected Kernel theKernel;
    protected DensityCalculator theDRC;
    protected DoubleUniform unif;
    
    protected int acceptanceCount;
    protected int proposalCount;
    
    /* Constructor */
    public MetropolisHastingsKernelWrapper(Kernel kern, DensityCalculator drc, String parameterName) {
        super(parameterName);
        theKernel = kern;
        theDRC = drc;
        unif = new DoubleUniform(simulation.Random.getEngine());
    }
    
    /* Constructor */
    public MetropolisHastingsKernelWrapper(Kernel kern, DensityCalculator priorDRC, DensityCalculator likelihoodDRC, String parameterName) {
        super(parameterName);
        theKernel = kern;
        theDRC = new PriorLikelihoodDensityRatioCalculator(priorDRC, likelihoodDRC);
        unif = new DoubleUniform(simulation.Random.getEngine());
    }
    
    public int getNumMovesPerSweep() { return 1; }
    
    public void incrementAcceptanceCount() {
        acceptanceCount++;
    }

    public void incrementProposalCount() {
        proposalCount++;
    }
    
    public String getName() {
        return theKernel.getName()+" for "+subStateName;
    }
    
    @Override
    public String printFinalSummary() {
        return String.format("#Acceptance rate for "+getName()+" = %7.7f%n",(double)acceptanceCount/proposalCount);
        //return String.format("#Acceptance rate for "+getName()+" = %d %d%n",acceptanceCount,proposalCount);
    }
    
    @Override
    public void sample(GlobalState theState, GlobalState theBackupState) throws treebase.AlgorithmException {
        try {
            sample(theState, theBackupState, 1.0, 1.0);
        } catch( java.lang.UnsupportedOperationException anEx) {
            // Work out target density before the move
            double targetDensityCurr = theDRC.logDensity(theState, subStateName);
            // Copy the variables to be changed to the backup state
            theKernel.copyToOrFromBackup(theState, theBackupState, subStateName, true);
        
            // Make a note of the log-likelihood, then set it to NaN so it will be
            // recomputed for the modified state
            double logLikelihoodCurr = theState.getLogLikelihood();
            theBackupState.setLogLikelihood(logLikelihoodCurr);
            theState.setLogLikelihood(Double.NaN);     // Note: if the log-likelihood does not
                                                       // change during this move (i.e. because
                                                       // it involves prior parameters), the
                                                       // log-likelihood for theState can 
                                                       // be copied from theBackupState during the move
            // Modify the state and compute the proposal ratio
            double proposalRatio = theKernel.sample(theState, theBackupState, subStateName);
            theBackupState.setLogLikelihood(Double.NaN);
        
            // Recompute the target density and compute the target density ratio
            double targetDensityProp = theDRC.logDensity(theState, subStateName);
            double targetDensityRatio = targetDensityProp - targetDensityCurr;
            //System.out.println("tDP "+ targetDensityProp+" tDC "+targetDensityCurr+" pR "+proposalRatio);
            
            // Compute the log acceptance probability       
            double logA = proposalRatio +targetDensityRatio;
            double logU = Math.log(unif.nextDouble());
            if(logU < logA) {
                // Accept
                incrementAcceptanceCount();
            } else {
                // Reject: copy the original contents of the state back from the back-up
                theKernel.copyToOrFromBackup(theState, theBackupState, subStateName, false);
                theState.setLogLikelihood(logLikelihoodCurr);
            }
            incrementProposalCount();
        }
    }

    @Override
    public void sample(GlobalState theState, GlobalState theBackupState, double priorTemperature, double likelihoodTemperature) throws treebase.AlgorithmException {
        
        // Work out target density before the move
        double targetDensityCurr = theDRC.logDensity(theState, subStateName, priorTemperature, likelihoodTemperature);
        // Copy the variables to be changed to the backup state
        theKernel.copyToOrFromBackup(theState, theBackupState, subStateName, true);
        
        // Make a note of the log-likelihood, then set it to NaN so it will be
        // recomputed for the modified state
        double logLikelihoodCurr = theState.getLogLikelihood();
        theBackupState.setLogLikelihood(logLikelihoodCurr);
        theState.setLogLikelihood(Double.NaN);     // Note: if the log-likelihood does not
                                                   // change during this move (i.e. because
                                                   // it involves prior parameters), the
                                                   // log-likelihood for theState can 
                                                   // be copied from theBackupState during the move
        // Modify the state and compute the proposal ratio
        double proposalRatio = theKernel.sample(theState, theBackupState, subStateName);
        theBackupState.setLogLikelihood(Double.NaN);
        
        // Recompute the target density and compute the target density ratio
        double targetDensityProp = theDRC.logDensity(theState, subStateName, priorTemperature, likelihoodTemperature);
        
        double targetDensityRatio = targetDensityProp - targetDensityCurr;
        // Compute the log acceptance probability
        //System.out.println(subStateName+", targetDensRat = "+targetDensityRatio+", proposalRat = "+proposalRatio);
        double logA = proposalRatio +targetDensityRatio;
        double logU = Math.log(unif.nextDouble());
        //System.out.print(getName()+": ");
        if(logU < logA) {
            // Accept
            incrementAcceptanceCount();
        } else {
            // Reject: copy the original contents of the state back from the back-up
            theKernel.copyToOrFromBackup(theState, theBackupState, subStateName, false);
            theState.setLogLikelihood(logLikelihoodCurr);
        }
        incrementProposalCount();

    }

    
    @Override
    public void resetRandomEngineSeed() {
        theKernel.resetRandomEngineSeed();
        unif = new DoubleUniform(simulation.Random.getEngine());
    }
    
    @Override
    public KernelWrapper makeCopy() {
        return new MetropolisHastingsKernelWrapper(theKernel.makeCopy(), theDRC.makeCopy(), subStateName);
    }
    
    @Override
    public void checkCompatibility(GlobalState x) throws AlgorithmException {
        theKernel.checkCompatibility(x, subStateName);
        theDRC.checkCompatibility(x, subStateName);
    }
    
    public static class PriorLikelihoodDensityRatioCalculator implements DensityCalculator {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1;
        
        protected DensityCalculator priorDRC;
        protected DensityCalculator likelihoodDRC;
        
        public PriorLikelihoodDensityRatioCalculator(DensityCalculator priDRC, DensityCalculator likeDRC) {
            priorDRC = priDRC;
            likelihoodDRC = likeDRC.makeCopy();
        }
        
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return priorDRC.logDensity(x, subStateName) + likelihoodDRC.logDensity(x, subStateName);
        }
        
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            return priorDRC.logDensity(x, subStateName, priorTemperature, likelihoodTemperature) + likelihoodDRC.logDensity(x, subStateName, priorTemperature, likelihoodTemperature);
        }

        @Override
        public DensityCalculator makeCopy() {
            return new PriorLikelihoodDensityRatioCalculator(priorDRC.makeCopy(), likelihoodDRC.makeCopy());
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            priorDRC.checkCompatibility(x, subStateName);
            likelihoodDRC.checkCompatibility(x, subStateName);
        }
        
    }
    
}
