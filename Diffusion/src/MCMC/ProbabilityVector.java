/*
 * ProbabilityVector
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
 * Class representing MCMC state for a probability vector.
 */

import java.util.Arrays;
import simulation.DirichletDistribution;
import treebase.AlgorithmException;

public class ProbabilityVector extends State {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /* Data: a vector of probabilities */
    protected double[] probVector;

    public ProbabilityVector(double[] x, String paramName) {
        super(paramName);
        probVector = new double[x.length];
        System.arraycopy(x,0,probVector,0,x.length);
    }
    
    public ProbabilityVector(int size, String paramName) {
        super(paramName);
        probVector = new double[size];
        Arrays.fill(probVector,1.0/size);
    }

    public int dimension() {
        return probVector.length;
    }

    /** Provide access to the parameter */
    public double[] getValue() {
        double[] x = new double[this.dimension()];
        System.arraycopy(probVector,0,x,0,x.length);
        return x;
    }
    public double getValue(int i) {
        return probVector[i];
    }
    public void setValue(double[] p) {
        System.arraycopy(p,0,probVector,0,probVector.length);
    }

    /** Output */
    public String toString() {
        String s = name+" = [";
        for (int i=0; i<probVector.length; i++) {
            s += String.format("%7.7f",probVector[i]);
            if (i<probVector.length-1) s += ",";
        }
        s += "]";
        return s;
    }
    @Override
    public String getValueString() {
        String s = new String();
        for (int i=0; i<probVector.length; i++) {
            s += String.format("%7.7f",probVector[i])+" ";
        }
        return s.trim();
    }

    @Override
    public String getHeader() {
        String str = new String();
        for(int i=0;i<probVector.length;i++){
            str += name+"["+(i+1)+"] ";
        }
        return str.trim();
    }
    
    /** Input */
    @Override
    public void setValueFromString(String str) throws AlgorithmException {
        String[] strArr = str.split(" ");
        if(strArr.length != probVector.length) throw new AlgorithmException("String does not represent a probability vector of the correct length.");
        double[] val = new double[strArr.length];
        double sum = 0.0;
        try {
            for(int i=0; i<strArr.length; i++) {
                val[i] = Double.parseDouble(strArr[i]);
                if(val[i]<0 || val[i]>1 ) throw new AlgorithmException("String does not represent a probability vector.");
                sum += val[i];
            }
            if(Math.abs(sum-1.0)>1e-5) throw new AlgorithmException("String does not represent a probability vector.");
        }
        catch(NumberFormatException anEx) {
            throw new AlgorithmException("Warning: string does not contain parsable double.");
        }
        System.arraycopy(val, 0, probVector, 0, probVector.length);
    }
    @Override
    public int getLengthStringRepr() {
        return probVector.length;
    }
    
    @Override
    public void updateRunningMean(State runningMean, int runningCount) throws AlgorithmException {
        ProbabilityVector x = (ProbabilityVector) runningMean;
        if(x.dimension() != dimension()) throw new AlgorithmException("Running mean not compatible with current probability vector.");
        for(int i=0; i<dimension(); i++) {
            x.probVector[i] = (x.probVector[i]*runningCount+probVector[i])/(runningCount+1);
        }
    }

    /**************************************************************************/
    
    // PRIORS

    /** Dirichlet prior */
    static public class DirichletPrior implements DensityCalculatorAndSampler {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1L;
        
        /* Data consist of mean probability vector and concentration parameter.
         These are fixed at the outset and are the same for all states. */
        protected DirichletDistribution sampler;

        /** Constructor */
        public DirichletPrior(double[] p, double conc) throws AlgorithmException {
            sampler = new DirichletDistribution(p, conc);
        }
        
        
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public DensityCalculator makeCopy() {
            DensityCalculator p;
            try {
                p = new DirichletPrior(sampler.getMean(), sampler.getConcParam());
            } catch(AlgorithmException anEx) {
                System.out.println("Problem making copy of DirichletPrior. This should "
                        + "not be possible.");
                p = null;
            }
            return p;
        }
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            ProbabilityVector xx = (ProbabilityVector) x.getSubState(subStateName);
            return sampler.logpdf(xx.getValue())*priorTemperature;
        }
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {       
            if (!(x.getSubState(subStateName) instanceof ProbabilityVector)) throw new AlgorithmException("Dirichlet prior not compatible with given state.");
            if (((ProbabilityVector)x.getSubState(subStateName)).dimension() != sampler.getMean().length) throw new AlgorithmException("Dimension of dirichlet prior is not compatible with given state.");
        }
        @Override
        public void sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ((ProbabilityVector)theState.getSubState(subStateName)).setValue(sampler.sample());
        }

    }

    /**************************************************************************/
    
    // PROPOSALS

    /** Dirichlet proposal */
    static public class DirichletProposal extends Kernel {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        protected DirichletDistribution sampler;
        protected double concParam;

        /** Constructor */
        public DirichletProposal(int n, double cp) throws AlgorithmException {
            name = "Dirichlet random walk proposal";
            double[] p = new double[n];
            double delta = 1.0/n;
            for (int i=0; i<n; i++) p[i] = delta;
            sampler = new DirichletDistribution(p, cp);
            concParam = cp;
        }
        
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof ProbabilityVector)) throw new AlgorithmException("Probability vector proposal not compatible with given state.");
            if (((ProbabilityVector)x.getSubState(subStateName)).dimension() != sampler.getMean().length) throw new AlgorithmException("Dimension of dirichlet proposal is not compatible with given state.");
        }
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            // whichWay = true if state->backup; otherwise false
            ProbabilityVector x = (ProbabilityVector) theState.getSubState(subStateName);
            ProbabilityVector b = (ProbabilityVector) theBackup.getSubState(subStateName);
            if(whichWay) b.setValue(x.getValue());
            else x.setValue(b.getValue());
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ProbabilityVector curr = (ProbabilityVector) theState.getSubState(subStateName);
            double[] currv = curr.getValue();
            double[] propv = sampler.sample(currv, concParam);
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            Kernel p;
            try {
                p = new DirichletProposal(sampler.getMean().length, concParam);
            } catch(AlgorithmException anEx) {
                System.out.println("Problem making copy of DirichletProposal. This should "
                        + "not be possible.");
                p = null;
            }
            return p;
        }
        public double logDensityRatio(double[] current, double[] proposed) throws AlgorithmException {
            double s = 0.0;
            for (int i=0; i<current.length; i++) {
                s += (concParam*proposed[i]-1.0)*Math.log(current[i]) + cern.jet.stat.tdouble.Gamma.logGamma(concParam*current[i]) - (concParam*current[i]-1.0)*Math.log(proposed[i]) - cern.jet.stat.tdouble.Gamma.logGamma(concParam*proposed[i]);
            }
            return s;
        }

    }
    
    /** Dirichlet proposal with nudge */
    static public class DirichletProposalWithNudge extends DirichletProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        protected double nudge;

        /** Constructors */        
        public DirichletProposalWithNudge(int n, double cp, double del) throws AlgorithmException {
            super(n, cp);
            name = "Dirichlet random walk proposal with nudge";
            nudge = del;
        }
        
        @Override
        public Kernel makeCopy() {
            Kernel p;
            try {
                p = new DirichletProposalWithNudge(sampler.getMean().length, concParam, nudge);
            } catch(AlgorithmException anEx) {
                System.out.println("Problem making copy of DirichletProposalWithNudge. This should "
                        + "not be possible.");
                p = null;
            }
            return p;
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ProbabilityVector curr = (ProbabilityVector) theState.getSubState(subStateName);
            double[] currv = curr.getValue();
            int n = currv.length;
            double[] p = new double[n];
            for (int i=0; i<n; i++){
                p[i] = (concParam*currv[i]+nudge)/(concParam + n*nudge);
            }
            double[] propv = sampler.sample(p, concParam + n*nudge);
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        @Override
        public double logDensityRatio(double[] current, double[] proposed) throws AlgorithmException {
            double s = 0.0;
            for (int i=0; i<current.length; i++) {
                s += (concParam*proposed[i]+nudge-1.0)*Math.log(current[i]) + cern.jet.stat.tdouble.Gamma.logGamma(concParam*current[i]+nudge);
                s += - (concParam*current[i]+nudge-1.0)*Math.log(proposed[i]) - cern.jet.stat.tdouble.Gamma.logGamma(concParam*proposed[i]+nudge);
            }
            return s;
        }
        
    }
    
    // INDEPENDENCE SAMPLER PROPOSALS

    /** Dirichlet independence sampler proposal */
    static public class DirichletIndependenceSamplerProposal extends Kernel {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of a dirichlet sampler. Its parameters cannot be changed
           during the chain */
        protected DirichletDistribution sampler;

        /** Constructor */
        public DirichletIndependenceSamplerProposal(double[] mean, double cp) throws AlgorithmException {
            name = "Dirichlet independence sampler proposal";
            sampler = new DirichletDistribution(mean, cp);
        }
        
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof ProbabilityVector)) throw new AlgorithmException("Probability vector proposal not compatible with given state.");
            if (((ProbabilityVector)x.getSubState(subStateName)).dimension() != sampler.getMean().length) throw new AlgorithmException("Dimension of dirichlet proposal is not compatible with given state.");
        }
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            // whichWay = true if state->backup; otherwise false
            ProbabilityVector x = (ProbabilityVector) theState.getSubState(subStateName);
            ProbabilityVector b = (ProbabilityVector) theBackup.getSubState(subStateName);
            if(whichWay) b.setValue(x.getValue());
            else x.setValue(b.getValue());
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ProbabilityVector curr = (ProbabilityVector) theState.getSubState(subStateName);
            double[] currv = curr.getValue();
            double[] propv = sampler.sample();
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            Kernel p;
            try {
                p = new DirichletIndependenceSamplerProposal(sampler.getMean(), sampler.getConcParam());
            } catch(AlgorithmException anEx) {
                System.out.println("Problem making copy of DirichletIndependenceSamplerProposal. This should "
                        + "not be possible.");
                p = null;
            }
            return p;
        }
        public double logDensityRatio(double[] current, double[] proposed) throws AlgorithmException {
            double s = 0.0;
            for (int i=0; i<current.length; i++) {
                s += (sampler.getConcParam(i)-1)*(Math.log(current[i])-Math.log(proposed[i]));
            }
            return s;
        }

    }


}
