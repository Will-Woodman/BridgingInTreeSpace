/*
 * PositiveParameter
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
 * Class representing MCMC state for a single positive real parameter.
 *
 * Priors: Gamma, Log-normal, Exponential, Inverse Gamma, Beta-prime
 *
 * Proposals: Log-normal random walk, multiplier, Gamma 
 */

import simulation.ExponentialDistribution;
import simulation.LogNormalDistribution;
import simulation.GammaDistribution;
import simulation.NormalDistribution;
import treebase.AlgorithmException;
import cern.jet.random.tdouble.DoubleUniform;
import java.util.logging.Level;
import java.util.logging.Logger;

public class PositiveParameter extends RealParameter {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /** Constructor */
    public PositiveParameter(double x, String paramName) {
        super(x, paramName);
        /* Error check: make sure parameter is positive.
         Don't throw an AlgorithmError: this makes everything clunky and means you constantly keep
         having to try and catch, even though you know (for example) your proposal is still positive.
         Prefer this simple warning to full-on exception handling. */
        if (x<0.0) System.out.println("AlgorithmError: negative value sent to constructor for a PositiveParameter");
    }
    
    /** Generate PositiveParameter based on string representation of its value */
    @Override
    public void setValueFromString(String str) throws AlgorithmException, OutputStringProcessing.InsufficientInformationException {
        try {
            double val = Double.parseDouble(str);
            if(val<0) throw new AlgorithmException("Warning: string does not contain positive number.");
            value = val;
        }
        catch(NumberFormatException anEx) {
            throw new AlgorithmException("Warning: string does not contain parsable double.");
        }
    }

    /* ---------------------------------------------------------------------- */
    
    // PRIORS
    
    
    /** Gamma prior */
    static public class GammaPrior implements DensityCalculatorAndSampler {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1L;
        
        /* Data consist of shape and scale of the prior.
         These are fixed at the outset and are the same for all states. */
        protected GammaDistribution sampler;

        /** Constructor */
        public GammaPrior(double shape, double scale) {
            sampler = new GammaDistribution(shape, scale);
        }
        @Override
        public DensityCalculator makeCopy() {
            return new GammaPrior(sampler.getShape(), sampler.getScale());
        }
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            PositiveParameter xx = (PositiveParameter) x.getSubState(subStateName);
            return sampler.logpdf(xx.getValue())*priorTemperature;
        }
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof PositiveParameter)) throw new AlgorithmException("Gamma prior not compatible with given state.");
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public void sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ((PositiveParameter)theState.getSubState(subStateName)).setValue(sampler.sample());
        }

    }
    
    
    /** Exponential prior */
    static public class ExponentialPrior implements DensityCalculatorAndSampler {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1L;
        
        /* Data consist of rate of the prior.
         This is fixed at the outset and is the same for all states. */
        protected ExponentialDistribution sampler;

        /** Constructor */
        public ExponentialPrior(double rate) {
            sampler = new ExponentialDistribution(rate);
        }
        @Override
        public DensityCalculator makeCopy() {
            return new ExponentialPrior(sampler.getRate());
        }
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            PositiveParameter xx = (PositiveParameter) x.getSubState(subStateName);
            return sampler.logpdf(xx.getValue())*priorTemperature;
        }
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof PositiveParameter)) throw new AlgorithmException("Exponential prior not compatible with given state.");
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public void sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ((PositiveParameter)theState.getSubState(subStateName)).setValue(sampler.sample());
        }

    }

    
    /** Log-normal prior */
    static public class LogNormalPrior implements DensityCalculatorAndSampler {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1L;
        
        /* Data consist of mu (log-scale) and sigma (squared shape) of the prior.
         These are fixed at the outset and are the same for all states. */
        protected LogNormalDistribution sampler;

        /** Constructor */
        public LogNormalPrior(double mu, double sigma) {
            sampler = new LogNormalDistribution(mu, sigma);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public DensityCalculator makeCopy() {
            return new LogNormalPrior(sampler.getMu(), sampler.getSigma());
        }
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            PositiveParameter xx = (PositiveParameter) x.getSubState(subStateName);
            return sampler.logpdf(xx.getValue())*priorTemperature;
        }
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            //System.out.println("I am lost here "+sampler.getMu()+" "+ sampler.getSigma());
            return logDensity(x, subStateName, 1.0, 1.0);
        }
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof PositiveParameter)) throw new AlgorithmException("Log-normal prior not compatible with given state.");
        }
        @Override
        public void sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ((PositiveParameter)theState.getSubState(subStateName)).setValue(sampler.sample());
        }

    }
    

    /* ---------------------------------------------------------------------- */
    
    // PROPOSALS

    /** Shared proposal behaviour */
    static abstract public class PositiveParameterProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof PositiveParameter)) throw new AlgorithmException("Positive parameter proposal not compatible with given state.");
        }
        //Proposal density methods used only in the calculation of marginal likelihood Chib
        public abstract double logProposalDensity(double Proposed, double Current) throws AlgorithmException;
        
    }
    
    /** Log normal proposal */
    static public class LogNormalProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of a standard deviation for the underlying normal distribution.
           This is a function of the coefficient of variation of the log-normal
           proposal distribution and can be changed during the chain. */
        protected LogNormalDistribution sampler;
        protected double sigma;

        /** Constructor */
        public LogNormalProposal(double sd) {
            name = "Log-normal random walk proposal";
            sigma = sd;
            sampler = new LogNormalDistribution(0.0,1.0);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new LogNormalProposal(sigma);
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            PositiveParameter curr = (PositiveParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = sampler.sample(Math.log(currv),sigma);
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            return Math.log(proposed) -Math.log(current);
        }
        
        public double logProposalDensity(double current, double proposed) throws AlgorithmException {
            //How do I check that is correct?
            return LogNormalDistribution.logpdf(proposed, Math.log(current), sigma);
        }
        
    }

    /** Log Uniform proposal a.k.a "multiplier"*/
    static public class LogUniformProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;

        protected DoubleUniform sampler;
        protected double sc;

        /** Constructor
         Proposal corresponds to scaling from 1/a to a ie. uniform on log scale between -log a and +log a*/
        public LogUniformProposal(double a) {
            name = "Log-uniform random walk (aka multipler) proposal";
            sc = Math.log(a);
            sampler = new DoubleUniform(-sc, sc, simulation.Random.getEngine());
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler = new DoubleUniform(-sc, sc, simulation.Random.getEngine());
        }
        @Override
        public Kernel makeCopy() {
            return new LogUniformProposal(Math.exp(sc));
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            //return new PositiveParameter(((PositiveParameter)x).getValue()*Math.exp(sampler.nextDouble()), x.name);
            PositiveParameter curr = (PositiveParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = Math.exp(Math.log(currv) + sampler.nextDouble());
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            //return ((PositiveParameter)proposedState).getValue() / ((PositiveParameter)currentState).getValue();
            return Math.log(proposed) - Math.log(current);
        }
         public double logProposalDensity(double current, double proposed) throws AlgorithmException {
            return Math.log(proposed) -Math.log(current);
        }
        
    }


    /** Gamma proposal */
    static public class GammaProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;

        protected double invCoeffVarSq; /* Reciprocal of the coefficient of variation of
                                         the proposal, squared. */
        protected GammaDistribution sampler;

        /** Constructor*/
        public GammaProposal(double a) {
            name = "Gamma random walk proposal";
            sampler = new GammaDistribution(1.0,1.0);
            invCoeffVarSq = a;
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new GammaProposal(invCoeffVarSq);
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            PositiveParameter curr = (PositiveParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = sampler.sample(invCoeffVarSq, curr.getValue()/invCoeffVarSq);
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            return (2.0*invCoeffVarSq-1.0)*(Math.log(current) - Math.log(proposed)) - invCoeffVarSq*(current/proposed - proposed/current);
        } 
        
         public double logProposalDensity(double current, double proposed) throws AlgorithmException {
            return GammaDistribution.logpdf(proposed,invCoeffVarSq,current/invCoeffVarSq);
        }
    }
    
    // INDEPENDENCE SAMPLER PROPOSALS
    
    /** Log normal independence sampler proposal */
    static public class LogNormalIndependenceSamplerProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of a lognormal sampler. Its parameters can be changed 
           during the chain. */
        protected LogNormalDistribution sampler;

        /** Constructor */
        public LogNormalIndependenceSamplerProposal(double mlog, double sdlog) {
            name = "Log-normal independence sampler proposal";
            sampler = new LogNormalDistribution(mlog,sdlog);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new LogNormalIndependenceSamplerProposal(sampler.getMu(), sampler.getSigma());
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            RealParameter curr = (RealParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = sampler.sample();
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            double logx = Math.log(current);
            double logy = Math.log(proposed);
            return -(logx - logy)*(1.0 + (logx + logy - 2.0*sampler.getMu())/(2.0*sampler.getSigma()*sampler.getSigma()));
        }
        
         public double logProposalDensity(double current, double proposed) throws AlgorithmException {
            return sampler.logpdf(proposed);
        }
        
    }
    
    /** Gamma independence sampler proposal */
    static public class GammaIndependenceSamplerProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of a gamma sampler. Its parameters can be changed 
           during the chain. */
        protected GammaDistribution sampler;

        /** Constructor*/
        public GammaIndependenceSamplerProposal(double shape, double scale) {
            name = "Gamma independence sampler proposal";
            sampler = new GammaDistribution(shape, scale);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new GammaIndependenceSamplerProposal(sampler.getShape(), sampler.getScale());
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            RealParameter curr = (RealParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = sampler.sample();
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            double shape = sampler.getShape();
            double scale = sampler.getScale();
            return (shape-1.0)*(Math.log(current)-Math.log(proposed))-(current-proposed)/scale;
        }
        
         public double logProposalDensity(double current, double proposed) throws AlgorithmException {
            return sampler.logpdf(proposed);
        }

    }
    
    /** Exponential independence sampler proposal */
    static public class ExponentialIndependenceSamplerProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of an exponential sampler. Its parameter can be changed 
           during the chain. */
        protected ExponentialDistribution sampler;

        /** Constructor */
        public ExponentialIndependenceSamplerProposal(double rate) {
            name = "Exponential independence sampler proposal";
            sampler = new ExponentialDistribution(rate);
        }

        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new ExponentialIndependenceSamplerProposal(sampler.getRate());
        }
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            RealParameter curr = (RealParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = sampler.sample();
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            double rate = sampler.getRate();
            return -rate*(current-proposed);
        }
        
         public double logProposalDensity(double current, double proposed) throws AlgorithmException {
            return sampler.logpdf(proposed);
        }
        
    }
    

}
