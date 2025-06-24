/*
 * RealParameter
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
 * Class representing a single real parameter
 */

import simulation.NormalDistribution;
import treebase.AlgorithmException;

public class RealParameter extends State{
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /* Data: well, a single real number! */
    protected double value;
 
    /** Constructor */
    public RealParameter(double x, String s) {
        super(s);
        value = x;
    }
    
    /** Provide access to the parameter */
    public double getValue() {
        return value;
    }
    
    /** Set the value
     This method is a bit naughty, but it's useful if you want to create a single 
     RealParameter Object and use it many times -- saves memory!
     Use with caution! */
    public void setValue(double x) {
        value = x;
    }

    /** Get string with parameter value */
    @Override
    public String getValueString() {
        return String.format("%7.7f", value);
    }
    
    /** Set the RealParameter based on string representation of its value */
    @Override
    public void setValueFromString(String str) throws AlgorithmException, OutputStringProcessing.InsufficientInformationException {
        try {
            value = Double.parseDouble(str);
        }
        catch(NumberFormatException anEx) {
            throw new AlgorithmException("Warning: string does not contain parsable double.");
        }
    }
    /** Return the number of unbroken blocks of characters required to represent a RealParameter through a string */
    @Override
    public int getLengthStringRepr() {
        return 1;
    }

    /** Output */
    public String toString() {
        return String.format(name+" = %7.7f", value);
    }
    
    @Override
    public void updateRunningMean(State runningMean, int runningCount) throws AlgorithmException {
        RealParameter x = (RealParameter) runningMean;
        x.setValue((x.getValue()*runningCount+value)/(runningCount+1));
    }

    /** Normal prior */
    static public class NormalPrior implements DensityCalculatorAndSampler {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 1L;
        
        /* Data consist of mean and variance of the prior.
         These are fixed at the outset and are the same for all states. */
        protected NormalDistribution sampler;

        /** Constructor */
        public NormalPrior(double mean, double sd) {
            sampler = new NormalDistribution(mean, sd);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public DensityCalculator makeCopy() {
            return new NormalPrior(sampler.getMu(), sampler.getSigma());
        }
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            RealParameter xx = (RealParameter) x.getSubState(subStateName);
            return sampler.logpdf(xx.getValue())*priorTemperature;
        }
        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {       
            if (!(x.getSubState(subStateName) instanceof RealParameter)) throw new AlgorithmException("Gaussian prior not compatible with given state.");
        }

        @Override
        public void sample(GlobalState theState, String subStateName) throws AlgorithmException {
            ((RealParameter)theState.getSubState(subStateName)).setValue(sampler.sample());
        }

    }
    
    static abstract public class RealParameterProposal extends Kernel {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof RealParameter)) throw new AlgorithmException("Real parameter proposal not compatible with given state.");
        }
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            // whichWay = true if state->backup; otherwise false
            RealParameter x = (RealParameter) theState.getSubState(subStateName);
            RealParameter b = (RealParameter) theBackup.getSubState(subStateName);
            if(whichWay) b.setValue(x.getValue());
            else x.setValue(b.getValue());
        }
        public abstract double logDensityRatio(double current, double proposed) throws AlgorithmException;
        
        public abstract double logProposalDensity(double current, double proposed) throws AlgorithmException;
        
    }

    /** Gaussian proposal */
    static public class NormalProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of a standard deviation.
         This can be changed during the chain. */
        protected NormalDistribution sampler;

        /** Constructor */
        public NormalProposal(double sd) {
            name = "Gaussian random walk proposal";
            sampler = new NormalDistribution(0.0, sd);
        }
        
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            RealParameter curr = (RealParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = currv+sampler.sample();
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new NormalProposal(sampler.getSigma());
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            return 0.0; // Symmetric proposal
        }
        public double logProposalDensity(double Proposed, double Current) throws AlgorithmException {
            return sampler.logpdf(Proposed - Current);
        }

    }
    
    /** Gaussian independence sampler proposal */
    static public class NormalIndependenceSamplerProposal extends RealParameterProposal {
        
        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 2L;
        
        /* Data consists of a mean and standard deviation.
           These can be changed during the chain. */
        protected NormalDistribution sampler;

        /** Constructor */
        public NormalIndependenceSamplerProposal(double mean, double sd) {
            name = "Gaussian independence sampler proposal";
            sampler = new NormalDistribution(mean, sd);
        }
        
        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            RealParameter curr = (RealParameter) theState.getSubState(subStateName);
            double currv = curr.getValue();
            double propv = sampler.sample();
            curr.setValue(propv);
            return logDensityRatio(currv, propv);
        }
        @Override
        public void resetRandomEngineSeed() {
            sampler.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new NormalIndependenceSamplerProposal(sampler.getMu(), sampler.getSigma());
        }
        public double logDensityRatio(double current, double proposed) throws AlgorithmException {
            double mean = sampler.getMu();
            double sd = sampler.getSigma();
            return -((current-mean)*(current-mean)-(proposed-mean)*(proposed-mean))/(2.0*sd*sd);
        }
        
         public double logProposalDensity(double Proposed, double Current) throws AlgorithmException {
            return sampler.logpdf(Proposed);
        }
        
    }

    public static void main(String[] args) {
        RealParameter x = new RealParameter(0.6,"realNumber");
        try {
            // Create byte-reader
            java.io.FileOutputStream fileOut = new java.io.FileOutputStream("/home/sarah/Desktop/realNumber.ser");
            // Wrap byte-reader to translate to object
            java.io.ObjectOutputStream out = new java.io.ObjectOutputStream(fileOut);
            out.writeObject(x);
            out.close();
            fileOut.close();
        } catch (java.io.IOException anEx) {
            System.out.println(anEx.getStackTrace());
        }
        
        RealParameter y = null;
        try {
            java.io.FileInputStream fileIn = new java.io.FileInputStream("/home/sarah/Desktop/realNumber.ser");
            java.io.ObjectInputStream in = new java.io.ObjectInputStream(fileIn);
            y = (RealParameter) in.readObject();
            in.close();
            fileIn.close();
        } catch (java.io.IOException anEx) {
            System.out.println(anEx.getStackTrace());
        } catch (ClassNotFoundException anEx) {
            System.out.println("RealParameter class not found.");
            System.out.println(anEx.getStackTrace());
        }
        System.out.println("Deserialized RealParameter...");
        System.out.println(y.getHeader());
        System.out.println(y.getValueString());
    }

}
