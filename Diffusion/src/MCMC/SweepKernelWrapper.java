/*
    SweepKernelWrapper
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
 * Represents a fixed sweep of KernelWrappers
 */

import java.util.ArrayList;
import treebase.AlgorithmError;
import treebase.AlgorithmException;

public class SweepKernelWrapper extends KernelWrapper {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    protected ArrayList<KernelWrapper> subKernelWrappers;
    protected ArrayList<Integer> repeats;
    protected int size;
    protected int proposalIndex;
    protected int repCounter;
    
    /* Constructor for use by extending classes */
    protected SweepKernelWrapper(String parameterName) {
        super(parameterName);
        subKernelWrappers = new ArrayList<KernelWrapper>();
        repeats = new ArrayList<Integer>();
    }
    
    /* Default constuctor: number of moves specified by each kernel-wrapper */
    public SweepKernelWrapper(String parameterName, ArrayList<KernelWrapper> theKernWrappers) {
        super(parameterName);
        subKernelWrappers = new ArrayList<KernelWrapper>();
        subKernelWrappers.addAll(theKernWrappers);
        size = subKernelWrappers.size();
        repeats = new ArrayList<Integer>();
        for(int i=0; i<size; i++) repeats.add(subKernelWrappers.get(i).getNumMovesPerSweep());

        proposalIndex = size-1;
        repCounter = repeats.get(size-1)-1;
    }
    
    /* Constuctor which allows you to specify more than one repeat of each kernel-wrapper */
    //20230221 Made this public for the t0 inference -- couldn't see why it was private
    public SweepKernelWrapper(String parameterName, ArrayList<KernelWrapper> theKernWrappers, int[] rpts) throws AlgorithmError {
        super(parameterName);
        if (rpts.length!=theKernWrappers.size()) throw new AlgorithmError("Mismatch in sweep wrapper between number of wrappers and number of repeats.");
        subKernelWrappers = new ArrayList<KernelWrapper>();
        subKernelWrappers.addAll(theKernWrappers);
        size = subKernelWrappers.size();
        repeats = new ArrayList<Integer>();
        for(int i=0; i<size; i++) repeats.add(rpts[i]);

        proposalIndex = size-1;
        repCounter = repeats.get(size-1)-1;
    }
    
    /** Number of moves specified by each kernel-wrapper */
    public void addKernelWrapper(KernelWrapper theKernWrapper) {
        subKernelWrappers.add(theKernWrapper);
        size++;
        repeats.add(theKernWrapper.getNumMovesPerSweep());
        
        proposalIndex = size-1;
        repCounter = repeats.get(size-1)-1;
    }
    
    /** Allows you to specify more than one repeat of each kernel-wrapper */
    public void addKernelWrapper(KernelWrapper theKernWrapper, int numRpts) {
        subKernelWrappers.add(theKernWrapper);
        size++;
        repeats.add(numRpts*theKernWrapper.getNumMovesPerSweep());
        
        proposalIndex = size-1;
        repCounter = repeats.get(size-1)-1;
    }
    
    
    
    @Override
    public int getNumMovesPerSweep() {
        int num = 0;
        for(int i=0; i<size; i++) num += repeats.get(i);
        return num;
    }
    
    @Override
    public void sample(GlobalState theState, GlobalState theBackupState) throws AlgorithmException {
        repCounter++;
        if (repCounter==repeats.get(proposalIndex)) {
            repCounter = 0;
            proposalIndex++;
        }
        if (proposalIndex==size) proposalIndex=0;

        subKernelWrappers.get(proposalIndex).sample(theState, theBackupState);
    }

    @Override
    public void sample(GlobalState theState, GlobalState theBackupState, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
        repCounter++;
        if (repCounter==repeats.get(proposalIndex)) {
            repCounter = 0;
            proposalIndex++;
        }
        if (proposalIndex==size) proposalIndex=0;

        subKernelWrappers.get(proposalIndex).sample(theState, theBackupState, priorTemperature, likelihoodTemperature);
    }
    
    @Override
    public String printFinalSummary() {
        String str = "";
        for(int i=0; i<subKernelWrappers.size(); i++) str += subKernelWrappers.get(i).printFinalSummary();
        return str;
    }

    @Override
    public void checkCompatibility(GlobalState x) throws AlgorithmException {
        for(int i=0; i<size; i++) {
            subKernelWrappers.get(i).checkCompatibility(x);
        }
    }

    @Override
    public void resetRandomEngineSeed() {
        for(int i=0; i<size; i++) {
            subKernelWrappers.get(i).resetRandomEngineSeed();
        }
    }

    @Override
    public KernelWrapper makeCopy() {
        ArrayList<KernelWrapper> kernWrappers = new ArrayList<KernelWrapper>();
        int[] rpts = new int[size];
        for(int i=0; i<size; i++) {
            kernWrappers.add(subKernelWrappers.get(i).makeCopy());
            rpts[i] = repeats.get(i);
        }
        SweepKernelWrapper ret;
        try {
            ret = new SweepKernelWrapper(subStateName, kernWrappers, rpts);
        } catch(AlgorithmException anEx) {
            System.out.println("This should not be possible");
            ret = null;
        }
        return ret;
    }
    
}
