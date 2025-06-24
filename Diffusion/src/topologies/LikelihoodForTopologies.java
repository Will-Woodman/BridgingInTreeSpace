
/*
   LikelihoodForTopologies
    Copyright (C) 2017  Tom M. W. Nye

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

package topologies;

import MCMC.DensityCalculator;
import MCMC.GlobalState;

import cern.jet.random.tdouble.DoubleUniform;
import diffbase.BrownianState;
import static diffbase.Simulator.randomWalkTree;
import static diffbase.Simulator.sampleRandomWalkParallel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import simulation.NormalDistribution;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeWithTopologicalOperations;

/**
 * Likelihood function for Brownian motion when only topologies of data are known
 */

public class LikelihoodForTopologies  implements DensityCalculator {

    /* Instance variables -- data */
    final private TopologicalData theData;
    
    /* Instance variables -- related to the walks */
    final private int numSteps, numParticles, numProcessors;
    
    final private double smeared; // This is ((2N-5)!!)^{-1} i.e. 1 / [number of (unrooted) othants]
    
    public LikelihoodForTopologies(TopologicalData td, int nStep, int nPart, int nProc) {
        numSteps = nStep;
        numParticles = nPart;
        numProcessors = nProc;
        theData = td;
        
        double s = 0.0;
        for (int i=1; i<=(td.getNumTaxa()-2); i++) {
            s = s - Math.log(2*i-1);
        }
        smeared = s;
    }

    public LikelihoodForTopologies(ArrayList<Tree> theTrees, int nStep, int nPart, int nProc) {
        numSteps = nStep;
        numParticles = nPart;
        numProcessors = nProc;
        theData = new TopologicalData(theTrees);
        
        double s = 0.0;
        for (int i=1; i<=(theData.getNumTaxa()-2); i++) {
            s = s - Math.log(2*i-1);
        }
        smeared = s;

    }

    public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
        
        if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
            throw new AlgorithmException("LikelihoodForTopologies not compatible with "+subStateName);
        }
        BrownianStateForTopologies bs = (BrownianStateForTopologies) x.getSubState(subStateName);
        Tree x0 = bs.getx0().getTree();
        double t0 = bs.gett0();
        
        /* Forward simulate particles */
        HashSet<Tree> particles = new HashSet();
        if (numProcessors>1) {
            particles = sampleRandomWalkParallel(x0, t0, numSteps, false, true, numParticles, numProcessors);
        }
        else {
            NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(t0/numSteps));
            DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
            for (int i=0; i<numParticles; i++) {
                TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(x0);
                randomWalkTree(t, numSteps, false, true, norm, unif);
                particles.add(t);
            }
        }
        
        /* Count topologies in the HashSet */
        HashMap<String, Integer> counts = new HashMap();
        Iterator<Tree> it = particles.iterator();
        while (it.hasNext()) {
            Tree t = it.next();
            int newCount = -1;
            String s = t.toTopologyString();
            
            for (Map.Entry<String, Integer> entry : counts.entrySet()) {
                if (s.equals(entry.getKey())) {
                    // This topology has been seen before
                    newCount = entry.getValue().intValue()+1;
                    break;
                }
            }
            
            if (newCount<0) {
                // new topology
                counts.put(s, 1);
            }
            else {
                // Existing topology
                counts.remove(s);
                counts.put(s, newCount);
            }
        }
        
        /* Loop thru' hash set */
        double loglike = 0.0;
        
        int dataPointsFound = 0;
        for (Map.Entry<String, Integer> entry : counts.entrySet()) {
//            double p = (Math.exp(smeared)+(double)entry.getValue().intValue())/((double)(numParticles+1)); // Probability for that topology
           double p = ((double)entry.getValue().intValue())/((double)(numParticles+1)); // Probability for that topology
            int n = theData.getCount(entry.getKey());
            loglike += n*Math.log(p);
            dataPointsFound += n;
        }
        
        // Add on term corresponding to "smeared" particle -- one particle in all othants
        loglike += (theData.totalCount-dataPointsFound)*(smeared-Math.log(numParticles+1));
                
//        x.setLogLikelihood(0);
//          return 0.0;
        x.setLogLikelihood(loglike);
        return loglike;
    }

    public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
        double l = logDensity(x,subStateName);
        return l*likelihoodTemperature;
    }

    public DensityCalculator makeCopy() {
        return new LikelihoodForTopologies(theData, numSteps, numParticles, numProcessors);
    }

    public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
        if (!(x.getSubState(subStateName) instanceof BrownianStateForTopologies)) {
            throw new AlgorithmException("LikelihoodForTopologies not compatible with given state.");
        }
    }
    
}
