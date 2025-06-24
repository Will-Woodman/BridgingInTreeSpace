/*
    DensityCalculator
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
 * Calculates log densities up to proportionality
 */

import treebase.AlgorithmException;

public interface DensityCalculator extends java.io.Serializable {
    
    public double logDensity(GlobalState x, String subStateName) throws treebase.AlgorithmException; // Returns log p(x) 
    public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws treebase.AlgorithmException; 
         // Returns the  logarithm of a powered version of p(x) in which any likelihood factors are raised to the power likelihoodTemperature
         // and any prior factors are raised to the power priorTemperature. 
         // Note: this method is not required by Chain. It is needed if you want to use the power posterior method to compute marginal 
         // likelihoods or if you want to use MCMCMC. If you don't want to use either of these two methods, throw an UnsupportedOperationException.

    public DensityCalculator makeCopy(); 
    public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException; // Check this is a valid DC for a particular state
    
}
