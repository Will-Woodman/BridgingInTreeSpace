/**
    LikelihoodCalculator

    Copyright (C) 2012  Tom M. W. Nye & Sarah E. Heaps

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

package sequencemodels;

/**
 * Methods for computing likelihoods.
 */

import treebase.AlgorithmError;

public interface LikelihoodCalculator {

    public double logLikelihood(Alignment a, SubstitutionModel s) throws AlgorithmError;
    public double siteLogLikelihood(Alignment a, SubstitutionModel s, int siteNum, double rate) throws AlgorithmError;

    /* Methods for rate heterogeneous calcs */
    public double logLikelihood(Alignment a, SubstitutionModel s, double[] rates) throws AlgorithmError; // Continuous rate calc.
    public double logLikelihood(Alignment a, SubstitutionModel s, double[] relativeRates, double[] rateProbabilities) throws AlgorithmError; // Rate mixture model

}
