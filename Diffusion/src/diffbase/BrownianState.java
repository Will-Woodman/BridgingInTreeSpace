/*
   BrownianStateForBridging
    Copyright (C) 2015  Tom M. W. Nye

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

package diffbase;

/**
    MCMC state representing the initial source tree and dispersion parameter for 
    a Brownian motion in tree-space.
    
    State parameters are 
    * x_0 -- tree;
    * t_0 -- dispersion
    
 */


import MCMC.State;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.OutputStringProcessing;
import MCMC.PositiveParameter;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;


public class BrownianState extends State {

    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    final private static boolean DEBUG = true;

    
    /* CONSTRUCTORS --------------------------------------------------------- */

    public BrownianState(TreeAsSplits x0, double t0)  {
        super("brownian_motion_parameters");
        TreeState theTreeState = new TreeState(x0.efficientClone(), "source");
        PositiveParameter sigSqu = new PositiveParameter(t0, "dispersion");
        subStates.put(theTreeState.getName(), theTreeState);
        subStates.put(sigSqu.getName(), sigSqu);
    }


    /* UTILITY METHODS ------------------------------------------------------ */

    protected TreeState getx0State() {
        return (TreeState)getSubState(TreeState.class);
    }
    protected PositiveParameter gett0State() {
        return (PositiveParameter)getSubState(PositiveParameter.class);
    }
    
    public TreeAsSplits getx0() {
        return ((TreeState)getSubState(TreeState.class)).getTree();
    }
    
    public double gett0() {
        return ((PositiveParameter)getSubState(PositiveParameter.class)).getValue();
    }
        

    /** Print Newick String */
    public String getValueString() {        
        return getx0().toString()+" "+gett0State().getValueString();
    }

    @Override
    /* Caution: :make sure number of diffusion steps is the same */
    public void setValueFromString(String s) throws AlgorithmException {
        String str = s.trim();
        int scInd = str.indexOf(";");
        String treeStr = str.substring(0, scInd);
        String sdStr = str.substring(scInd+1);

        getx0State().setValueFromString(treeStr);
        try {
            gett0State().setValueFromString(sdStr.trim());
        } catch (OutputStringProcessing.InsufficientInformationException ex) {
            throw new AlgorithmException("Nonsense error relating to setting t0 state.");
        }
    }

    @Override
    public int getLengthStringRepr() {
        return 2; // Text = tree string + sd value, so two blocks
    }

    @Override
    public void updateRunningMean(State runningMean, int runningCount) throws AlgorithmException {
        throw new UnsupportedOperationException("Not supported.");
    }
   
    public String getHeader() {
        return getx0State().getName()+" "+gett0State().getName();
    }
    
    protected void sett0(double x) {
        ((PositiveParameter)getSubState(PositiveParameter.class)).setValue(x);
    }

    
    protected void setx0(TreeAsSplits x) {
        ((TreeState)getSubState(TreeState.class)).setTree(x);
    }
    
    protected int getNumTaxa() {
        return ((TreeState)getSubState(TreeState.class)).getTree().getNumTaxa();
    }
    
    /* PRIORS  -------------------------------------------------------------- */

    /** Prior based on combining indept priors on tree and variance param */
    static public class SimpleBrownianStatePrior implements DensityCalculator {

        /* Instance variables: a prior for the TreeState and a Prior for the standard deviation */
        private DensityCalculator treePrior, variancePrior;

        /* Constructors */

        /** Construct from a prior for TreeState and a prior for PositiveParameter */
        public SimpleBrownianStatePrior(DensityCalculator t, DensityCalculator v) {
            treePrior = t;
            variancePrior = v;
        }

        /** Construct from a prior for TreeState and provide rate for an exp prior on variance.
         See paper for details. */
        public SimpleBrownianStatePrior(DensityCalculator t, double expRate) {
            treePrior = t;
            if (expRate>0)
                variancePrior = new PositiveParameter.ExponentialPrior(expRate);
            else {
                System.out.println("Warning: setting prior on diffusion variance to be exponential with rate 18.5");
                variancePrior = new PositiveParameter.ExponentialPrior(18.5);
            }
        }

        @Override
        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            return logDensity(x, subStateName, 1.0, 1.0);
        }

        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            State subst = x.getSubState(subStateName);
            double l = treePrior.logDensity(x, subst.getSubState(TreeState.class).getName())+variancePrior.logDensity(x,subst.getSubState(PositiveParameter.class).getName());
            return l*priorTemperature;
        }

        @Override
        public DensityCalculator makeCopy() {
            return new SimpleBrownianStatePrior(treePrior, variancePrior);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if(!(x.getSubState(subStateName) instanceof BrownianState)) {
                throw new AlgorithmException("SimpleBrownianStatePrior not compatible with given state.");
            }
            treePrior.checkCompatibility(x,x.getSubState(subStateName).getSubState(TreeState.class).getName());
            variancePrior.checkCompatibility(x,x.getSubState(subStateName).getSubState(PositiveParameter.class).getName());
        }


    }


}
    

 