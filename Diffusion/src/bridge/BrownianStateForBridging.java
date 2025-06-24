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

package bridge;

/**
    MCMC state representing the initial source tree and dispersion parameter for 
    a Brownian motion in tree-space.
    
    State parameters are 
    * x_0 -- tree;
    * t_0 -- dispersion
    * Array of bridges from x_0 to data points
 */
import MCMC.DensityCalculator;
import MCMC.Kernel;
import MCMC.GlobalState;
import MCMC.KernelWrapper;
import MCMC.MetropolisHastingsKernelWrapper;
import MCMC.PositiveParameter;
import MCMC.SweepKernelWrapper;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.BrownianState;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import simulation.CategoricalDistribution;
import simulation.NormalDistribution;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;
import treebase.Tree;
import treebase.Split;


public class BrownianStateForBridging extends BrownianState {

    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    /** Fixed number of steps -- not inferred */
    public final int numSteps;
    public final int numBridges;
    
    final private static boolean DEBUG = false;

    
    /* CONSTRUCTORS --------------------------------------------------------- */

    public BrownianStateForBridging(TreeAsSplits x0, double t0, ArrayList<TreeAsSplits> dataPoints, boolean constructBridge, ForwardStepBridge template) throws AlgorithmException {
        super(x0,t0);
//        TreeState theTreeState = new TreeState(x0.efficientClone(), "source");
//        PositiveParameter sigSqu = new PositiveParameter(t0, "dispersion");
//        subStates.put(theTreeState.getName(), theTreeState);
//        subStates.put(sigSqu.getName(), sigSqu);
        numSteps = template.numSteps;
        /* Add bridges */
        numBridges = dataPoints.size();
        StarTree.getInstance().setTree(x0);
        
        for (int i=0; i<dataPoints.size(); i++) {
            if (constructBridge) System.out.print("Initializing bridge "+i+"... ");
            BridgeState yi=null;
            int count = 1;
            while (count>0) {
                try {
                    ForwardStepBridge b = template.makeFreshBridge(x0, dataPoints.get(i));
                    yi = new BridgeState(b, getBridgeParamName(i));
                    if (constructBridge) {
                        yi.theBridge.independenceProposal(t0);
                        System.out.println("done in "+count+" attempts.");
                    }
                    break;
                }
                catch (AlgorithmException ex) {
                    count++;
                }
            }
            subStates.put(yi.getName(), yi);
        }
    }

    public BrownianStateForBridging(TreeAsSplits x0, double t0, int m, ArrayList<BridgeState> bridges) {
        super(x0,t0);
//        TreeState theTreeState = new TreeState(x0.efficientClone(), "source");
//        PositiveParameter sigSqu = new PositiveParameter(t0, "dispersion");
//        subStates.put(theTreeState.getName(), theTreeState);
//        subStates.put(sigSqu.getName(), sigSqu);
        numSteps = bridges.get(0).theBridge.getNumSteps();
        /* Add bridges */
        numBridges = bridges.size();
        
        try {
            StarTree.getInstance().setTree(x0);
        } catch (AlgorithmError ex) {
            System.out.println("Problem initializing star tree... Might lead to null pointer crash.");
        }
        System.out.println(numSteps);
        String bridgeName;
        for (int i=0; i<bridges.size(); i++) {
            bridgeName = "bridge_"+String.format("%03d", (i+1)); // 3-digit string for integer
            BridgeState yi = bridges.get(i);
            subStates.put(bridgeName, yi);
        }
    }

    /* UTILITY METHODS ------------------------------------------------------ */
    
    protected BridgeState getBridgeState(int k) {
        return (BridgeState) getSubState(getBridgeParamName(k));
    }

    public int getNumSteps() {
        return numSteps;
    }
    
    public double getLogLike() {
        double l = 0.0;
        for (int i=0; i<numBridges; i++) {
            l += getBridgeState(i).theBridge.getTotalLogLike();
        }
        return l;
    }

    /** Print Newick String
     * @return  */
    @Override
    public String getValueString() {
        if (DEBUG) {
            String s = getx0().toString()+" "+gett0State().getValueString();
            for (int i=0; i<numBridges; i++) {
                //output all the trees on the bridge
                //s += " "+getBridgeState(i).getValueString();
                //some other outputs that could be useful for debugging:
                //s += " "+1;
                //s+=" "+getBridgeState(i).theBridge.getTotalLogLike();
               // s += " "+getBridgeState(i).getTreeList()[getBridgeState(i).getTreeList().length-1];
                //s+= " "+getBridgeState(i).getNumberOfBoundaryCrosses();
                
                
                for(int j=0; j<numSteps;j++){
                    try {
                        s+=" "+getBridgeState(i).theBridge.theTrees[j].getTree().toTopologyString();
                    } catch (AlgorithmError ex) {
                        Logger.getLogger(BrownianStateForBridging.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                
                
            }
            return s;
        }
        
        return getx0().toString()+" "+gett0State().getValueString();
    }


    @Override
    public int getLengthStringRepr() {
        if (DEBUG) return 2+numBridges;
        return 2; // Text = tree string + sd value, so two blocks
    }

    private String getBridgeParamName(int k) {
        return "bridge_"+String.format("%03d", (k+1)); // 3-digit string for integer
    }
    
    public String getHeader() {
        if (DEBUG) {
            String s = getx0State().getName()+" "+gett0State().getName();
            for (int i=0; i<numBridges; i++) {
                s += " "+getBridgeState(i).getName();
            }
            return s;
        }
        return getx0State().getName()+" "+gett0State().getName();
    }
    
    protected void setBridge(BridgeState bs, int k) {
        this.getBridgeState(k).theBridge = bs.theBridge;
    }
    
     
    /* LIKELIHOOD CALCULATOR ------------------------------------------------ */
        
    static public class BrownianStateLikelihoodCalculator implements DensityCalculator {
        
        public BrownianStateLikelihoodCalculator() {
            // Nothing to do
        }
        
        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("BrownianState likelihood calculator not compatible with given state.");
            }
        }
       
        @Override
        public DensityCalculator makeCopy() {
            return new BrownianStateLikelihoodCalculator();
        }
        
        @Override
        public double logDensity(GlobalState x, String subStateName, double priorTemperature, double likelihoodTemperature) throws AlgorithmException {
            double l = logDensity(x,subStateName);
            return l*likelihoodTemperature;
        }

        public double logDensity(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("BrownianStateForBridging likelihood calculator not compatible with given state.");
            }
            
            return ((BrownianStateForBridging)x.getSubState(subStateName)).getLogLike();
        }

    }
        /* PROPOSALS ------------------------------------------------------------ */
    
    /** Sweep thru' independence proposal for each bridge */
    public SweepKernelWrapper makeIndependenceBridgeSweepProposal() {
        // Construct a kernel wrapper for each bridge
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        for (int i=0; i<numBridges; i++) {
            Kernel prop = new BridgeState.IndependenceProposal();
            BridgeState theBridgeState = this.getBridgeState(i);
            KernelWrapper theWrapper = new MetropolisHastingsKernelWrapper(prop, new BridgeState.BridgeLikelihoodCalculator(), theBridgeState.getName());
            kernelWrapperList.add(theWrapper);
        }
        return new SweepKernelWrapper("Bridge independence sweep", kernelWrapperList);
    }

    /** Sweep thru' partial proposal for each bridge.
     Pass catDist=null if you want length of partial bridge to be uniform */
    public SweepKernelWrapper makePartialBridgeSweepProposal(CategoricalDistribution catDist) {
        // Construct a kernel wrapper for each bridge
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        for (int i=0; i<numBridges; i++) {
            Kernel prop = new BridgeState.PartialBridgeProposal(catDist);
            BridgeState theBridgeState = this.getBridgeState(i);
            KernelWrapper theWrapper = new MetropolisHastingsKernelWrapper(prop, new BridgeState.BridgeLikelihoodCalculator(), theBridgeState.getName());
            kernelWrapperList.add(theWrapper);
        }
        return new SweepKernelWrapper("Partial bridge sweep", kernelWrapperList);
    }

    /** Sweep thru' partial proposal for each bridge.  */
    public SweepKernelWrapper makePartialBridgeSweepProposal(int fixedLen) {
        // Construct a kernel wrapper for each bridge
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        for (int i=0; i<numBridges; i++) {
            Kernel prop = new BridgeState.PartialBridgeProposal(fixedLen);
            BridgeState theBridgeState = this.getBridgeState(i);
            KernelWrapper theWrapper = new MetropolisHastingsKernelWrapper(prop, new BridgeState.BridgeLikelihoodCalculator(), theBridgeState.getName());
            kernelWrapperList.add(theWrapper);
        }
        return new SweepKernelWrapper("Partial bridge sweep", kernelWrapperList);
    }

   /** Sweep thru' partial proposal for each bridge.  */
    public SweepKernelWrapper makePartialBridgeSweepProposal(int a, int b) {
        // Construct a kernel wrapper for each bridge
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();
        for (int i=0; i<numBridges; i++) {
            Kernel prop = new BridgeState.PartialBridgeProposal(a, b);
            BridgeState theBridgeState = this.getBridgeState(i);
            KernelWrapper theWrapper = new MetropolisHastingsKernelWrapper(prop, new BridgeState.BridgeLikelihoodCalculator(), theBridgeState.getName());
            kernelWrapperList.add(theWrapper);
        }
        return new SweepKernelWrapper("Partial bridge sweep", kernelWrapperList);
    }


    /** Log normal t0 proposal, retain existing bridges */
    public class T0RandomWalkProposal extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        private Kernel t0Proposal;
        private double propLogLike;
       
        /** Constructor  */
        public T0RandomWalkProposal(Kernel t0Prop) {
            name = "Brownian state t0 proposal retaining bridges";
            t0Proposal = t0Prop;
        }

        public T0RandomWalkProposal(double sd) {
            name = "Brownian state t0 proposal retaining bridges";
            t0Proposal = new PositiveParameter.LogNormalProposal(sd);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("Bridge independence proposal not compatible with given state in proposal for t0 with indept bridges.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            t0Proposal.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new T0RandomWalkProposal(t0Proposal);
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            BrownianStateForBridging currentBrownianState = (BrownianStateForBridging) theState.getSubState(subStateName);
                       
            /* Get log prop density for existing bridges */
            double cpd = 0.0;
            double cll = 0.0;

            for (int i=0; i<numBridges; i++) {
                ForwardStepBridge b = currentBrownianState.getBridgeState(i).theBridge;
                cpd += b.sumStepLogPropDensity();
                cll += b.getTotalLogLike();
            }
            
            double t0PropRatio = t0Proposal.sample(theState, "dispersion");
            double newt0 = ((PositiveParameter)currentBrownianState.getSubState(PositiveParameter.class)).getValue();
                                    
            /* Update bridges with new value of t0 */
            double ppd = 0.0;
            double pll = 0.0;

            for (int i=0; i<numBridges; i++) {
                ForwardStepBridge propBridge= currentBrownianState.getBridgeState(i).theBridge;
                propBridge.updateForNewt0(newt0);
                ppd += propBridge.sumStepLogPropDensity();
                pll += propBridge.getTotalLogLike();       
            }
            
             
            propLogLike = pll;
                       
            return t0PropRatio+cpd-ppd;    
        }
            
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(propLogLike);
            return ratio;
        }


        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForBridging curr = (BrownianStateForBridging) theState.getSubState(subStateName);
            BrownianStateForBridging back = (BrownianStateForBridging) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store t0 and bridges */
                back.sett0(curr.gett0());
                for (int i=0; i<curr.numBridges; i++){
                    curr.getBridgeState(i).theBridge.copyTo(back.getBridgeState(i).theBridge);
                }
            }
            else {
                /* Restore t0 and bridges */
                curr.sett0(back.gett0());
                for (int i=0; i<curr.numBridges; i++){
                    back.getBridgeState(i).theBridge.copyTo(curr.getBridgeState(i).theBridge);
                }

            }
        }
        
    }
    
    

    /** Direct MVN sample for x0, with bridge hierarchy subsample */
    public class X0DirectMVNProposalWithHierarchySubsample extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        private double propLogLike;
        private double sigma; // SD for a single "simple" perturbation
        private int fixedLength=-1;
        
        /* Random number generators */
        private DoubleUniform unifDist;
        private NormalDistribution normDist;
        private CategoricalDistribution catDist; // Should be supported on 0 .. numSteps-1
        

        /** Constructor  */
        public X0DirectMVNProposalWithHierarchySubsample(double sd, CategoricalDistribution cd) {
            name = "x0 proposal (direct MVN with hierarchy subsample)";
            sigma = sd; 
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            catDist = cd;
            normDist = new NormalDistribution(0.0, 1.0);
        }

        /** Make proposal which always resamples same length. 
         Set length = numSteps to achieve an indept proposal. */
        public X0DirectMVNProposalWithHierarchySubsample(double sd, int theLen) {
            name = "x0 proposal (direct MVN with hierarchy subsample)";
            sigma = sd; 
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            normDist = new NormalDistribution(0.0, 1.0);
            fixedLength = theLen;
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("Bridge independence proposal not compatible with given state in proposal for t0 with indept bridges.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            normDist.resetRandomEngineSeed();
            unifDist = new DoubleUniform(simulation.Random.getEngine());
        }
        @Override
        public Kernel makeCopy() {
            throw new UnsupportedOperationException("Request to copy a BrownianState proposal. This wasn't expected.");
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {
            
          

            BrownianStateForBridging currentBrownianState = (BrownianStateForBridging) theState.getSubState(subStateName);
            double t0 = ((PositiveParameter)currentBrownianState.getSubState(PositiveParameter.class)).getValue();
                                    
            //print out some info for understanding the MCMC
            TreeAsSplits currentx0= currentBrownianState.getx0().clone();
            //back to actual code    
            /* How long to resample the bridge? */
            int bridgeLength = fixedLength;
            if (bridgeLength<0) bridgeLength = catDist.sample();
            //System.out.println(bridgeLength);
            /* Sample new x0 */
            
             TreeAsSplits x0 = currentBrownianState.getx0();

            TreeDistributions.sampleDirectMVN(x0, sigma, normDist, unifDist, new double[3]); // Discard logpropdensity -- not needed, since proposal is symmetric
            // currentBrownianState.setx0(x0); // Unnecessary, since x0 *is* current tree
            
            
            /* Make new bridges, conditional on new value of x0 */
            double ppd = 0.0;
            double pll = 0.0;
            double cpd = 0.0;

            for (int i=0; i<numBridges; i++) {
                ForwardStepBridge theBridge = currentBrownianState.getBridgeState(i).theBridge;
                
                cpd += theBridge.computePartialBridgeLogDensity(t0, 1, bridgeLength+1);
                
                try {
                    ppd += theBridge.updateForNewx0(x0, t0, bridgeLength);
                }
                catch (AlgorithmException ex) {
                    // Reject -- bridge failed
                    propLogLike = theState.getLogLikelihood();
                    return Double.NEGATIVE_INFINITY;
                }
                
                pll += theBridge.getTotalLogLike();  
                //System.out.println(theBridge.getNumberOfBoundaryCrosses());
                
                
                //print out some info for understanding the MCMC
                /*
                ArrayList<TreeAsSplits> a = theBridge.getTreePathList();
        String res = "";
        for (int j=0; j<a.size(); j++) {
            res = res + a.get(j).getTree().toTopologyString()+" ";
        }
                
                if(!currentx0.matchesTopology(x0)){
                System.out.println(currentx0.toString()+" "+x0.toString()+" "+currentx0.getTree().toTopologyString()+" "+x0.getTree().toTopologyString()+" "+ bridgeLength+" "+res);
                }
                
                */
                /*
                ArrayList<TreeAsSplits> a = theBridge.getTreePathList();
        String res = "";
        for (int j=0; j<a.size(); j++) {
            res = res + a.get(j).getTree().toTopologyString()+" ";
        }
                
                
                System.out.println(currentx0.toString()+" "+x0.toString()+" "+currentx0.getTree().toTopologyString()+" "+x0.getTree().toTopologyString()+" "+ bridgeLength+" "+res);
                
                
                */
                
                //replace the above with the exact distribution on the 3 spider(just for 4 taxa).
                /*
                double toAdd =0;
                String endpointString= currentBrownianState.getBridgeState(i).getTreeList()[currentBrownianState.getBridgeState(i).getTreeList().length-1];
                TreeAsSplits endpoint= new TreeAsSplits( new Tree(endpointString));
                if(x0.matchesTopology(endpoint)){
                    Split theSplit = x0.getNonTrivialSplits().iterator().next();
                    double theLength1=x0.getSplitLength(theSplit)-endpoint.getSplitLength(theSplit);
                    double theLength2=x0.getSplitLength(theSplit)-endpoint.getSplitLength(theSplit);
                    pll += Math.log(normDist.pdf(theLength1,0,t0)-0.333*normDist.pdf(theLength2,0,t0));
                }
                else{
                    Split theSplit1 = x0.getNonTrivialSplits().iterator().next();
                    Split theSplit2=endpoint.getNonTrivialSplits().iterator().next();
                    double theLength=x0.getSplitLength(theSplit1)+endpoint.getSplitLength(theSplit2);
                    pll += Math.log(0.6667)+normDist.logpdf(theLength,0,t0);
                }    
                
                    */
                
               
            }
            
            propLogLike = pll;
            
            return cpd - ppd; // + x0PropRatio    
            
        }
            
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(propLogLike);
            return ratio;
        }

        
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForBridging curr = (BrownianStateForBridging) theState.getSubState(subStateName);
            BrownianStateForBridging back = (BrownianStateForBridging) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store x0 and bridges */
                back.setx0(curr.getx0().efficientClone());
                for (int i=0; i<curr.numBridges; i++){
                    curr.getBridgeState(i).theBridge.copyTo(back.getBridgeState(i).theBridge);
                }
 
            }
            else {
                /* Restore x0 and bridges */
                curr.setx0(back.getx0());
                for (int i=0; i<curr.numBridges; i++){
                    back.getBridgeState(i).theBridge.copyTo(curr.getBridgeState(i).theBridge);
                }
            }
        }

    }
    
/** Direct MVN sample for x0, with bridge hierarchy subsample -- do not recommend using as this gives
very poor mixing */
    public class X0CGDProposalWithHierarchySubsample extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        private double propLogLike;
        private double sigma; // SD for a single "simple" perturbation
        private int fixedLength=-1;
        
        /* Random number generators */
        private DoubleUniform unifDist;
        private NormalDistribution normDist;
        private CategoricalDistribution catDist; // Should be supported on 0 .. numSteps-1
        

        /** Constructor  */
        public X0CGDProposalWithHierarchySubsample(double sd, CategoricalDistribution cd) {
            name = "x0 proposal (direct MVN with hierarchy subsample)";
            sigma = sd; 
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            catDist = cd;
            normDist = new NormalDistribution(0.0, 1.0);
        }

        /** Make proposal which always resamples same length. 
         Set length = numSteps to achieve an indept proposal. */
        public X0CGDProposalWithHierarchySubsample(double sd, int theLen) {
            name = "x0 proposal (direct MVN with hierarchy subsample)";
            sigma = sd; 
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            normDist = new NormalDistribution(0.0, 1.0);
            fixedLength = theLen;
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("Bridge independence proposal not compatible with given state in proposal for t0 with indept bridges.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            normDist.resetRandomEngineSeed();
            unifDist = new DoubleUniform(simulation.Random.getEngine());
        }
        @Override
        public Kernel makeCopy() {
            throw new UnsupportedOperationException("Request to copy a BrownianState proposal. This wasn't expected.");
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            BrownianStateForBridging currentBrownianState = (BrownianStateForBridging) theState.getSubState(subStateName);
            double t0 = ((PositiveParameter)currentBrownianState.getSubState(PositiveParameter.class)).getValue();
                                    
            /* How long to resample the bridge? */
            int bridgeLength = fixedLength;
            //Will taking out to try sampling the catDist within the for loop
            //if (bridgeLength<0) bridgeLength = catDist.sample();
            
            //Sample new x0 and get x0 proposal densities
            double ppd = 0.0;
            double cpd = 0.0;
            // - get proposal density of current x0 - cpd+= currentBrownianState.getx0State().
            
            /* Sample new x0 */
            cpd+= logCentralDensity(currentBrownianState.getx0(),sigma);
            TreeAsSplits x0 = ForwardStepBridge.sampleCentralGaussian(sigma, normDist, unifDist);
            ppd+=logCentralDensity(x0,sigma);
          
            
            
             // Discard logpropdensity -- not needed, since proposal is symmetric - not true for this one
            // currentBrownianState.setx0(x0); // Unnecessary, since x0 *is* current tree
       
            /* Make new bridges, conditional on new value of x0 */
            
            double pll = 0.0;
            

            for (int i=0; i<numBridges; i++) {
                
                ForwardStepBridge theBridge = currentBrownianState.getBridgeState(i).theBridge;
                cpd += theBridge.computePartialBridgeLogDensity(t0, 1, bridgeLength+1);
                
                try {
                    ppd += theBridge.updateForNewx0(x0, t0, bridgeLength);
                }
                catch (AlgorithmException ex) {
                    // Reject -- bridge failed
                    propLogLike = theState.getLogLikelihood();
                    return Double.NEGATIVE_INFINITY;
                }
                
                pll += theBridge.getTotalLogLike();         
               
            }
            
            propLogLike = pll;
            
            return cpd - ppd; // + x0PropRatio    
        }
            
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(propLogLike);
            return ratio;
        }

        
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForBridging curr = (BrownianStateForBridging) theState.getSubState(subStateName);
            BrownianStateForBridging back = (BrownianStateForBridging) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store x0 and bridges */
                back.setx0(curr.getx0().efficientClone());
                for (int i=0; i<curr.numBridges; i++){
                    curr.getBridgeState(i).theBridge.copyTo(back.getBridgeState(i).theBridge);
                }
 
            }
            else {
                /* Restore x0 and bridges */
                curr.setx0(back.getx0());
                for (int i=0; i<curr.numBridges; i++){
                    back.getBridgeState(i).theBridge.copyTo(curr.getBridgeState(i).theBridge);
                }
            }
        }

    }
    
   
    //do not recommend using as mixing is very poor
    public class X0CGDProposalWithNewBridges extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        private double propLogLike;
        private double sigma; // SD for a single "simple" perturbation
        private int fixedLength=-1;
        
        /* Random number generators */
        private DoubleUniform unifDist;
        private NormalDistribution normDist;
        

        /** Constructor  */
        public X0CGDProposalWithNewBridges(double sd) {
            name = "x0 proposal (direct MVN with hierarchy subsample)";
            sigma = sd; 
            unifDist = new DoubleUniform(simulation.Random.getEngine());
            normDist = new NormalDistribution(0.0, 1.0);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("Bridge independence proposal not compatible with given state in proposal for t0 with indept bridges.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            normDist.resetRandomEngineSeed();
            unifDist = new DoubleUniform(simulation.Random.getEngine());
        }
        @Override
        public Kernel makeCopy() {
            throw new UnsupportedOperationException("Request to copy a BrownianState proposal. This wasn't expected.");
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            BrownianStateForBridging currentBrownianState = (BrownianStateForBridging) theState.getSubState(subStateName);
            double t0 = ((PositiveParameter)currentBrownianState.getSubState(PositiveParameter.class)).getValue();
                                    
            /* How long to resample the bridge? */
            //Sample new x0 and get x0 proposal densities
            double ppd = 0.0;
            double cpd = 0.0;
            // - get proposal density of current x0 - cpd+= currentBrownianState.getx0State().
            
            /* Sample new x0 */
            cpd+= logCentralDensity(currentBrownianState.getx0(),sigma);
            TreeAsSplits x0 = ForwardStepBridge.sampleCentralGaussian(sigma, normDist, unifDist);
            ppd+=logCentralDensity(x0,sigma);
          
            
            
             // Discard logpropdensity -- not needed, since proposal is symmetric - not true for this one
            // currentBrownianState.setx0(x0); // Unnecessary, since x0 *is* current tree
       
            /* Make new bridges, conditional on new value of x0 */
            
            double pll = 0.0;
            

            for (int i=0; i<numBridges; i++) {
                
                ForwardStepBridge theBridge = currentBrownianState.getBridgeState(i).theBridge;
                cpd += theBridge.computeIndependenceLogDensity(t0);
                
                try {
                    ppd += theBridge.independenceProposal(t0);
                }
                catch (AlgorithmException ex) {
                    // Reject -- bridge failed
                    propLogLike = theState.getLogLikelihood();
                    return Double.NEGATIVE_INFINITY;
                }
                
                pll += theBridge.getTotalLogLike();         
               
            }
            
            propLogLike = pll;
            
            return cpd - ppd; // + x0PropRatio    
        }
            
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(propLogLike);
            return ratio;
        }

        
        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForBridging curr = (BrownianStateForBridging) theState.getSubState(subStateName);
            BrownianStateForBridging back = (BrownianStateForBridging) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store x0 and bridges */
                back.setx0(curr.getx0().efficientClone());
                for (int i=0; i<curr.numBridges; i++){
                    curr.getBridgeState(i).theBridge.copyTo(back.getBridgeState(i).theBridge);
                }
 
            }
            else {
                /* Restore x0 and bridges */
                curr.setx0(back.getx0());
                for (int i=0; i<curr.numBridges; i++){
                    back.getBridgeState(i).theBridge.copyTo(curr.getBridgeState(i).theBridge);
                }
            }
        }

    }
    
    
    
    private static double logCentralDensity(TreeAsSplits x, double sd) {
        int nprime = x.getNumTaxa()-3;
        double l = 0.0;
        for (int k=1; k<=nprime-2; k++) {
            l -= Math.log(2*k-1);
        }
        l += nprime*Math.log(2.0);
        l -= nprime*Math.log(sd);
        l -= 0.5*x.sumSquaredLengths(true)/sd/sd;
        return l;
    }
    /* OLD CODE ---------------------------------------------------------------*/
    
    /** Propose t0 based on bridges, retain bridges 
     DO NOT USE!!! 
    */
    public class T0ProposalViaBridgeSteps extends Kernel {

        /** Version number for serialization - increment when structural changes are made */
        private static final long serialVersionUID = 0L;
        
        private double propLogLike;
        private NormalDistribution normDist;
       
        /** Constructor  */
        public T0ProposalViaBridgeSteps(double sigma) {
            name = "Brownian state t0 proposal via bridge steps";
            normDist = new NormalDistribution(0.0, sigma);
        }

        @Override
        public void checkCompatibility(GlobalState x, String subStateName) throws AlgorithmException {
            if (!(x.getSubState(subStateName) instanceof BrownianStateForBridging)) {
                throw new AlgorithmException("t0 proposal with partial bridges not compatible with given state.");
            }
        }

        @Override
        public void resetRandomEngineSeed() {
            normDist.resetRandomEngineSeed();
        }
        @Override
        public Kernel makeCopy() {
            return new T0RandomWalkProposal(normDist.getSigma());
        }

        @Override
        public double sample(GlobalState theState, String subStateName) throws AlgorithmException {

            BrownianStateForBridging currentBrownianState = (BrownianStateForBridging) theState.getSubState(subStateName);
                       
            /* Get log prop density for existing bridges */
            double cpd = 0.0;
            double dsqu = 0.0;
            double cll = 0.0;
            double oldt0 = ((PositiveParameter)currentBrownianState.getSubState(PositiveParameter.class)).getValue();

            for (int i=0; i<numBridges; i++) {
                ForwardStepBridge b = currentBrownianState.getBridgeState(i).theBridge;
                cpd += b.sumStepLogPropDensity();
                dsqu += b.sumSquaredDistances();
                cll += b.getTotalLogLike();
            }
            
            double estimatedt0 = dsqu/numBridges/(currentBrownianState.getNumTaxa()-3);
            double proposedt0 = estimatedt0*Math.exp(normDist.sample());
            ((PositiveParameter)currentBrownianState.getSubState("dispersion")).setValue(proposedt0);
                                    
            /* Update bridges with new value of t0 */
            double ppd = 0.0;
            double pll = 0.0;
            for (int i=0; i<numBridges; i++) {
                ForwardStepBridge propBridge= currentBrownianState.getBridgeState(i).theBridge;
                propBridge.updateForNewt0(proposedt0);
                ppd += propBridge.sumStepLogPropDensity();
                pll += propBridge.getTotalLogLike();       
            }

            propLogLike = pll;
            double t0PropRatio = Math.log(proposedt0) - Math.log(oldt0)+normDist.logpdf(Math.log(oldt0/estimatedt0))-normDist.logpdf(Math.log(proposedt0/estimatedt0));
            return t0PropRatio+cpd-ppd;    
        }
            
        
        public double sample(GlobalState theState, GlobalState theBackup, String subStateName) throws treebase.AlgorithmException {
            double ratio =  sample(theState, subStateName);
            theState.setLogLikelihood(propLogLike);
            return ratio;
        }


        @Override
        public void copyToOrFromBackup(GlobalState theState, GlobalState theBackup, String subStateName, boolean whichWay) {
            BrownianStateForBridging curr = (BrownianStateForBridging) theState.getSubState(subStateName);
            BrownianStateForBridging back = (BrownianStateForBridging) theBackup.getSubState(subStateName);
            if(whichWay) {
                /* Store t0 and bridges */
                back.sett0(curr.gett0());
                for (int i=0; i<curr.numBridges; i++){
                    curr.getBridgeState(i).theBridge.copyTo(back.getBridgeState(i).theBridge);
                }
 
            }
            else {
                /* Restore t0 and bridges */
                curr.sett0(back.gett0());
                for (int i=0; i<curr.numBridges; i++){
                    back.getBridgeState(i).theBridge.copyTo(curr.getBridgeState(i).theBridge);
                }

            }
        }
        
    }
    
   


}
