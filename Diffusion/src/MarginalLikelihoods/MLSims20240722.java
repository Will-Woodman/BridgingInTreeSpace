/*
    InferBrownianParamsMCMC
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

package MarginalLikelihoods;

/**
 *  
 */

import bridge.*;
import MCMC.Chain;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import MCMC.KernelWrapper;
import MCMC.MetropolisHastingsKernelWrapper;
import MCMC.PositiveParameter;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.TreeState;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import simulation.CategoricalDistribution;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;
import bridge.BridgeState;

/*
Run the MCMC needed to compute a stepping stone sampling estimate of the marginal likelihood when t0 is fixed.
Use the quasistatic method so burnin is only needed to be run once. First distribution (indpendence proposals) is
simulated from directly
*/

public class MLSims20240722 {
    

    public static void main(String[] args) throws AlgorithmException {
        
double[] t0={0.125,0.5};
Integer[] treeNames ={1,2,3};
//double[] t0={0.25};
//Integer[] treeNames ={3};
String FolderName="/data/ww24/MarginalLikelihoods/10TaxaOhne20240912/";
String TreeFolderName="/home/c1032934/Documents/Netbeans/TopInf20240115/IndepPropTuning/10Taxa20240323/";

Integer StartPoint = Integer.parseInt(args[0]);
Integer EndPoint = Integer.parseInt(args[1]);

String StartTreeFilename="";
String EndTreeFilename="";

String ChibFilename = "Chib/";
String StepStoneFilename="SteppingStone/";
String Filename="";


for(Integer treeName:treeNames){
for(Double t0Val:t0){
    StartTreeFilename=TreeFolderName+"SourceTree"+treeName+".txt";
    EndTreeFilename=TreeFolderName+"TargetTree"+treeName+".txt";
    Filename= FolderName+ChibFilename+"Trees"+treeName+"Disp"+t0Val+"r";
    String[] args2 = {Filename,StartTreeFilename,EndTreeFilename,StartPoint.toString(),EndPoint.toString(),t0Val.toString()};
    ChibSampler.main(args2);
}
      
    }



for(Integer treeName:treeNames){
for(Double t0Val:t0){
    StartTreeFilename=TreeFolderName+"SourceTree"+treeName+".txt";
    EndTreeFilename=TreeFolderName+"TargetTree"+treeName+".txt";
    Filename= FolderName+StepStoneFilename+"Trees"+treeName+"/Disp"+t0Val+"r";
    System.out.println(Filename);
    String[] args2 = {Filename,StartTreeFilename,EndTreeFilename,StartPoint.toString(),EndPoint.toString(),t0Val.toString()};
    SteppingStoneSampler.main(args2);
}
       
    }


    }

}