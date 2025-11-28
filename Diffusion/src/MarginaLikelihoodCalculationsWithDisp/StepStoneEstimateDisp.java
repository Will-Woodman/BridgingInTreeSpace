/*
StepStoneEstimateDisp
    Copyright (C) 2025  William M Woodman

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

    Contact the author at:  <w.m.woodman2@ncl.ac.uk>
                           
 */


package MarginaLikelihoodCalculationsWithDisp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class StepStoneEstimateDisp {

    /**
     Main class to calculate the stepping stone estimate of marginal likelihood when t0 is unknown, after
     * running all the simulations
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Step Stone");
        //Data files: 
        String PropSimpleFileName = args[0];
        String PosteriorFileName= args[1];
        File PropSimpleFile = new File(PropSimpleFileName);
        Double PropSimple=getPropSimpleFromFile(PropSimpleFile); //get the estimate for the normalising constant of the ref dist (see thesis)   
        //number of beta_k's sampled
        int numOfFiles=100;
        //calculate and print the estimate

        Double est = overallEst(PosteriorFileName,numOfFiles);
        System.out.println(est+Math.log(PropSimple));
        
    }
    
    private static double overallEst(String posteriorFilename,int numFiles) throws IOException{
        Double[] Est= new Double[numFiles];
        double est=0;
        for(int j=0;j<numFiles;j++){
        File PostFile = new File(posteriorFilename+j+".txt");
        ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile);
        //post data values here are (log like - log propdens)*(beta_j-beta_{j-1}) -- see thesis
        ArrayList<Double> PostLogLike = PostData.get(1);
        Est[j]=0.0;

        double eta=Double.NEGATIVE_INFINITY;
        //Take care of numerical stability by subtracting the highest value
        for(int i=0;i<PostLogLike.size();i++){
        eta=Double.NEGATIVE_INFINITY;
        if(PostLogLike.get(i)>eta){
            eta= PostLogLike.get(i);
        }
        
        }
        //average over the sampled values
        for(int i=0;i<PostLogLike.size();i++){
        Est[j]=Est[j]+Math.exp(PostLogLike.get(i)-eta);
        }
        
        Est[j]=Est[j]/PostLogLike.size();
        est=eta+Math.log(Est[j])+est;//add eta back in

        }
        
        return est;
    }
    
    //read in the information from the files
    public static ArrayList<ArrayList<Double>> readTheFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double>> out = new ArrayList();
        
        
        ArrayList<Double> propDens = new ArrayList();
        ArrayList<Double> logLike = new ArrayList();
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){//stop when we get to the acceptance rates
            propDens.add(Double.parseDouble(ss[0]));
            logLike.add(Double.parseDouble(ss[1]));
             
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        
        return(out);
       
    }
    
    //class to calculate the normalising constant of the ref dist -- just need to sum up the values in the file
    private static Double getPropSimpleFromFile(File theFile) throws FileNotFoundException, IOException{
        Double out = 0.0;
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {

            String[] ss = s.split(" ");
            out+=Double.parseDouble(ss[0]);
            
            }

        
        return(out);
       
    }
    
}

/* //old main class for looping over runs to compute Monte Carlo error
    public static void main(String[] args) throws IOException {
        System.out.println("Step Stone");
        //first read in the data
        for(int k=0;k<100;k++){//loop over the separate stepping stone runs to see the variance of the estimator
        
          //get the estimate for the normalising constant of the ref dist (see thesis)  
            
        //String PropSimpleFileName = "/data/ww24/MarginalLikelihoods/10Taxa20240422/SteppingStone/PropSimpleOut10DataPointsx0"+k+".txt";
        //not needed for 4 taxa:
        //String PropSimpleFileName = "/media/c1032934/3DA3-29D2/MarginalLikelihoods/5Taxa20240503/StepStoneFD/Test"+k+"PropSimple.txt";
        //File PropSimpleFile = new File(PropSimpleFileName);
        //Double PropSimple=getPropSimpleFromFile(PropSimpleFile);    
        Double PropSimple=1.0;    
        //number of beta_k's sampled
        int numOfFiles=100;
        //calculate and print the estimate
        String PosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/StepStone/x1_"+k+"r";
        Double est = overallEst(PosteriorFileName,numOfFiles);
        System.out.println(Math.exp(est+Math.log(PropSimple)));
        }
    }

*/