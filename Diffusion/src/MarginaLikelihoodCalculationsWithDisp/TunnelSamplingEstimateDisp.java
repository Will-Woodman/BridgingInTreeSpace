/*
TunnelSamplingEstimateDisp
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
import java.util.Collections;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 */
public class TunnelSamplingEstimateDisp {

    /**
     Main class to calculate the tunnel sampling estimate of the marginal likelihood after running the 
     * sampler, in the unknown dispersion case
     */

    public static void main(String[] args) throws IOException {

        double initialVal =0.0001;//initial value for the iterative estimate
        int numIts  =100;//number of iterations        
        System.out.println("Tunnel");
        String PropSimpleFileName = args[0];
        String PosteriorFileName = args[1];
        String PropsFileName = args[2];
        String AdditionalT0FileName = args[3];
        //get the estimate
        System.out.println(Math.exp(overallEst(PosteriorFileName,PropsFileName,PropSimpleFileName,AdditionalT0FileName,initialVal,numIts)));

    }
    
    //if there is an addional t0 file specified
    private static double overallEst(String posteriorFilename,String propsFilename,String PropSimpleFilename,String additionalT0Filename,double initial,int its) throws IOException{
        return overallEst(posteriorFilename,propsFilename,PropSimpleFilename,true,additionalT0Filename,initial,its);
    
    }
    
     //if there is no addional t0 file specified
     private static double overallEst(String posteriorFilename,String propsFilename,String PropSimpleFilename,double initial,int its) throws IOException{
        return overallEst(posteriorFilename,propsFilename,PropSimpleFilename,false,"",initial, its);
    
    }
    
   //get the tunnel sampling estimate from the files
   private static double overallEst(String posteriorFilename,String propsFilename,String PropSimpleFilename,boolean hasAdditionalT0File,String additionalT0Filename,double initial,int its) throws IOException{
        //get the estimate of the normalising constant of the reference distribution
        //File PropSimpleFile = new File(PropSimpleFilename);
       // Double PropSimple=getPropSimpleFromFile(PropSimpleFile);
        Double PropSimple=1.0;
        //Double PropSimple =1.0;
        File PostFile = new File(posteriorFilename);
        File PropsFile = new File(propsFilename);
        ArrayList<ArrayList<Double>> PropData= readTheFile(PropsFile);
        ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile);
        ArrayList<Double> PropsPropDens = PropData.get(1);
        ArrayList<Double> PropsLogLike = PropData.get(0);
        ArrayList<Double> PropsT0RefDens = PropData.get(2);
        ArrayList<Double> PropsT0PriorDens = PropData.get(3);
        
        ArrayList<Double> PostPropDens = PostData.get(1);
        ArrayList<Double> PostLogLike = PostData.get(0);
        ArrayList<Double> PostT0RefDens = PostData.get(2);
        ArrayList<Double> PostT0PriorDens = PostData.get(3);
        
        if(hasAdditionalT0File){

            File AdditionalT0File = new File(additionalT0Filename);
            PostT0RefDens = readAdditionalT0File(AdditionalT0File);
        }
        
        //combine the prop density, t0 ref density, log likelihood and t0 prior for the estimate        
        ArrayList<Double> PropsValues = normaliseValuesDisp(PropsPropDens,PropsLogLike,PropsT0RefDens,PropsT0PriorDens);
        ArrayList<Double> PostValues = normaliseValuesDisp(PostPropDens,PostLogLike,PostT0RefDens,PostT0PriorDens);
        
        //subtract the median value from the posterior for numerical stability
        Double lStar = returnMedian(PostValues);
        
        for(int i =0; i<PropsValues.size();i++){
            Double temp = PropsValues.get(i)-lStar;
            PropsValues.set(i, temp);
        }
        
        for(int i =0; i<PostValues.size();i++){
            Double temp = PostValues.get(i)-lStar;
            PostValues.set(i, temp);
        }
        
        double out=initial;
        for(int i=0 ; i<its;i++){
        out=bridgeSamplingIt(PropsValues,PostValues,out);
        }
        //add back in the median posterior sample and the normalising constant of the ref dist
        return Math.log(out)+Math.log(PropSimple)+lStar;
   } 
    
    
    
  //calculate one iteration of the bridgeSamplingEst
     public static double bridgeSamplingIt(ArrayList<Double> ProposalVals,ArrayList<Double> PosteriorVals, double current){
        double M1 = PosteriorVals.size();
        double M2 = ProposalVals.size();
        double s1 = M1 / (M1+ M2);
        double s2 = M2/(M1+M2);

        double denom=0;
        //calculate the denominator from the posterior values
        for(int i=0; i<M1;i++){
        double Div=1/(s1*Math.exp(PosteriorVals.get(i))+s2*current);
        if(!(Double.isNaN(Div))){
        denom+=Div;
        }
        else{
           System.out.println("NaN where there shouldn't be in the posterior values");
        }
        }   
        //calculate the numerator from the proposal values
        double numerator=0;
        for(int i=0; i<M2;i++){
        double Div=Math.exp(ProposalVals.get(i))/(s1*Math.exp(ProposalVals.get(i))+s2*current);
        if(!(Double.isNaN(Div))){
        numerator+=Div;
        }
    }

        return((double)(M1)/((double)(M2))*Math.exp(Math.log(numerator)-Math.log(denom)));
    }
     
    //get the median for numerical stability 
    private static Double returnMedian(ArrayList<Double> vals){
        int numSamples = vals.size();
        Double median=0.0;
        Collections.sort(vals);
        if(numSamples%2==0){
            int Ind1 = numSamples / 2 - 1;
            int Ind2 = numSamples / 2;
            return (vals.get(Ind1) + vals.get(Ind2)) / 2.0;
        }
        else{
            median = vals.get(numSamples/2);
        }
        
        
        return median;
        
    } 
    
     //divide the log likelihoods*prior by the independence prop densities*t0 ref densities:
     private static ArrayList<Double> normaliseValuesDisp(ArrayList<Double> PropDensity, ArrayList<Double> LogLike,ArrayList<Double> t0RefDist,ArrayList<Double> t0PriorDist){
        ArrayList<Double> normalised = new ArrayList();
        int numSamples = PropDensity.size();
        for(int i=0; i<numSamples;i++ ){
            normalised.add(LogLike.get(i)+t0PriorDist.get(i)-PropDensity.get(i)-t0RefDist.get(i));
            
        }
        
        return normalised;
        
    } 

    
    /*
     read in the information from the main file: bridge prop densities, bridge log likes
     t0 ref dist densities and t0 prior densities
    */
    public static ArrayList<ArrayList<Double>> readTheFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double>> out = new ArrayList();
        
        ArrayList<Double> propDens = new ArrayList();
        ArrayList<Double> logLike = new ArrayList();
        ArrayList<Double> T0RefDens = new ArrayList();
        ArrayList<Double> T0PriorDens = new ArrayList();
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){//stop when we reach lines with the acceptance rates

            propDens.add(Double.parseDouble(ss[0]));
            logLike.add(Double.parseDouble(ss[1]));
            T0RefDens.add(Double.parseDouble(ss[2]));
            T0PriorDens.add(Double.parseDouble(ss[3]));

             
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        out.add(T0RefDens);
        out.add(T0PriorDens);
        
        return(out);
       
    }
    
    //read the info from the additional t0 file -- has the t0 ref dens for the sampled t0 values in posterior
    public static ArrayList<Double> readAdditionalT0File(File theFile) throws FileNotFoundException, IOException{

        ArrayList<Double> T0RefDens = new ArrayList();
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {
            T0RefDens.add(Double.parseDouble(s));
            }
        
        return T0RefDens;
       
    }
    
        
    //class to calculate the normalising constant of the ref dist -- just need to sum up the values in the file
    private static Double getPropSimpleFromFile(File theFile) throws FileNotFoundException, IOException{
        Double out = 0.0;
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;

        while ((s = br.readLine())!=null) {

            String[] ss = s.split(" ");
            out+=Double.parseDouble(ss[0]);
            
            }

        
        return(out);
       
    }
    
}  

/* old main class for looping through many estimates
public static void main(String[] args) throws IOException {

        double initialVal =0.0001;//initial value for the iterative estimate
        int numIts  =100;//number of iterations        
        System.out.println("Tunnel");
        for(int j=0;j<100;j++){
        //loop over the different files to estimate the variance of the estimator

        
        String PropSimpleFileName = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_PropSimple"+j+".txt";
        String PosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_Post"+j+".txt";
        String PropsFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_Props"+j+".txt";
        String AdditionalT0FileName = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_NewT0RefDist"+j+".txt";
        //get the estimate
        System.out.println(Math.exp(overallEst(PosteriorFileName,PropsFileName,PropSimpleFileName,AdditionalT0FileName,initialVal,numIts)));

        }
    }

*/