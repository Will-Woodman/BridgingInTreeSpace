/*
ChibJeliEstimate
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

package MarginalLikelihoodCalculations;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 */
public class ChibJeliEstimate {

    /**
     Main class to calculate the Chib one block estimate of the marginal likelihood after running the 
     * sampler, in the fixed dispersion case
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Chib estimate");
        //work out the number of data points
        String PropsFileName = args[0];
        File PropsFile = new File(PropsFileName);
        String PosteriorFileName=args[1];
        int n = getNumDataPoints(PropsFile);
        int propYstars=500;//estimate the log marginal likelihood at every propYstars-th bridge in the posterior
        
        //calculate the estimate
        double theEst = overallEst(PosteriorFileName,PropsFileName,n,propYstars);
        System.out.println(theEst);
    }
    
    /*
    Calculate the Chib one block estimate of the marginal likelihood after running sampler. Input the
    name of the sample posterior and proposals file as well as the number of data points and the interval
    for the number of Y^*'s to use in the estimate (see paper)
    */
    public static double overallEst(String posteriorFilename,String PropsFilename,int n,int prop) throws IOException{
        double theEst =0;//to store the overall estimate
        Double[] theEsts = new Double[n];//to store the individual estimates for each data point
        
        File PostFile = new File(posteriorFilename);
        File PropsFile = new File(PropsFilename);
        for( int k=0; k<n;k++){//loop over the data points
        theEsts[k]=0.0;
        ArrayList<ArrayList<Double>> PropData= readTheFile(PropsFile,k);
        ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile,k);
        
        ArrayList<Double> PropsPropDens = PropData.get(1);//proposal density of bridges simmed under the prop
        ArrayList<Double> PropsLogLike = PropData.get(0);//log likelihood of bridges simmed under the prop
        ArrayList<Double> PostPropDens = PostData.get(1);//proposal density of bridges simmed under the posterior
        ArrayList<Double> PostLogLike = PostData.get(0);//log likelihood of bridges simmed under the posterior
        
        //get an estimate of the posterior density for each of the posterior Y's
        double postSize =PostPropDens.size();
        double numYStars=0.0;//how many bridges at which we estimate the posterior density
        int numNans =0;//need to check if this still occurs -- it doesn't after taking care of numerical stability
        for(int i=0 ; i<postSize;i++){
           if(i%prop==0){ //get the estimate for the ith value of Y^*
           numYStars++;
        double PostEst=ChibJeliIndEst(i,PropsPropDens, PropsLogLike, PostPropDens, PostLogLike );
        if(Double.isNaN(PostEst)){
            numNans+=1;//this does not appear to have an affect after taking care of numerical stability
        }
        else{
        theEsts[k]+=Math.exp(PostLogLike.get(i)-PostEst);
        }
        }
        }
        theEsts[k]=theEsts[k]/(numYStars);//average over the estimates
        theEst+=Math.log(theEsts[k]);//add on the estimated log marginal likelihood for the current data point
        //System.out.println(theEsts[k]);
        }
        return theEst;
    }
    
    //calculate the proposal acceptance ratio:
    public static double propRatio( double YPropDens,double YLogLike,double YStarPropDens, double YStarLogLike){
        //acceptance ratio of moving to YStar given Y
        double alpha = YStarLogLike+YPropDens-YLogLike-YStarPropDens;
        //System.out.println(alpha);
        if(YStarPropDens==Double.NEGATIVE_INFINITY){
            return(Double.NEGATIVE_INFINITY);
        }
        return(Math.min(0.0,alpha));
    }
    
    //estimate of the posterior density at a given Y^*, th ith bridge
    public static double ChibJeliIndEst(int i,ArrayList<Double> PropPropDen,ArrayList<Double> PropLogLik,ArrayList<Double> PostPropDen,ArrayList<Double> PostLogLik ){
        double qYStar = PostPropDen.get(i);
        double numerator = 0;
        double M1 = PostPropDen.size();
        //calculate the numerator by summing up acceptance rates for Y^* posterior samples.
        for(int j=0; j<M1;j++){
            numerator+=Math.exp(propRatio(PostPropDen.get(j),PostLogLik.get(j),PostPropDen.get(i),PostLogLik.get(i)));
        }
        
        double M2 = PropPropDen.size();
        double denom = 0;
        //calculate the numerator by summing up acceptance rates for independence proposal samples against Y^*.
        for(int j=0; j<M2;j++){
            if(PropLogLik.get(j)==Double.NEGATIVE_INFINITY){
            }
            else{
            denom+=Math.exp(propRatio(PostPropDen.get(i),PostLogLik.get(i),PropPropDen.get(j),PropLogLik.get(j)));
            }
        }
        
        return(qYStar-Math.log(M1)+Math.log(M2)+Math.log(numerator)-Math.log(denom));
    }
    
    
    //read in the information from the files -- only read for each data point at a time because
    //files could be very large for a high number of data points
    public static ArrayList<ArrayList<Double>> readTheFile(File theFile,int k) throws FileNotFoundException, IOException{
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
            if(!(ss[1].equals("rate"))){
            propDens.add(Double.parseDouble(ss[2*k]));//add in the IP density
            logLike.add(Double.parseDouble(ss[2*k+1]));//add in the log likelihood
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        
        return(out);
       
    }
    
      /*auxilliary method to get the number of data points from the direct proposals file -- just counts
    the number of columns in the header and divides by two
    */ 
   private static int getNumDataPoints(File theFile) throws FileNotFoundException, IOException{
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        //introduce t so we don't end up reading last line with the acceptance rates
        int i=0;
        s = br.readLine();
            String[] ss = s.split(" ");
           
        Double out = ss.length/2.0;
        
        return(out.intValue());
       
    } 
    
   /*
   //This is to loop through estimates when repeatedly estimating MLs to find the variance 
   public static void main(String[] args) throws IOException {
        System.out.println("Chib");
        //work out the number of data points
String firstFileName = "/data/ww24/MarginalLikelihoods/10Taxa20240529/Tunnel/PropDataTrees1"+0+".txt";
        File firstFile = new File(firstFileName);
        int n = getNumDataPoints(firstFile);
        for(int j=20;j<21;j++){
        //loop over the different files to estimate the variance of the estimator
        String PosteriorFileName="/data/ww24/MarginalLikelihoods/10Taxa20240529/Tunnel/PostOutTrees1"+j+".txt";
        String PropsFileName="/data/ww24/MarginalLikelihoods/10Taxa20240529/Tunnel/PropDataTrees1"+j+".txt";
        int propYstars=500;//estimate the log marginal likelihood at every propYstars-th bridge in the posterior
        
        //calculate the estimate
        double theEst = overallEst(PosteriorFileName,PropsFileName,n,propYstars);
        System.out.println(Math.exp(theEst));
        }
    }

*/
    
    
}
