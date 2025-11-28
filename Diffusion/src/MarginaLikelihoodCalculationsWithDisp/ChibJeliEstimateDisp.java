/*
ChibJeliEstimateDisp
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


/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 */
public class ChibJeliEstimateDisp {

      /**
     Main class to calculate the Chib one block estimate of the marginal likelihood after running the 
     * sampler, in the unknown dispersion case
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Chib one block");
        
        String PropSimpleFileName = args[0];
        String PosteriorFileName = args[1];
        String PropsFileName = args[2];
        String AdditionalT0FileName = args[3];
        
        int propYStars=100;//estimate the log marginal likelihood at every propYstars-th bridge in the posterior
        //calculate the estimate
        double theEst = overallEst(PosteriorFileName,PropsFileName,PropSimpleFileName,AdditionalT0FileName ,propYStars);
        System.out.println(theEst);
        
    }
    
     //if there is an addional t0 file specified
    private static double overallEst(String posteriorFilename,String propsFilename,String PropSimpleFilename,String additionalT0Filename,int prop) throws IOException{
        return overallEst(posteriorFilename,propsFilename,PropSimpleFilename,true,additionalT0Filename, prop);
    
    }
    
     //if there is no addional t0 file specified
     private static double overallEst(String posteriorFilename,String propsFilename,String PropSimpleFilename,int prop) throws IOException{
        return overallEst(posteriorFilename,propsFilename,PropSimpleFilename,false,"", prop);
    
    }
    
      /*
    Calculate the Chib one block estimate of the marginal likelihood after running sampler for unknown t0. Input the
    name of the posterior sample, proposals file, proportion of simple file and additional T0 file
    as well as the interval for the number of (t_0^*,Y^*)'s to use in the estimate (see paper)
    */
   private static double overallEst(String posteriorFilename,String propsFilename,String propSimpleFilename,boolean hasAdditionalT0File,String additionalT0Filename,int prop) throws IOException{
        //get the estimate of the normalising constant of the reference distribution 
        //File PropSimpleFile = new File(propSimpleFilename);
        //Double PropSimple=getPropSimpleFromFile(PropSimpleFile);
        Double PropSimple =1.0;
        //read in the data
        File PostFile = new File(posteriorFilename);
        File PropsFile = new File(propsFilename);
        ArrayList<ArrayList<Double>> PropData= readTheFile(PropsFile);
        ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile);
        ArrayList<Double> PropsPropDens = PropData.get(1);
        ArrayList<Double> PropsLogLike = PropData.get(0);
        ArrayList<Double> PropsT0PropDens = PropData.get(2);
        ArrayList<Double> PropsT0PriorDens = PropData.get(3);

        ArrayList<Double> PostPropDens = PostData.get(1);
        ArrayList<Double> PostLogLike = PostData.get(0);
        ArrayList<Double> PostT0PropDens = PostData.get(2);
        ArrayList<Double> PostT0PriorDens = PostData.get(3);

            if(hasAdditionalT0File){
            File AdditionalT0File = new File(additionalT0Filename);
            PostT0PropDens =TunnelSamplingEstimateDisp.readAdditionalT0File(AdditionalT0File);
        }

       
        
        //get an estimate of the posterior density for each of the posterior Y's
        double postSize =PostPropDens.size();
        double numYStars =0.0;
        double theEsts = 0;
        for(int i=0 ; i<postSize;i++){
        if(i%prop==0){
          numYStars++;
        //compute the estimate of the posterior density at the ith point in the posterior  
        double PostEst=ChibJeliIndEst(i,PropsPropDens, PropsLogLike,PropsT0PriorDens, PropsT0PropDens,PostPropDens, PostLogLike,PostT0PriorDens, PostT0PropDens);
        if(Double.isNaN(PostEst)){
            System.out.println("Estimate should not produce NaNs in this context");
        }
        else{
        //add on the estimate of the marginal likelihood at the ith point in the posterior  
        theEsts+=Math.exp(PostLogLike.get(i)+PostT0PriorDens.get(i)-PostEst);
        }
        }
        }
        theEsts=theEsts/(numYStars);
        return Math.log(theEsts) +Math.log(PropSimple);
   
   }
    
    public static double propRatio( double YPropDens,double YLogLike,double t0PriorDens,double t0PropDens,double YStarPropDens, double YStarLogLike,double t0StarPriorDens,double t0StarPropDens){
        //acceptance ratio of moving to (t_0^*,Y^*) given (t_0,Y)
        double alpha = YStarLogLike+t0StarPriorDens+YPropDens+t0PropDens-YLogLike-t0PriorDens-YStarPropDens-t0StarPropDens;
        if(Double.isNaN(alpha)){
           System.out.println("Error calculating propRatio--investigate this");
        }
        
        if(YStarPropDens==Double.NEGATIVE_INFINITY){
            System.out.println("Should not have zero proposal density in this context -- investigate this");
        }
        return(Math.min(0.0,alpha));
    }
    
    
    public static double ChibJeliIndEst(int i,ArrayList<Double> PropPropDen,ArrayList<Double> PropLogLik,ArrayList<Double> PropT0PriorDen,ArrayList<Double> PropT0PropDen,ArrayList<Double> PostPropDen,ArrayList<Double> PostLogLik ,ArrayList<Double> PostT0PriorDen,ArrayList<Double> PostT0PropDen){
        //compute the Chib one block estimate at an individual value (t_0^*,Y^*)
        double qYStar = PostPropDen.get(i);
        double numerator = 0;
        double M1 = PostPropDen.size();
        for(int j=0; j<M1;j++){
            numerator+=Math.exp(propRatio(PostPropDen.get(j),PostLogLik.get(j),PostT0PriorDen.get(j),PostT0PropDen.get(j),PostPropDen.get(i),PostLogLik.get(i),PostT0PriorDen.get(i),PostT0PropDen.get(i)));
        }
        
        double M2 = PropPropDen.size();

        double denom = 0;
        for(int j=0; j<M2;j++){
            if(PropLogLik.get(j)==Double.NEGATIVE_INFINITY){
             System.out.println("Should not have zero log likelihood in this context -- investigate this");
            }
            else{
            denom+=Math.exp(propRatio(PostPropDen.get(i),PostLogLik.get(i),PostT0PriorDen.get(i),PostT0PropDen.get(i),PropPropDen.get(j),PropLogLik.get(j),PropT0PriorDen.get(j),PropT0PropDen.get(j)));
            }
        }
        

        return(PostT0PropDen.get(i)+qYStar+Math.log(M2)-Math.log(M1)+Math.log(numerator)-Math.log(denom));
    }
    
    
    //read in the information from the files
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
    
    //class to calculate the normalising constant of the ref dist -- just need to sum up the values in the file
    public static Double getPropSimpleFromFile(File theFile) throws FileNotFoundException, IOException{
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

/* //old main class for computing many estimates for Monte Carlo error
    public static void main(String[] args) throws IOException {
        System.out.println("Chib");
        boolean AdditionalT0RefDistFile =true;
        for(int j=0;j<100;j++){
        //loop over the different output files to estimate the variance of the estimator

        
        String PropSimpleFileName = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_PropSimple"+j+".txt";
        String PosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_Post"+j+".txt";
        String PropsFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_Props"+j+".txt";
        String AdditionalT0FileName = "/data/ww24/MarginalLikelihoods/4TaxaWith20250212/Tunnel/x1_NewT0RefDist"+j+".txt";
        
        int propYStars=100;//estimate the log marginal likelihood at every propYstars-th bridge in the posterior
        //calculate the estimate
        double theEst = overallEst(PosteriorFileName,PropsFileName,PropSimpleFileName,AdditionalT0FileName ,propYStars);
        System.out.println(Math.exp(theEst));
        }
    }
*/