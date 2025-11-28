/*
ChibTwoBlockEstimate
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


public class ChibTwoBlockEstimate {
    /**
     Main class to calculate the Chib two block estimate of the marginal likelihood after running the 
     * sampler
     */
    public static void main(String[] args) throws IOException {
        
        System.out.println("Chib two block");
        //loop over the different output files to estimate the variance of the estimator
        String FixedT0PosteriorFileName = args[0];
        String FixedT0PropsFileName = args[1];
        String FullPosteriorFileName = args[2];
        int propYStars =100;//estimate the fixed t0 part of the log marginal likelihood at every propYstars-th bridge in the posterior
        
        //calculate the estimate
        double theEst=overallEst(FullPosteriorFileName,FixedT0PosteriorFileName,FixedT0PropsFileName,propYStars);
        System.out.println(theEst);
           
        
    }
    
    private static double overallEst(String fullPosteriorFilename,String fixedPostFilename,String fixedPropsFilename,int prop) throws IOException{
        //first read in the data    
        File FixedT0PostFile = new File(fixedPostFilename);
        File FixedT0PropsFile = new File(fixedPropsFilename);
        File FullPostFile = new File(fullPosteriorFilename);
        
        ArrayList<ArrayList<Double[]>> FixedT0PropData= readTheFile(FixedT0PropsFile);
        ArrayList<ArrayList<Double[]>> FixedT0PostData= readTheFile(FixedT0PostFile);
        ArrayList<ArrayList<Double>> FullPostData= readTheFullPosteriorFile(FullPostFile);
        
        ArrayList<ArrayList<Double>> OtherDistData= readOtherDistFile(FixedT0PostFile);
        
        
        ArrayList<Double[]> FixedT0PropsPropDens = FixedT0PropData.get(1);
        ArrayList<Double[]> FixedT0PropsLogLike = FixedT0PropData.get(0);
        ArrayList<Double[]> FixedT0PostPropDens = FixedT0PostData.get(1);
        ArrayList<Double[]> FixedT0PostLogLike = FixedT0PostData.get(0);
        
        ArrayList<Double> FullPostLogLikeStar = FullPostData.get(1);
        ArrayList<Double> FullPostLogLike = FullPostData.get(0);
        ArrayList<Double> FullPostT0RefDens= FullPostData.get(2);
        ArrayList<Double> FullPostT0PriorDens = FullPostData.get(3);
        ArrayList<Double> FullPostT0StarRefDens = FullPostData.get(4);
        
        ArrayList<Double> NewPropDens = OtherDistData.get(1);
        ArrayList<Double> NewLogLike = OtherDistData.get(0);
        ArrayList<Double> PropDensRatio = OtherDistData.get(2);
        ArrayList<Double> PriorDensT0 = OtherDistData.get(3);
        Double PriorDensT0Star = OtherDistData.get(4).get(0);
        //System.out.println(PriorDensT0Star+"pi(t0*)");
        
        
        
        Integer numBridges= FixedT0PropsPropDens.get(0).length;
        
        //numerator part of the first block of the estimate:
        double FullPostEst = ChibJeliFullPostEst(FullPostLogLikeStar,FullPostLogLike,FullPostT0RefDens,FullPostT0PriorDens,FullPostT0StarRefDens,PriorDensT0Star);
        
        //second block of the estimator -- marginal likelihoods for each data point for the fixed t_0^*
        double postSize =FixedT0PostPropDens.size();
        double numYStars=0.0;//how many bridges at which we estimate the posterior density
        double[] theEsts = new double[numBridges];
        for(int i=0 ; i<postSize;i++){
        if(i%prop==0){ //get the estimate for the ith value of Y_k^* for each k
           numYStars++;
        double[] PostEst=ChibJeliIndEst(i,FixedT0PropsPropDens, FixedT0PropsLogLike, FixedT0PostPropDens, FixedT0PostLogLike, numBridges );
        for(int k=0;k<numBridges;k++){
        if(Double.isNaN(PostEst[k])){
            System.out.println("Individidual estimate for the fixed t0 ML produced NaN");
        }
        else{
        theEsts[k]+=Math.exp(FixedT0PostLogLike.get(i)[k]-Math.log(PostEst[k]));
        }
        }
        }
        }
        double theEst=0;
        for(int k=0;k<numBridges;k++){
        theEsts[k]=theEsts[k]/(numYStars);
        theEst+=Math.log(theEsts[k]);
        }
        
        //denominator part of the first block of the estimate:
        double OtherDistEst= ChibJeliOtherDistEst(FixedT0PostLogLike,NewLogLike,FixedT0PostPropDens ,NewPropDens,PropDensRatio,PriorDensT0,PriorDensT0Star);
        //System.out.println(theEst+" "+FullPostEst+" "+OtherDistEst+" "+PriorDensT0Star);
        theEst+=-FullPostEst+OtherDistEst+PriorDensT0Star;//combine all parts of the estimator
        return theEst;
    }
    
    public static double propRatio( double YPropDens,double YLogLike,double YStarPropDens, double YStarLogLike){
        //acceptance ratio of moving to YStar given Y for fixed t_0
        double alpha = YStarLogLike+YPropDens-YLogLike-YStarPropDens;
        if(Double.isNaN(alpha)){
           // System.out.println(YStarLogLike+" "+YPropDens+" "+YLogLike+" "+YStarPropDens);
        }
        
        if(YStarPropDens==Double.NEGATIVE_INFINITY){//what is this
            return(Double.NEGATIVE_INFINITY);
        }
        return(Math.min(0.0,alpha));
    }
    
    public static double propRatioFullPosterior( double LogLikeStar,double LogLike,double t0PropDens,double t0PriorDens,double t0StarPropDens, double t0StarPriorDens){
        //acceptance ratio of moving to (t_0^*,Y^*) given current value (t_0,Y)
        double alpha = LogLikeStar-LogLike+t0PropDens-t0StarPropDens+t0StarPriorDens-t0PriorDens;
        if(Double.isNaN(alpha)){
           System.out.println("Should not have NaN acceptance ratio in full posterior calcs -- investigate");
        }
       
        return(Math.min(0.0,alpha));
    }
    
    public static double propRatioOtherDist( Double[] LogLikeStar,Double LogLike, Double t0PropDensRatio,Double t0PriorDens, Double t0StarPriorDens){
        //acceptance ratio of moving from t_0^* to t_0 given fixed bridges
        //sum up the log likelihoods of all the bridges with dispersion t_0^*
        double LogLikeStarSum = 0;
        for(int i=0;i<LogLikeStar.length;i++){
            LogLikeStarSum+=LogLikeStar[i] ;
          
        }

        double alpha = LogLike-LogLikeStarSum+t0PropDensRatio-t0StarPriorDens+t0PriorDens;
        if(Double.isNaN(alpha)){
           System.out.println("Should not have NaN acceptance ratio in middle dist calcs -- investigate");
        }
       
        return(Math.min(0.0,alpha));
    }
    
    
     public static double ChibJeliFullPostEst(ArrayList<Double> LogLikeStar,ArrayList<Double> LogLike,ArrayList<Double> t0PropDens,ArrayList<Double> t0PriorDens,ArrayList<Double> t0StarPropDens ,double t0StarPriorDens){
        /*calculate the numberator part of the estimate of pi(t_0^*|x_0,x) from the full posterior sample
         */
        double numerator = 0;
        double M1 = LogLike.size();
        for(int j=0; j<M1;j++){
            numerator+=Math.exp(t0StarPropDens.get(j)+propRatioFullPosterior(LogLikeStar.get(j),LogLike.get(j),t0PropDens.get(j),t0PriorDens.get(j),t0StarPropDens.get(j),t0StarPriorDens));
   
        }
        numerator=Math.log(numerator)-Math.log(M1);
        
        
        return numerator;
    }
     
     public static double ChibJeliOtherDistEst(ArrayList<Double[]> LogLikeStar,ArrayList<Double> LogLike,ArrayList<Double[]> PropDensStar,ArrayList<Double> PropDens,ArrayList<Double> t0PropDensRatio,ArrayList<Double> t0PriorDens ,Double t0StarPriorDens){
                /*calculate the denominator part of the estimate of pi(t_0^*|x_0,x) from the middle dist sample
         */
        double numerator = 0;
        double M1 = LogLike.size();
        for(int j=0; j<M1;j++){
            numerator+=Math.exp(propRatioOtherDist(LogLikeStar.get(j),LogLike.get(j),t0PropDensRatio.get(j),t0PriorDens.get(j),t0StarPriorDens));
   
        }
        numerator=Math.log(numerator)-Math.log(M1);
        
        
        return numerator;
    }
    
    /*
     Estimate the marginal likelihood for fixed t_0^* using the ith set of bridges in the marginal posterior sample
     as Y^*
     */
    public static double[] ChibJeliIndEst(int i,ArrayList<Double[]> PropPropDen,ArrayList<Double[]> PropLogLik,ArrayList<Double[]> PostPropDen,ArrayList<Double[]> PostLogLik ,Integer numBridges){
        double[] out = new double[numBridges];//will output the individual marginal likelihood for each data point
        for(int k=0;k<numBridges;k++){
        double qYStar = PostPropDen.get(i)[k];
        double numerator = 0;
        double M1 = PostPropDen.size();
        for(int j=0; j<M1;j++){
            numerator+=Math.exp(propRatio(PostPropDen.get(j)[k],PostLogLik.get(j)[k],PostPropDen.get(i)[k],PostLogLik.get(i)[k]));
        }
        int numToIgnore=0;
        double M2 = PropPropDen.size();
        double denom = 0;
        for(int j=0; j<M2;j++){
            if(PropLogLik.get(j)[k]==Double.NEGATIVE_INFINITY){
                numToIgnore++;
            }
            else{
            denom+=Math.exp(propRatio(PostPropDen.get(i)[k],PostLogLik.get(i)[k],PropPropDen.get(j)[k],PropLogLik.get(j)[k]));
            }
        }
        
        out[k] =Math.exp(qYStar+Math.log(M2)-Math.log(M1)+Math.log(numerator)-Math.log(denom));
        }
        return out;
    }
    
    
    
    
    //read in the information for the part of the estimate calculated at fixed t_0
    public static ArrayList<ArrayList<Double[]>> readTheFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double[]>> out = new ArrayList();
        
        
        ArrayList<Double[]> propDens = new ArrayList();
        ArrayList<Double[]> logLike = new ArrayList();
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        Integer numBridges=1;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                numBridges = CountBridges(s);
                i=i+1;
               
            }
            else{
            Double[] propDensVec= new Double[numBridges];//store the proposal densities for each of the bridges
            Double[] logLikeVec= new Double[numBridges];//store the log likelihood for each of the bridges
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){
            for(int j=0; j<numBridges;j++){
            propDensVec[j]=Double.parseDouble(ss[2*j]);
            logLikeVec[j]=Double.parseDouble(ss[2*j+1]); 
            }
            propDens.add(propDensVec);
            logLike.add(logLikeVec);
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        
        return(out);
       
    }
    
    //Read the information from the 'middle' distribution -- see thesis
     public static ArrayList<ArrayList<Double>> readOtherDistFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double>> out = new ArrayList();
        
        
        ArrayList<Double> NewPropDens = new ArrayList();
        ArrayList<Double> NewLogLike = new ArrayList();
        ArrayList<Double> PropDensRatio = new ArrayList();
        ArrayList<Double> PriorDensT0 = new ArrayList();
        ArrayList<Double> PriorDensT0Star = new ArrayList();
        
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        Integer numBridges=1;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                numBridges = CountBridges(s);//get the number of bridges - middle dist info starts after ind bridge outputs
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){
            NewPropDens.add(Double.parseDouble(ss[numBridges*2]));
            NewLogLike.add(Double.parseDouble(ss[numBridges*2+1]));
            PropDensRatio.add(Double.parseDouble(ss[numBridges*2+2]));
            PriorDensT0.add(Double.parseDouble(ss[numBridges*2+3]));
            PriorDensT0Star.add(Double.parseDouble(ss[numBridges*2+4]));
            }
           
            
            }
        }
        
        out.add(NewLogLike);
        out.add(NewPropDens);
        out.add(PropDensRatio);
        out.add(PriorDensT0);
        out.add(PriorDensT0Star);
        
        return(out);
       
    }
    
    //read the information from the sample from the full posterior 
   public static ArrayList<ArrayList<Double>> readTheFullPosteriorFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double>> out = new ArrayList();
        
        
        ArrayList<Double> propDens = new ArrayList();
        ArrayList<Double> logLike = new ArrayList();
        
        ArrayList<Double> refDens = new ArrayList();
        ArrayList<Double> priorDens = new ArrayList();
        
        ArrayList<Double> t0StarRefDens = new ArrayList();
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        //introduce t so we don't end up reading last line with the acceptance rates
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){

            propDens.add(Double.parseDouble(ss[0]));
            logLike.add(Double.parseDouble(ss[1]));
            refDens.add(Double.parseDouble(ss[2]));
            priorDens.add(Double.parseDouble(ss[3]));
            t0StarRefDens.add(Double.parseDouble(ss[4]));
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        out.add(refDens);
        out.add(priorDens);
        out.add(t0StarRefDens);
        
        return(out);
       
    }
   
         /*auxilliary method to get the number of data points from the fixed t0 direct proposals file
    */ 
    private static Integer CountBridges(String header){
        int numBridges =0 ;
        boolean foundEnd=false;
        String[] ss = header.split(" ");
        while(foundEnd==false&&numBridges<ss.length){
            if(ss[numBridges].startsWith("bridge", 0)){
               numBridges +=1; 
            }
            else foundEnd=true;
        }
        
        //System.out.println(numBridges+" numBridges");
        return numBridges/2;
    }
   
}

/* //old main class to loop over many estimates to estimate Monte Carlo error
    public static void main(String[] args) throws IOException {
        
        System.out.println("Chib two block");
        //loop over the different output files to estimate the variance of the estimator
        for(int j=0;j<100;j++){
        String FixedT0PosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlock/x1_FixedDispPost"+j+".txt";
        String FixedT0PropsFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlock/x1_FixedDispPropData"+j+".txt";
        String FullPosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlock/x1_Post"+j+".txt";
        int propYStars =100;//estimate the fixed t0 part of the log marginal likelihood at every propYstars-th bridge in the posterior
        
        //calculate the estimate
        double theEst=overallEst(FullPosteriorFileName,FixedT0PosteriorFileName,FixedT0PropsFileName,propYStars);
        System.out.println(Math.exp(theEst));
        }
        

        
        
    }
   */ 