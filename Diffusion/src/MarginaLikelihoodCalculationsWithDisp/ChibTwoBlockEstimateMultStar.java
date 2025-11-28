/*
ChibTwoBlockEstimateMultStar
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
public class ChibTwoBlockEstimateMultStar {

     /**
     Main class to calculate the Chib two block estimate of the marginal likelihood after running the 
     * sampler, where we use multiple values of the fixed t_0^* rather than just one
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Chib");
        int numT0Stars =4;
        int numT0Props =2;
        int propYStars =100;//estimate the fixed t0 part of the log marginal likelihood at every propYstars-th bridge in the marginal posterior
        
        //loop over the different output files to estimate the variance of the estimator
        for(int j=0;j<100;j++){
        String FixedT0PosteriorFileNamept1="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_FixedDisp";
        String FixedT0PosteriorFileNamept2 ="PostOut"+j+".txt";
        String FixedT0PropsFileNamept1="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_FixedDisp";
        String FixedT0PropsFileNamept2="PropData"+j+".txt";
        String FullPosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaWith20250212/ChibTwoBlockStar/x1_PostOut"+j+".txt";
        
        double theEst = overallEst(FullPosteriorFileName,FixedT0PosteriorFileNamept1,FixedT0PosteriorFileNamept2,FixedT0PropsFileNamept1,FixedT0PropsFileNamept2,numT0Stars,numT0Props,propYStars);
        System.out.println(Math.exp(theEst));
        }
        
        
        
    }
    
        private static double overallEst(String fullPosteriorFilename,String fixedPostFilenamept1,String fixedPostFilenamept2,String fixedPropsFilenamept1,String fixedPropsFilenamept2,int numStars,int numProps,int prop) throws IOException{
        double[] FixedT0PostEsts = new double[numStars];
        double[] OtherDistEst = new double[numStars];
        double[] t0StarPrior = new double[numStars];
        double[] FixedT0Ests = new double[3];
  
        for(int k=1; k<=numStars;k++){
        String FixedT0PosteriorFileName=fixedPostFilenamept1+k+fixedPostFilenamept2;
        String FixedT0PropsFileName=fixedPropsFilenamept1+k+fixedPropsFilenamept2;
        FixedT0Ests = FixedT0Est(FixedT0PosteriorFileName, FixedT0PropsFileName,prop);
        FixedT0PostEsts[k-1] = FixedT0Ests[0];//ML for fixed t0
        OtherDistEst[k-1] = FixedT0Ests[1];//denominator part of the first block of the estimate
        t0StarPrior[k-1] = FixedT0Ests[2];
        }
        
        //numerator part of the first block of the estimate
        File FullPostFile = new File(fullPosteriorFilename);
        ArrayList<ArrayList<Double[]>> FullPostData= readTheFullPosteriorFile(FullPostFile);

        ArrayList<Double[]> FullPostLogLikeStar = FullPostData.get(1);
        ArrayList<Double[]> FullPostLogLike = FullPostData.get(0);
        ArrayList<Double[]> FullPostT0RefDens= FullPostData.get(2);
        ArrayList<Double[]> FullPostT0PriorDens = FullPostData.get(3);
        ArrayList<Double[]> FullPostT0StarRefDens = FullPostData.get(4); 
        
        double fullEst=0;
        double[] FullPostEst = ChibJeliFullPostEst(FullPostLogLikeStar,FullPostLogLike,FullPostT0RefDens,FullPostT0PriorDens,FullPostT0StarRefDens,t0StarPrior);
        //combine all the estimates
        for(int k=0; k<numStars;k++){
        FixedT0PostEsts[k]+=-FullPostEst[k]+OtherDistEst[k]+t0StarPrior[k];

        //System.out.println(Math.exp(FixedT0PostEsts[k]));
        fullEst+=FixedT0PostEsts[k];
        }
        return fullEst/numStars;
        }
    
    /*estimate marginal lieklihood over bridges for a fixed value of t0Star
      and the numerator part of the first block of the estimate from the middle distribution
        */
    public static double[] FixedT0Est( String PostFileName,String PropsFileName,int prop) throws IOException{
       //read in the files:
        File PostFile = new File(PostFileName);
        File PropsFile = new File(PropsFileName);
        
        ArrayList<ArrayList<Double[]>> FixedT0PropData= readTheFile(PropsFile);
        ArrayList<ArrayList<Double[]>> FixedT0PostData= readTheFile(PostFile);
        ArrayList<ArrayList<Double[]>> OtherDistData= readOtherDistFile(PostFile);
        
        ArrayList<Double[]> PropsPropDens = FixedT0PropData.get(1);
        ArrayList<Double[]> PropsLogLike = FixedT0PropData.get(0);
        ArrayList<Double[]> PostPropDens = FixedT0PostData.get(1);
        ArrayList<Double[]> PostLogLike = FixedT0PostData.get(0);
        
        ArrayList<Double[]> NewLogLike = OtherDistData.get(0);
        ArrayList<Double[]> PropDensRatio = OtherDistData.get(1);
        ArrayList<Double[]> PriorDensT0 = OtherDistData.get(2);
        Double PriorDensT0Star = OtherDistData.get(3).get(0)[0];
        
        Integer numBridges= PropsPropDens.get(0).length;//or equivalently the number of data points
        double postSize =PostPropDens.size();
        double numYStars=0.0;//how many bridges at which we estimate the posterior density
        double[] theEsts = new double[numBridges];
        int[] numNans =new int[numBridges];
        for(int i=0 ; i<postSize;i++){
        if(i%prop==0){
            numYStars++;
        double[] PostEst=ChibJeliIndEst(i,PropsPropDens, PropsLogLike, PostPropDens, PostLogLike, numBridges );
        for(int k=0;k<numBridges;k++){
        if(Double.isNaN(PostEst[k])){
            System.out.println("Individidual estimate for the fixed t0 ML produced NaN");
        }
        else{
        theEsts[k]+=Math.exp(PostLogLike.get(i)[k]-Math.log(PostEst[k]));
        }
        }
        }
        }
        double theEst=0;
        for(int k=0;k<numBridges;k++){
        theEsts[k]=theEsts[k]/(numYStars);
        theEst+=Math.log(theEsts[k]);
        
    }
        double OtherDistEst= ChibJeliOtherDistEst(PostLogLike,NewLogLike,PropDensRatio,PriorDensT0,PriorDensT0Star);
        double[] Output = {theEst,OtherDistEst,PriorDensT0Star};
        return Output;
    }
    
    public static double propRatio( double YPropDens,double YLogLike,double YStarPropDens, double YStarLogLike){
        //acceptance ratio of moving to YStar given Y for fixed t_0
        double alpha = YStarLogLike+YPropDens-YLogLike-YStarPropDens;
        if(Double.isNaN(alpha)){
           // System.out.println(YStarLogLike+" "+YPropDens+" "+YLogLike+" "+YStarPropDens);
        }
        
        if(YStarPropDens==Double.NEGATIVE_INFINITY){
            return(Double.NEGATIVE_INFINITY);
        }
        return(Math.min(0.0,alpha));
    }
    
    public static double propRatioFullPosterior( double LogLikeStar,double LogLike,double t0PropDens,double t0PriorDens,double t0StarPropDens, double t0StarPriorDens){
        //acceptance ratio of moving to (t_0^*,Y) given current value (t_0,Y)
        double alpha = LogLikeStar-LogLike+t0PropDens-t0StarPropDens+t0StarPriorDens-t0PriorDens;
        if(Double.isNaN(alpha)){
           System.out.println("Should not have NaN acceptance ratio in full posterior calcs -- investigate");
        }
       
        return(Math.min(0.0,alpha));
    }
    
    public static double propRatioOtherDist( Double[] LogLikeStar,Double LogLike,Double t0PropDensRatio,Double t0PriorDens, Double t0StarPriorDens){
        //acceptance ratio of moving from t_0^* to t_0 given fixed bridges
        //sum up the log likelihoods of all the bridges with dispersion t_0^*
        double LogLikeStarSum = 0;
        for(int i=0;i<LogLikeStar.length;i++){
            LogLikeStarSum+=LogLikeStar[i] ; 
        }
        

        double alpha = LogLike-LogLikeStarSum-t0PropDensRatio-t0StarPriorDens+t0PriorDens;
        if(Double.isNaN(alpha)){
           System.out.println("Should not have NaN acceptance ratio in middle dist calcs -- investigate");
        }
       
        return(Math.min(0.0,alpha));
    }
    
    
     public static double[] ChibJeliFullPostEst(ArrayList<Double[]> LogLikeStar,ArrayList<Double[]> LogLike,ArrayList<Double[]> t0PropDens,ArrayList<Double[]> t0PriorDens,ArrayList<Double[]> t0StarPropDens ,double[] t0StarPriorDens){
        /*calculate the numberator part of the estimate of pi(t_0^*|x_0,x) from the full posterior sample for each
         fixed value of t_0^*
         */
        int numStars = LogLikeStar.get(0).length;
        double[] numerator = new double[numStars];
        double M1 = LogLike.size();
        for(int j=0; j<M1;j++){
            for(int k=0; k<numStars;k++){
            numerator[k]+=Math.exp(t0StarPropDens.get(j)[k]+propRatioFullPosterior(LogLikeStar.get(j)[k],LogLike.get(j)[0],t0PropDens.get(j)[k],t0PriorDens.get(j)[0],t0StarPropDens.get(j)[k],t0StarPriorDens[k]));
            }
        }
        for(int k=0; k<numStars;k++){
        numerator[k]=Math.log(numerator[k])-Math.log(M1);
        }
        
        return numerator;
    }
     
     public static double ChibJeliOtherDistEst(ArrayList<Double[]> LogLikeStar,ArrayList<Double[]> LogLike,ArrayList<Double[]> t0PropDensRatio,ArrayList<Double[]> t0PriorDens ,Double t0StarPriorDens){
         /*calculate the denominator part of the estimate of pi(t_0^*|x_0,x) from the middle dist sample
         */
        int numProps = LogLike.get(0).length;
        double numerator = 0;
        double M1 = LogLike.size();
        for(int j=0; j<M1;j++){
            for(int k=0;k<numProps;k++){
            numerator+=Math.exp(propRatioOtherDist(LogLikeStar.get(j),LogLike.get(j)[k],t0PropDensRatio.get(j)[k],t0PriorDens.get(j)[k],t0StarPriorDens));
            }
        }
        numerator=Math.log(numerator)-Math.log(M1)-Math.log(numProps);
        
        
        return numerator;
    }
     /*
     Estimate the marginal likelihood for a fixed t_0^* using the ith set of bridges in the marginal posterior sample
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
        
        double M2 = PropPropDen.size();
        double denom = 0;
        for(int j=0; j<M2;j++){
            if(PropLogLik.get(j)[k]==Double.NEGATIVE_INFINITY){
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
            Double[] propDensVec= new Double[numBridges];
            Double[] logLikeVec= new Double[numBridges];
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
     public static ArrayList<ArrayList<Double[]>> readOtherDistFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double[]>> out = new ArrayList();
        
        
        ArrayList<Double[]> NewLogLike = new ArrayList();
        ArrayList<Double[]> PropDensRatio = new ArrayList();
        ArrayList<Double[]> PriorDensT0 = new ArrayList();
        ArrayList<Double[]> PriorDensT0Star = new ArrayList();
        
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        Integer numBridges=1;
        Integer numProps=1;
        Double[] PriorDensStar = new Double[1];
        //introduce t so we don't end up reading last line with the acceptance rates
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                numBridges = CountBridges(s);
                numProps = CountProps(s);
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            Double[] LogLike = new Double[numProps];
            Double[] PropRatio = new Double[numProps];
            Double[] PriorDens = new Double[numProps];
            if(!(ss[1].equals("rate"))){
            for(int j=0;j<numProps;j++){
                LogLike[j]=Double.parseDouble(ss[numBridges*2+j*3]);
                PropRatio[j]=Double.parseDouble(ss[numBridges*2+j*3+1]);
                PriorDens[j]=Double.parseDouble(ss[numBridges*2+j*3+2]);
            }
            PriorDensStar[0]=Double.parseDouble(ss[numBridges*2+numProps*3]);
            
            NewLogLike.add(LogLike);
            PropDensRatio.add(PropRatio);
            PriorDensT0.add(PriorDens);
            PriorDensT0Star.add(PriorDensStar);
            }
           
            
            }
        }
        
        out.add(NewLogLike);
        out.add(PropDensRatio);
        out.add(PriorDensT0);
        out.add(PriorDensT0Star);
        
        return(out);
       
    }
    
       //read the information from the sample from the full posterior   
   public static ArrayList<ArrayList<Double[]>> readTheFullPosteriorFile(File theFile) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double[]>> out = new ArrayList();
        int numStars= 1;
        
        ArrayList<Double[]> logLikeStar = new ArrayList();
        ArrayList<Double[]> logLike = new ArrayList();
        
        ArrayList<Double[]> refDens = new ArrayList();
        ArrayList<Double[]> priorDens = new ArrayList();
        
        ArrayList<Double[]> t0StarRefDens = new ArrayList();
        
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        //introduce t so we don't end up reading last line with the acceptance rates
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                i=i+1;
               numStars = countStars(s);
            }
            else{
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){
              Double[] likeStar = new Double[numStars];
              Double[] like = new Double[1];
              Double[] prior = new Double[1];
              Double[] propDen = new Double[numStars];//t0 proposal density from each fixed value to the current value of t_0
              Double[] starPropDen = new Double[numStars];//t0 proposal density from the current value of t_0 to each fixed value 
              like[0]=Double.parseDouble(ss[numStars*3]);
              prior[0]=Double.parseDouble(ss[numStars*3+1]);
              
              for(int k=0;k<numStars;k++){
                  likeStar[k] = Double.parseDouble(ss[k*3]);
                  propDen[k] = Double.parseDouble(ss[k*3+1]);
                  starPropDen[k] = Double.parseDouble(ss[k*3+2]);
              }
            logLikeStar.add(likeStar);
            logLike.add(like);
            refDens.add(propDen);
            priorDens.add(prior);
            t0StarRefDens.add(starPropDen);
             
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(logLikeStar);
        out.add(refDens);
        out.add(priorDens);
        out.add(t0StarRefDens);
        
        return(out);
       
    }
   
            /*auxilliary method to get the number of data points from fixed t0 the direct proposals file
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
        

        return numBridges/2;
    }
    
    /*function to count how many times t0 was proposed from the pi*rho distribution for each set of bridges
    */
    private static Integer CountProps(String header){
        int numProps =0 ;
        int j=0;
        boolean foundEnd=false;
        String[] ss = header.split(" ");
        while(foundEnd==false&&numProps<ss.length){
            if(ss[j].startsWith("TotalLogLikeNewDisp", 0)||ss[j].startsWith("DispPropRatio", 0)||ss[j].startsWith("T0Prior", 0)){
               numProps +=1; 
               j+=1;
            }
            else if(ss[j].startsWith("bridge", 0)){
                j+=1;
            }
            else foundEnd=true;
        }

        return numProps/3;
    }
    
   /*function to count how many distinct values of the fixed t_0^* were used in the simulations
    */
    private static Integer countStars(String header){
        int numStars =0 ;
        boolean foundEnd=false;
        String[] ss = header.split(" ");
        while(foundEnd==false&&numStars<ss.length){
            if(ss[numStars].startsWith("LogLikeStar", 0)||ss[numStars].startsWith("t0PropDens", 0)||ss[numStars].startsWith("t0StarPropDens", 0)){
               numStars +=1; 
            }
            else foundEnd=true;
        }
        
        return numStars/3;
    }
   
}
