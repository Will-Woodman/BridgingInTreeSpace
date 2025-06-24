/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package simulateMultipleTopologies;
import simulation.RandomTreeSampler;
import treebase.Tree;
import topologies.TopologicalData;
import cern.jet.random.tdouble.DoubleUniform;
import diffbase.Simulator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import treebase.AlgorithmException;
import treebase.Graph;
import treebase.Graph.Edge;
import treebase.NewickStringUtils;
import treebase.RootedTree;
import treebase.TreeAsSplits;
import java.util.LinkedHashMap;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
/**
 *
 * @author will
 */
public class investigateNoOfWalkStepsLoop {

    /**
     Generate a random source tree. Then loop over a given 
     * number of steps and simulate a number of random walks started at the source tree
     Then count the topologies in the set of endpoints in the random walks for each number of steps
     * to assess convergence.
     */
    
    
    public static void main(String[] args) throws AlgorithmException {
 
        double t_0 =1;
        int[] numTaxa={25};
        int[] numSteps={10000,10,20,100,250,500,1000,2000};
        int numParticles= 10000;
        //number of the top topologies to consider
        int m = 5;
        String[] theTopologies = new String[m];
        HashMap<String, int[]> topMTopologies = new HashMap();
        for(int k=0; k< numTaxa.length; k++){
        RandomTreeSampler theTreeSampler = new RandomTreeSampler();
        Tree theTree = theTreeSampler.sampleCoalescent(numTaxa[k],2.0,1.4,false);
        System.out.println(theTree);
        

        for(int j=0; j<numSteps.length; j++){
        Tree[] treeSample = Simulator.sampleRandomWalk(theTree,t_0, numSteps[j], false, false, numParticles);
        int t=1;
        //Convert the treeSample into an ArrayList
        ArrayList<Tree> treeSampleArray = new ArrayList();
        for(Tree tree : treeSample){
            treeSampleArray.add(tree);
        }
       treeSample= null;
       //get the topology data from the sample
       TopologicalData treeTopologyData = new TopologicalData(treeSampleArray);
       //The number of topologies in the data set
       System.out.print(numTaxa[k]);
       System.out.print(" ");
       System.out.print(numSteps[j]);
       System.out.print(" ");
       System.out.print(" ");
       System.out.print(treeTopologyData.getTheCounts().size());
       System.out.print(" ");
       System.out.print(treeTopologyData.getTotalCount());
       System.out.print(" ");
       System.out.println(treeTopologyData.getTheSplitCounts().size());
       
       
       //Order the topology data in terms of how many times the topology appears
       HashMap<String, Integer> theCounts = treeTopologyData.getTheCounts();
       // Now sort the hashmap
        List<Map.Entry<String, Integer>> list = new LinkedList<Map.Entry<String, Integer>>(theCounts.entrySet());
        Collections.sort(list, new MyComparator2());
        LinkedHashMap<String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
        for (Map.Entry<String, Integer> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());   
         }
     
        if(j==0){
               Iterator theIterator = sortedMap.keySet().iterator();
               int z=0;
            for(int l=0 ; l<m ;l++){

                theTopologies[l]= (String) theIterator.next();
                int[] currentInt= new int[numSteps.length];
                currentInt[j]=sortedMap.get(theTopologies[l]);
                topMTopologies.put(theTopologies[l],currentInt);
            }
        }
        else{
            for(String top : theTopologies){
                int[] currentInts = topMTopologies.get(top);
                if(theCounts.get(top)!=null){
                currentInts[j]=theCounts.get(top);
                topMTopologies.replace(top, currentInts);
            }
                else{
                        currentInts[j]=0;
                topMTopologies.replace(top, currentInts);
                        }
            }
        }
        }
        for(int l=0; l<m ; l++){
            System.out.print(theTopologies[l]);
            for(int j=0 ; j<numSteps.length;j++){
            System.out.print(topMTopologies.get(theTopologies[l])[j]);
            System.out.print(" ");
            }
            System.out.println(" ");
       
        }
        }

    
        }
   
 private static class MyComparator2 implements Comparator<Map.Entry<String, Integer>> {
        public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
            return -o1.getValue().compareTo(o2.getValue()); // Want big integers first!
        }
    }
    
}
