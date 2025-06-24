/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package simulateMultipleTopologies;
import MCMC.PosteriorAnalysis;
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
public class investigateNoOfWalkSteps {

    /**
     Generate a random source tree and a number of random walks started at the source tree
     * Then count the topologies in the set of endpoints in the random walks
     */
    
    
    public static void main(String[] args) throws AlgorithmException {

        double t_0 =1;
        int numTaxa=10;
        int numSteps=500;
        int numParticles= 10000;
        
        //simulate the source tree
        RandomTreeSampler theTreeSampler = new RandomTreeSampler();
        Tree theTree = theTreeSampler.sampleCoalescent(numTaxa,2.0,0.05,false);
        //run the random walks
        Tree[] treeSample = Simulator.sampleRandomWalk(theTree,t_0, numSteps, false, true, numParticles);

        //Convert the treeSample into an ArrayList
        ArrayList<Tree> treeSampleArray = new ArrayList();
        for(Tree tree : treeSample){
            treeSampleArray.add(tree);
        }
       //get the topology data from the sample
       TopologicalData treeTopologyData = new TopologicalData(treeSampleArray);
       //Output the number of topologies in the data set
       System.out.println(treeTopologyData.getTheCounts().size());
       System.out.println(treeTopologyData.getTotalCount());
       System.out.println(treeTopologyData.getTheSplitCounts().size());
       //Order the topology data in terms of how many times the topology appears
       HashMap<String, Integer> theCounts = treeTopologyData.getTheCounts();
       
       /* Now sort the hashmap */
        List<Map.Entry<String, Integer>> list = new LinkedList<Map.Entry<String, Integer>>(theCounts.entrySet());
        Collections.sort(list, new MyComparator2());
        LinkedHashMap<String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
        for (Map.Entry<String, Integer> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        
}

    private static class MyComparator2 implements Comparator<Map.Entry<String, Integer>> {
        public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
            return -o1.getValue().compareTo(o2.getValue()); // Want big integers first!
        }
    }
    
}
