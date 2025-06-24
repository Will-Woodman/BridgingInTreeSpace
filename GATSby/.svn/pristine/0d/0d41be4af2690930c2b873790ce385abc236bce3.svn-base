/*
 * RandomTreeSampler
    Copyright (C) 2013  Tom M. W. Nye

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

package simulation;

import cern.jet.random.tdouble.DoubleUniform;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import treebase.AlgorithmException;
import treebase.Graph;
import treebase.Graph.Edge;
import treebase.NewickStringUtils;
import treebase.RootedTree;
import treebase.Tree;

/**
 * Random tree generation. Can generate:
 * (i)  a rooted tree using the Yule model of speciation OR
 * (ii) an unrooted tree by randomly resolving an appropriate star tree, or a 
 *      rooted tree by rooting a randomly resolved star tree on a randomly
 *      selected edge
 */

public class RandomTreeSampler implements java.io.Serializable, RandomEngineSeedSetter  {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /* Data */
    private GammaDistribution gamDistrib;
    private DoubleUniform u;
    
    /** Constructor */
    public RandomTreeSampler() {
        gamDistrib = new GammaDistribution(1.0, 1.0);
        u = new DoubleUniform(Random.getEngine());
    }
    
    @Override
    public void resetRandomEngineSeed() {
        gamDistrib.resetRandomEngineSeed();
        u = new DoubleUniform(Random.getEngine());
    }
    
    /** YULE MODEL OF SPECIATION -------------------------------------------- */
    
    /** Samples a N-taxa tree from the Yule distribution with the given birth 
        rate and taxa names t1, t2, ..., tN */
    public RootedTree sampleYule(double birthRate, int numTaxa) 
            throws AlgorithmException {
        
        HashSet<String> taxaNames = new HashSet<String>();
        for(int i=0; i<numTaxa; i++) taxaNames.add("t"+(i+1));
        return sampleYule(birthRate, taxaNames);
        
    }
    
    /** As above but resamples the edge lengths from a Gamma(shape, scale) 
        distribution*/
    public RootedTree sampleYule(int numTaxa, double shape, 
            double scale) throws AlgorithmException {
        
        RootedTree tree = sampleYule(1.0, numTaxa);
        HashSet<Graph.Edge> edges = tree.getEdges();
        for(Graph.Edge e : edges) e.setLength(gamDistrib.sample(shape, scale));
        return tree;
        
    }
    
    /** Samples a tree from the Yule distribution with the given birth rate and 
        taxa names*/
    public RootedTree sampleYule(double birthRate, HashSet<String> taxaNames) 
            throws AlgorithmException {
        
        int numSpecies = taxaNames.size();
        // Create tree with one taxon
        RootedTree tree = new RootedTree(taxaNames.iterator().next()); 
                                              // Trick needed to ensure taxa
                                              // hashset is correct at end
        tree.getRoot().label = "";
        // Store leaves and external edges
        HashSet<Graph.Vertex> extant = new HashSet<Graph.Vertex>();
        extant.add(tree.getRoot());
        HashSet<Graph.Edge> externalEdges = new HashSet<Graph.Edge>();
        int N;
        double dt;
        do {
            N = extant.size();
            // Sample interevent time
            dt = gamDistrib.sample(1.0, 1.0/(birthRate*N));
            // Sample extant species uniformly at random
            Graph.Vertex v = (Graph.Vertex) CategoricalDistribution.
                    sampleFromSetWithReplacement(extant, u);
            // Get ancestor
            Graph.Vertex w0 = null;
            if(!tree.getRoot().equals(v)) w0 = v.neighboursIterator().next();
            // Split sampled node
            Graph.Vertex w1 = tree.addNewVertex("");
            Graph.Vertex w2 = tree.addNewVertex("");
            Graph.Edge e1 = tree.connect(v, w1, 0.0);
            Graph.Edge e2 = tree.connect(v, w2, 0.0);
            // Update set of external edges and increment associated edge lengths
            externalEdges.add(e1);
            externalEdges.add(e2);
            if(w0 != null) externalEdges.remove(v.getEdge(w0));
            for(Graph.Edge e : externalEdges) {
                e.setLength(e.getLength() + dt);
            }
            // Update set of extant species
            extant.remove(v);
            extant.add(w1);
            extant.add(w2);
        } while(N<(numSpecies-1));
        // Label leaves and update hashset of taxa
        HashSet<String> taxa = new HashSet<String>();
        taxa.addAll(taxaNames);
        for(Graph.Vertex v : extant) {
            String label = (String) CategoricalDistribution.
                    sampleFromSetWithoutReplacement(taxa, u);
            Graph.Vertex w = v.neighboursIterator().next();
            double length = w.getNeighbourDistance(v);
            tree.remove(v);
            Graph.Vertex vnew = tree.addNewVertex(label);
            tree.connect(vnew, w, length);
        }
        // Calculate ancestry
        tree.recalculateAncestry();
        return tree;
    }
    
    /** As above but resamples the edge lengths from a Gamma(shape, scale) 
        distribution*/
    public RootedTree sampleYule(HashSet<String> taxaNames, 
            double shape, double scale) throws AlgorithmException {
        
        RootedTree tree = sampleYule(1.0, taxaNames);
        HashSet<Graph.Edge> edges = tree.getEdges();
        for(Graph.Edge e : edges) e.setLength(gamDistrib.sample(shape, scale));
        return tree;
        
    }
    
    /** RANDOMLY RESOLVING STAR TREE ---------------------------------------- */
    
    /** Given the specified number of taxa, samples a tree with unit branch lengths. 
        The default is to sample an unrooted tree. The taxa are labelled t1,t2,... */
    public Tree sampleCoalescent(int numTaxa) throws AlgorithmException {
        return sampleCoalescent(numTaxa, false);
    }
    public Tree sampleCoalescent(int numTaxa, boolean rootedTree) throws AlgorithmException {
        if(numTaxa<=0) throw new AlgorithmException("Random tree must have at least 1 taxon.");
        HashSet<String> taxonNames = new HashSet<String>();
        for(int i=0; i<numTaxa; i++) {
            taxonNames.add("t"+(i+1));
        }
        if(numTaxa==1) {
            if(rootedTree) throw new AlgorithmException("Random rooted tree must have at least 2 taxa.");
            else return sampleCoalescent(null, taxonNames, rootedTree);
        }
        double[] edgeLengths;
        if(rootedTree) edgeLengths = new double[2*numTaxa-2];
        else edgeLengths = new double[2*numTaxa-3];
        for(int i=0; i<edgeLengths.length; i++) {
            edgeLengths[i] = 1.0;
        }
        return sampleCoalescent(edgeLengths, taxonNames, rootedTree);
    }
    
    /** As above but the branch lengths are sampled from a Gamma(shape, scale) 
        distribution. */
    public Tree sampleCoalescent(int numTaxa, double shape, double scale) throws AlgorithmException {
        return sampleCoalescent(numTaxa, shape, scale, false);
    }
    public Tree sampleCoalescent(int numTaxa, double shape, double scale, boolean rootedTree) throws AlgorithmException {
        if(numTaxa <=0) throw new AlgorithmException("Random tree must have at least 1 taxon.");
        HashSet<String> taxonNames = new HashSet<String>();
        for(int i=0; i<numTaxa; i++) {
            taxonNames.add("t"+(i+1));
        }
        if(numTaxa==1) {
            if(rootedTree) throw new AlgorithmException("Random rooted tree must have at least 2 taxa.");
            else return sampleCoalescent(null, taxonNames, rootedTree);
        }
        double[] edgeLengths;
        if(rootedTree) edgeLengths = new double[2*numTaxa-2];
        else edgeLengths = new double[2*numTaxa-3];
        for(int i=0; i<edgeLengths.length; i++) {
            edgeLengths[i] = gamDistrib.sample(shape, scale);
        }
        return sampleCoalescent(edgeLengths, taxonNames, rootedTree);
    }
    
    /** Given the specified set of taxon names, samples a tree with unit branch lengths. 
        The default is to sample an unrooted tree. */
    public Tree sampleCoalescent(HashSet<String> taxonNames) throws AlgorithmException {
        return sampleCoalescent(taxonNames, false);
    }
    public Tree sampleCoalescent(HashSet<String> taxonNames, boolean rootedTree) throws AlgorithmException {
        if(taxonNames==null) throw new AlgorithmException("Random tree must have at least 1 taxon.");
        int numTaxa = taxonNames.size();
        if(numTaxa==1) {
            if(rootedTree) throw new AlgorithmException("Random rooted tree must have at least 2 taxa.");
            else return sampleCoalescent(null, taxonNames, rootedTree);
        }
        double[] edgeLengths;
        if(rootedTree) edgeLengths = new double[2*numTaxa-2];
        else edgeLengths = new double[2*numTaxa-3];
        for(int i=0; i<edgeLengths.length; i++) {
            edgeLengths[i] = 1.0;
        }
        return sampleCoalescent(edgeLengths, taxonNames, rootedTree);
    }
    
    /** As above but the branch lengths are sampled from a Gamma(shape, scale) 
        distribution. */
    public Tree sampleCoalescent(HashSet<String> taxonNames, double shape, double scale) throws AlgorithmException {
        return sampleCoalescent(taxonNames, shape, scale, false);
    }
    public Tree sampleCoalescent(HashSet<String> taxonNames, double shape, double scale, boolean rootedTree) throws AlgorithmException {
        if(taxonNames==null) throw new AlgorithmException("Random tree must have at least 1 taxon.");
        int numTaxa = taxonNames.size();
        if(numTaxa==1) {
            if(rootedTree) throw new AlgorithmException("Random rooted tree must have at least 2 taxa.");
            else return sampleCoalescent(null, taxonNames, rootedTree);
        }
        double[] edgeLengths;
        if(rootedTree) edgeLengths = new double[2*numTaxa-2];
        else edgeLengths = new double[2*numTaxa-3];
        for(int i=0; i<edgeLengths.length; i++) {
            edgeLengths[i] = gamDistrib.sample(shape, scale);
        }
        return sampleCoalescent(edgeLengths, taxonNames, rootedTree);
    }
    
    /** Samples a rooted or unrooted tree with the specified branch lengths
        and taxon names. */
    private Tree sampleCoalescent(double[] edgeLengths, HashSet<String> taxonNames, boolean rootedTree) throws AlgorithmException {
        int numTaxa = taxonNames.size();
        Iterator<String> itS = taxonNames.iterator();
        if(numTaxa==1) {
            if(rootedTree) throw new AlgorithmException("Random rooted tree must have at least 2 taxa. This "
                    + "error should not be possible.");
            else return new Tree(itS.next()+";");
        }
        else if(numTaxa==2) {
            String name1 = itS.next();
            String name2 = itS.next();
            if(rootedTree) return new treebase.RootedTree("("+name1+":"+(edgeLengths[0])+","+name2+":"+(edgeLengths[1])+");");
            else {
                Tree theTree = new Tree("("+name1+":"+(edgeLengths[0]/2)+","+name2+":"+(edgeLengths[0]/2)+");");
                theTree.removeDegreeTwoVertices();
                return theTree;
            }
        } else if(numTaxa==3) {
            String name1 = itS.next();
            String name2 = itS.next();
            String name3 = itS.next();
            Tree theTree = new Tree("("+name1+":"+edgeLengths[0]+","+name2+":"+edgeLengths[1]+","+name3+":"+edgeLengths[2]+");");
            if(rootedTree) {
                Graph.Edge e = (Graph.Edge) CategoricalDistribution.sampleFromSetWithReplacement(theTree.getEdges(), u);
                return new treebase.RootedTree(theTree, e);
            }
            else return theTree;
        } else {
            // Make a star tree
            String s = "(";
            for (int i = 0; i < numTaxa; i++) {
                String name = itS.next();
                s += name + ":1.0,";
            }
            s = s.substring(0, s.length() - 1);
            s += ");";
            // Resolve it
            Tree t = randomlyResolve(s, rootedTree);
            // Set new edge lengths
            Iterator<Edge> itE = t.getEdgeIterator();
            for (int i = 0; i < t.numEdges(); i++) {
                itE.next().setLength(edgeLengths[i]);
            }
            return t;
        }
    }
    
    /** Randomly reroots a rooted tree */
    public RootedTree randomlyReroot(RootedTree tree) throws AlgorithmException {
        if(tree.getRoot().degree()!=2) throw new AlgorithmException("Rerooting only possible for tree with degree 2 root.");
        Tree t = Tree.copy(tree);
        // Find old root
        Iterator<Graph.Vertex> itV = t.getVertexIterator();
        Graph.Vertex oldRoot = null;
        while(itV.hasNext()) {
            Graph.Vertex v = itV.next();
            if(v.degree()==2) {
                oldRoot = v;
                break;
            }
        }
        // Choose new edge for rooting
        HashSet<Edge> edges = t.getEdges();  
        Graph.Edge newRootEdge = (Graph.Edge)simulation.CategoricalDistribution.sampleFromSetWithReplacement(edges, u);
        // Remove old root
        Iterator<Graph.Vertex> itN = oldRoot.neighboursIterator();
        Graph.Vertex vA = (Graph.Vertex)itN.next(); Graph.Vertex vB = (Graph.Vertex)itN.next();
        double mergedLength = vA.getEdge(oldRoot).getLength() + vB.getEdge(oldRoot).getLength();
        t.remove(oldRoot);
        t.connect(vA, vB, mergedLength);
        // Subdivide the new root edge and label the new vertex as the root
        Graph.Vertex newRoot = t.divide(newRootEdge, "");
        return new RootedTree(t, newRoot);
    }
    

    /** Build a tree from the Newick string newickStr, then resolve it by randomly 
        resolving degree > 3 vertices. Return the so-created resolved tree.
     */
     private Tree randomlyResolve(String newickStr, boolean rootedTree) {

         // Build tree
         Tree theTree;
         try {
             theTree = new Tree(newickStr);
         } catch (AlgorithmException anErr) {
             System.out.println("Error building star tree. "+anErr.getMessage());
             theTree = null;
         }

         // Find a degree>3 vertex
         Graph.Vertex v;
         do {
             v = null;
             ArrayList<Graph.Vertex> vertices = theTree.getOrderedListOfVertices();
             for (int i=0; i<theTree.numVertices(); i++) {
                 Graph.Vertex w = vertices.get(i);
                 if (w.degree()>3) v=w;
             }
             if (v!=null) {
                 // OK, found a vertex to resolve
                 try {
                     // Resolve this vertex.

                     // Start by adding two new vertices
                     Graph.Vertex vA = theTree.addNewVertex("");
                     Graph.Vertex vB = theTree.addNewVertex("");
                     theTree.connect(vA, vB, 1.0); // New branch

                     // Get two random neighbours & join to vA:
                     // Start by getting an ordered list of neighbours
                     ArrayList<Graph.Vertex> orderedNeighbours = new ArrayList();
                     for (int i=0; i<vertices.size(); i++) {
                         Graph.Vertex w = vertices.get(i);
                         if (w.isConnectedTo(v)) orderedNeighbours.add(w);
                     }
                     // Sanity check
                     if (orderedNeighbours.size()!=v.degree()) System.out.println("Warning: algorithm for resolving a vertex randomly has failed. This should not be possible.");

                     // Get two integers from 0 to deg(v)-1
                     int[] k = new int[2];
                     k[0] = u.nextIntFromTo(0, (v.degree()-1));
                     k[1] = u.nextIntFromTo(0, (v.degree()-2));
                     k[1] = (k[1]<k[0]) ? k[1] : k[1]+1;

                     // Get two random neighbours & join to vA
                     for (int i=0; i<2; i++) {
                         Graph.Vertex w = orderedNeighbours.get(k[i]);
                         Graph.Edge e  = v.getEdge(w);
                         double z = e.getLength();
                         theTree.cut(e);
                         Graph.Edge f = theTree.connect(vA, w, z);
                     }
                     // Join remaining two vertices to vB
                     for (int i=0; i<orderedNeighbours.size(); i++) {
                         if ((i!=k[0])&&(i!=k[1])) {
                             // This vertex should be joined to vB
                             Graph.Vertex w = orderedNeighbours.get(i);
                             Graph.Edge e  = v.getEdge(w);
                             double z = e.getLength();
                             theTree.cut(e);
                             Graph.Edge f = theTree.connect(vB, w, z);
                         }
                     }
                     if ((vA.degree()<3)||(vB.degree()<3)) {
                         System.out.println("Warning: when resolving a tree new vertex has degree less than 3. Trouble ahead...");
                     }
                     theTree.remove(v);

                 }
                 catch (AlgorithmException anErr) {
                     System.out.println("Error during tree resolution. "+anErr.getMessage());
                 }

             }
         }
         while (v!=null); // If you found an unresolved vertex, then loop thru' again

         if(rootedTree) {
             // Randomly choose edge on which to root
             Graph.Edge e = (Graph.Edge) CategoricalDistribution.sampleFromSetWithReplacement(theTree.getEdges(), u);
             try {
                 theTree = new treebase.RootedTree(theTree, e);
             } catch (AlgorithmException anErr) {
                 System.out.println("Error during tree rooting. "+anErr.getMessage());
                 theTree = null;
             }
             
         }

         return theTree;
     }     
     
     /* ---------------------------------------------------------------------- */

    
}
