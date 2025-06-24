/**
    Split
    Class representing a bipartition of a set of taxa.

    Copyright (C) 2011  Tom M. W. Nye

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

package treebase;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.ListIterator;

/**
 * Class representing a bipartition of a set of taxa.
 * The taxa are referred to by Strings (names)
 */


public class Split implements Comparable, java.io.Serializable {

    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    // Store two hash sets of leaves
    private HashSet<String> setA;
    private HashSet<String> setB;

    // Hashcode
    private final int hc;
    /* In future if you write any methods which change the nature of the split, you will need
     to recompute the hash code cf. doing a topological move on a tree means you need to recompute splits. */

    private String eve;

     /** Constructor based on two HashSets */
    public Split(HashSet<String> a, HashSet<String> b) throws AlgorithmError {
        HashSet x = intersection(a, b);
        if (x.size()>0) throw new AlgorithmError("Attempt to create a Split with intersecting taxa.");
        setA = new HashSet<String>();
        setB = new HashSet<String>();

        String s = Collections.min(a);
        String t = Collections.min(b);
        if (s.compareTo(t)<0) {
            eve = s;
            setA.addAll(a);
            setB.addAll(b);
        }
        else {
            eve = t;
            setA.addAll(b);
            setB.addAll(a);
        }
        hc = computeHashCode();
    }


    /** Constructor based on csv strings */
    public Split(String sa, String sb) throws AlgorithmError {
        HashSet<String> a = convertCSVtoSet(sa);
        HashSet<String> b = convertCSVtoSet(sb);
        HashSet x = intersection(a, b);
        if (x.size()>0) throw new AlgorithmError("Attempt to create a Split with intersecting taxa.");
        setA = new HashSet<String>();
        setB = new HashSet<String>();

        String s = Collections.min(a);
        String t = Collections.min(b);
        if (s.compareTo(t)<0) {
            eve = s;
            setA.addAll(a);
            setB.addAll(b);
        }
        else {
            eve = t;
            setA.addAll(b);
            setB.addAll(a);
        }
        hc = computeHashCode();
    }


    /** Get first taxon */
    public String getEve() {
        return eve;
    }


    /** Get all the taxa */
    public HashSet<String> allTaxa() {
        HashSet<String> x = new HashSet<String>();
        x.addAll(setA);
        x.addAll(setB);
        return x;
    }

    public int numTaxa() {
        return setA.size()+setB.size();
    }

    /** Do the sets of all taxa match up? */
    public Boolean checkTaxaMatch(Split anotherSplit) {
        HashSet<String> x = this.allTaxa();
        HashSet<String> y = anotherSplit.allTaxa();
        return ((x.containsAll(y))&&(y.containsAll(x)));
    }

    public int hashCode() {
        return hc;
    }
    private int computeHashCode() {
        ArrayList<String> oa = new ArrayList();
        oa.addAll(setA);
        Collections.sort(oa);
        ListIterator<String> itA = oa.listIterator();
        int r = 17;
        while (itA.hasNext()) {
            r = 31*r+itA.next().hashCode();
        }
        return r;
    }

    /* Test whether two bipartitions are equal
     NB: you might have (A1 = A2 and B1 = B2) or (A1 = B2 and B1 = A2) when the bipartitions are equal */
    public boolean equals(Object p) {
        Split s = (Split) p;
        return equals(s);
    }
    public boolean equals(Split p) {
        return p.containsAll(setA,setB)&&this.containedIn(p);
    }
    private boolean containsAll(HashSet a, HashSet b) {
        return (setA.containsAll(a)&&setB.containsAll(b))||(setA.containsAll(b)&&setB.containsAll(a));
    }
    private boolean containedIn(Split p) {
        return (p.containsAll(setA, setB));
    }

    /** Implement Comparable */
    public int compareTo(Object o) {
        Split p = (Split) o;
        if (hashCode()<p.hashCode()) return -1;
        if (hashCode()>p.hashCode()) return 1;
        return -p.compareSet(setA);
     }
    private int compareSet(HashSet<String> setX) {

        if (setA.size()<setX.size()) return -1;
        if (setA.size()>setX.size()) return 1;

        // Sets same size. Could either be equal, or one contains a different element
        ArrayList<String> oA = new ArrayList();
        oA.addAll(setA);
        ArrayList<String> oX = new ArrayList();
        oX.addAll(setX);
        Collections.sort(oA);
        Collections.sort(oX);
        for (int i=0; i<oA.size(); i++) {
            int d = oA.get(i).compareTo(oX.get(i));
            if (d!=0) return d;
        }
        return 0; 
    }

    /* Test whether the bipartition corresponds to a terminal branch (ie. one ending in a leaf)
     Return the single leaf if this is the case, otherwise return null */
    public String isTerminal() {
        if (setA.size()==1) {
            // Get the contents of set A
            Iterator<String> i = setA.iterator();
            return i.next();
        }
        if (setB.size()==1) {
            // Get the contents of set B
            Iterator<String> i = setB.iterator();
            return i.next();
        }
        return null;
    }


    public String toString() {
        String s = new String();
        s += setA.toString();
        s += "|";
        s += setB.toString();
        return s;
    }

    public String toShortString() {
        String s = new String();
        if (setA.size()<=setB.size()) {
            s += setA.toString();
        }
        else {
            s += setB.toString();
        }
        return s;

    }

    public String toOrderedString() {
        String s = new String();
        ArrayList<String> namesA = new ArrayList<String>();
        ArrayList<String> namesB = new ArrayList<String>();
        namesA.addAll(setA);
        namesB.addAll(setB);
        Collections.sort(namesA);
        Collections.sort(namesB);
        s += namesA.toString();
        s += "|";
        s += namesB.toString();
        return s;
    }

    public String toOrderedShortString() {
        String s = new String();
        ArrayList<String> names = new ArrayList<String>();
        if (setA.size()<=setB.size()) names.addAll(setA);
        else names.addAll(setB);
        Collections.sort(names);
        s += names.toString();
        return s;
    }

    /** Test whether bipartition separates two species */
    public boolean separatesTaxa(String x, String y) throws AlgorithmException {
        if ((!setA.contains(x))&&(!setB.contains(x))) throw new AlgorithmException("Taxon not found in split.");
        if ((!setA.contains(y))&&(!setB.contains(y))) throw new AlgorithmException("Taxon not found in split.");
        if ((setA.contains(x))&&(setB.contains(y))) return true;
        if ((setA.contains(y))&&(setB.contains(x))) return true;
        return false;
    }

    /** Test whether a split is compatible with another in the sense of tree-like topology. */
    public boolean isCompatible(Split p) {
        return p.isCompatible(setA, setB);
    }
    private boolean isCompatible(HashSet a, HashSet b) {
        if ((setA.containsAll(a))&&(b.containsAll(setB))) return true;
        if ((setA.containsAll(b))&&(a.containsAll(setB))) return true;
        if ((setB.containsAll(a))&&(b.containsAll(setA))) return true;
        if ((setB.containsAll(b))&&(a.containsAll(setA))) return true;
        return false;
    }

    /** Return taxa corresponding to one half of the split */
    public HashSet<String> getTaxonSubset(HashSet<String> envelope) {

        HashSet<String> h = new HashSet<String>();
        if (envelope==null) {
            h.addAll(setA);
            return h;
        }

        // If envelope non-null, find a set that is contained in envelope
        if (envelope.containsAll(setA)) {
            h.addAll(setA);
            return h;
        }
        else if (envelope.containsAll(setB)) {
            h.addAll(setB);
            return h;
        }
        else return null;
    }

    /** Partial order induced by a root.
     Leaves should be assigned either setA or setB from the top-level bi-partition
     being used to root a Dendrogram. */
    public boolean partialOrder(Split p, HashSet<String> rootSet) throws AlgorithmException {
        // q = this split
        /* We return true if p<q.
           This happens only when hanging taxa p < hanging taxa q < rootSet

           Return false if q<p
           i.e. hanging taxa q < hanging taxa p < rootSet
         */
        HashSet<String> hp = p.getTaxonSubset(rootSet);
        HashSet<String> hq = this.getTaxonSubset(rootSet);
        if ((hp==null)||(hq==null)) throw new AlgorithmException("Incomparable splits");

        if (hq.containsAll(hp)) return true;
        if (hp.containsAll(hq)) return false;

        throw new AlgorithmException("Incomparable splits");
    }

    /** Test whether a set of splits are <pairwise> compatible */
    public static boolean testPairwiseCompatibility(ArrayList<Split> s) {
        Split p,q;
        for (int i=0; i<s.size(); i++) {
            p = s.get(i);
            for (int j=(i+1); j<s.size(); j++) {
                q = s.get(j);
                //good = (good&&(p.isCompatible(q)));
                if (!p.isCompatible(q)) return false;
            }
        }
        return true;
    }

    /** Test compatibilty with a set of splits */
    public boolean testCompatibilityWithSetofSplits(HashSet<Split> s) {
        Iterator<Split> it = s.iterator();
        Split p;
        for (int i=0; i<s.size(); i++) {
            p = it.next();
            if (!(this.isCompatible(p))) return false;
        }
        return true;
    }

    /** translate the taxa according to some HashMap */
    public Split translate(HashMap<String,String> translationInfo) throws AlgorithmError {
        HashSet<String> a = new HashSet<String>();
        String o;
        Iterator<String> itA = setA.iterator();
        for (int i=0; i<setA.size(); i++) {
            o = itA.next();
            if (!(translationInfo.containsKey(o))) throw new AlgorithmError("Translation failed for a split. Missing taxon was "+o.toString());
            a.add(translationInfo.get(o));
        }
        // Repeat for b
        HashSet<String> b = new HashSet<String>();
        Iterator<String> itB = setB.iterator();
        for (int i=0; i<setB.size(); i++) {
            o = itB.next();
            if (!(translationInfo.containsKey(o))) throw new AlgorithmError("Translation failed for a split. Missing taxon was "+o.toString());
            b.add(translationInfo.get(o));
        }
        return new Split(a,b);
    }

    /** Score two Splits -- see notes from 16-18 may 2005 and Bioinformatics paper */
    public double scoreAgainst(Split p) throws AlgorithmError {
        // Check compatibility of the bipartitions
        HashSet x = this.allTaxa();
        HashSet y = p.allTaxa();
        if(!((x.containsAll(y))&&(y.containsAll(x))))
            throw new AlgorithmError("Error comparing bipartitions.");

        // New score
        double n = x.size();
        HashSet<String> z = p.getTaxonSubset(null);
        HashSet<String> a11 = intersection(setA, z);
        HashSet<String> a21 = intersection(setB, z);
        double n11 = a11.size();
        double n12 = setA.size()-n11;
        double n21 = a21.size();
        double n22 = setB.size()-n21;
        double p11 = 1.0 - n11/(n11+n12+n21);
        double p12 = 1.0 - n12/(n11+n12+n22);
        double p21 = 1.0 - n21/(n11+n21+n22);
        double p22 = 1.0 - n22/(n21+n12+n22);

        double d1 = (p11>p22 ? p11 : p22);
        double d2 = (p12>p21 ? p12 : p21);
        double d = (d1>d2 ? d2 : d1);

        double s = 1.0 - d;

        return s;
    }

    /** Split a string up into taxa */
    public static HashSet<String> convertCSVtoSet(String theString) {
        HashSet<String> h = new HashSet();
        while (!theString.isEmpty()) {
            int ind = theString.indexOf(",");
            if (ind<0) {
                h.add(theString);
                theString = "";
            }
            else {
                String taxon = theString.substring(0, ind);
                theString = theString.substring(ind+1);
                h.add(taxon);
            }
        }
        return h;
    }


    /** Find the intersection of two sets */
    static public HashSet intersection(HashSet a, HashSet b) {
        HashSet x = new HashSet();
        Iterator it = a.iterator();
        for (int i=0; i<a.size(); i++) {
            Object o = it.next();
            if (b.contains(o)) x.add(o);
        }
        return x;
    }

}
