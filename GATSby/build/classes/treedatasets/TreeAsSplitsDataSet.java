/*
 * TreeAsSplitsDataSet.java

    Copyright (C) 2012  Tom M. W. Nye

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

package treedatasets;


import treebase.TreeAsSplits;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.NewickStringUtils;
import treebase.Split;


/**
 * Class representing a set of trees as weighted splits
 */


public class TreeAsSplitsDataSet {

    /* Instance data */
    public ArrayList<TreeAsSplits> theTrees; // Set of trees stored as a map

    public HashSet<Split> allSplits; // The collection of all splits in the data set
    public HashSet<Split> allNonTrivialSplits;
    public HashSet<Split> pendantSplits;

    public ArrayList<String> treeNames;     // Tree names
    public HashSet<String> taxa;            // All taxa

    public HashMap<Split,Double> branchTransFac; // Factor by which branches are transformed

    public int numTrees;
    public int numTaxa;

    public String fileName;


    /* Constructors --------------------------------------------------------- */

    /**
    * Creates a new instance of TreeAsSplitsDataSet from a file
    */
    public TreeAsSplitsDataSet(File theFile) throws IOException {
        // Read until line ends ";"
        BufferedReader br = new BufferedReader(new FileReader(theFile));

        String firstLine = br.readLine();
        firstLine = firstLine.trim();
        if (firstLine.equalsIgnoreCase("#nexus")) {
            buildFromNexusFile(theFile, true);
        }
        else {
            buildFromTreeListFile(theFile, 0, true);
        }

//        if (numTrees<2) {
//            throw new IOException("Not enough trees in input file. ");
//        }

        fileName = theFile.getName();
        branchTransFac = null;
    }
    public TreeAsSplitsDataSet(File theFile, boolean verb) throws IOException {
        // Read until line ends ";"
        BufferedReader br = new BufferedReader(new FileReader(theFile));

        String firstLine = br.readLine();
        firstLine = firstLine.trim();
       
        if (firstLine.equalsIgnoreCase("#nexus")) {
            buildFromNexusFile(theFile, verb);
        }
        else {
            buildFromTreeListFile(theFile, 0, verb);
        }

//        if (numTrees<2) {
//            throw new IOException("Not enough trees in input file. ");
//        }

        fileName = theFile.getName();
        branchTransFac = null;
    }
    public TreeAsSplitsDataSet(File theFile, boolean verb, boolean nexus) throws IOException {
        
        if (nexus) {
            buildFromNexusFile(theFile, verb);
        }
        else {
            buildFromTreeListFile(theFile, 0, verb);
        }

        fileName = theFile.getName();
        branchTransFac = null;
    }


    /* Normalization -------------------------------------------------------- */

    /** Replace branch lengths with normalized version: unit mean length */
    public void branchNormalize() {
        branchTransFac = new HashMap(); // Storage for the transformation factor
        double f;
        // Loop thru each split and calculate mean length
        Iterator<Split> itS = allSplits.iterator();
        for (int i=0; i<allSplits.size(); i++) {
            Split s = itS.next();
            int numTreesWithSplit = 0;
            double totalLength = 0.0;
            // Loop thru trees
            for (int j=0; j<numTrees; j++) {
                TreeAsSplits t = theTrees.get(j);
                if (t.contains(s)) {
                    numTreesWithSplit++;
                    totalLength += t.getSplitLength(s);
                }
            }
            if (totalLength>1.0E-10) {
                f = (1.0*numTreesWithSplit)/totalLength;
                branchTransFac.put(s, new Double(f));
            }
            else {
                f = -1.0;
                branchTransFac.put(s, new Double(-1.0));
            }
            

            // Loop thru trees a second time and transform
            for (int j=0; j<numTrees; j++) {
                TreeAsSplits t = theTrees.get(j);
                if (t.contains(s)) {
                    try {
                        if (f>0.0)
                            t.scaleSplitLength(s,f);
                        else
                            t.setSplitLength(s, 1.0);
                    }
                    catch (AlgorithmError anErr) {
                        System.out.println("Error scaling a split. This should not be possible. ");
                    }
                }
            }
        }
    }

    /** Invert the normalisation */
    public void inverseBranchNormalize() {
        if (branchTransFac!=null) {
            Iterator<Split> itS = allSplits.iterator();
            for (int i=0; i<allSplits.size(); i++) {
                Split s = itS.next();
                double f = branchTransFac.get(s).doubleValue();
                // Loop thru trees and transform
                for (int j=0; j<numTrees; j++) {
                    TreeAsSplits t = theTrees.get(j);
                    if (t.contains(s)) {
                        try {
                            t.scaleSplitLength(s,1.0/f);
                        }
                        catch (AlgorithmError anErr) {
                            System.out.println("Error scaling a split. This should not be possible. ");
                        }
                    }
                }
            }

        }
    }


    /** Scale each tree to unit Euclidean norm */
    public void treeNormalize() {
        // Storage
        double x,u,w;
        // Loop thru' trees
        for (int i=0; i<numTrees; i++) {
            // Work out norm of this tree
            u = 0.0;
            TreeAsSplits t = theTrees.get(i);
            HashSet<Split> splits = t.getSplits();
            Iterator<Split> it = splits.iterator();
            for (int j=0; j<splits.size(); j++) {
                Split p = it.next();
                x = t.getSplitLength(p);
                u += x*x;
            }
            w = Math.sqrt(u);
            // Then scale
            t.scaleAllLengths(1/w);
        }

    }
    
    /* Utils ---------------------------------------------------------------- */

    /* General */

    public TreeAsSplits getTree(int k) {
        return theTrees.get(k);
    }

    /* Counting methods */
    public int countTreeTopologies() {

        ArrayList<HashSet<Split>> theTopologies = new ArrayList();

        for (int i=0; i<numTrees; i++) {
            HashSet<Split> newTop = theTrees.get(i).getSplits();
            boolean found = false;
            for (int j=0; j<theTopologies.size(); j++) {
                HashSet<Split> oldTop = theTopologies.get(j);
                if ((oldTop.containsAll(newTop))&&(newTop.containsAll(oldTop))) {
                    found = true;
                }
            }
            if (!found) theTopologies.add(newTop);
        }
        return theTopologies.size();
    }
    
    
    /* Read in from files --------------------------------------------------- */
    
       /** Build from a file containing a list of trees.
        This is also used by the nexus input function: hence this function includes a
        starting line (so nexus header can be ignored), nexus comment removal
        (blocks in []), and check for lines "end;" */
    private void buildFromTreeListFile(File theFile, int startingLine, boolean verb) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(theFile));
        // Skip lines as required
        for (int i=0; i<startingLine; i++) {
            String dummy = br.readLine();
        }

        // Read until line ends ";"
        treeNames = new ArrayList();
        String firstTreeStr = readSingleTree(br);
        String tempTreeStr;
        try {
            tempTreeStr = getTreeNameAndTidyString(firstTreeStr, 0);
            taxa = extractTaxonNamesFromString(tempTreeStr);
        }
        catch (NewickStringUtils.ParsingException anError) {
            throw new IOException(anError.getMessage());
        }

    /* ---------------------------------------------------------------------- */

    /*  PASS THRU' FILE:
        Extract:
        1. Tree names
        2. Complete set of splits
        3. ArrayList of TreeAsMaps
    */

        // Initialize
        treeNames.clear(); // There might already be a name if the first line contained one
        allSplits = new HashSet();
        theTrees = new ArrayList();
        allNonTrivialSplits = new HashSet();
        pendantSplits = new HashSet();
        numTrees = 0;
        numTaxa = taxa.size();

        // Now parse each string
        String treeString = firstTreeStr;
        while (treeString != null) {
            // Get tree name
            try {
                 treeString = getTreeNameAndTidyString(treeString, treeNames.size());
// System.out.println("Extracting splits from "+treeNames.get(numTrees));
            }
            catch (NewickStringUtils.ParsingException anError) {
                throw new IOException(anError.getMessage());
            }

            TreeAsSplits splitsAndLengths = new TreeAsSplits(numTaxa);
            try {
                extractSplitsAndLengthsFromString(splitsAndLengths, treeString, taxa);
            }
            catch (NewickStringUtils.ParsingException anError) {
                throw new IOException(anError.getMessage());
            }
            
            theTrees.add(splitsAndLengths);

            // Add all splits to the collection of splits.
            // NB: no check on positivity
            allSplits.addAll(splitsAndLengths.getSplits());

            numTrees++;

            treeString = readSingleTree(br);
            // Quit if treeString is "end;"
            if (treeString!=null) {
                if (treeString.startsWith("end;")||treeString.startsWith("END;")||treeString.startsWith("End;")) treeString=null;
            }
        } // End while loop
        
    /* END OF PASS THRU THE FILE */

    /* ---------------------------------------------------------------------- */

        /* Get non trivial splits */
        
        allNonTrivialSplits = new HashSet();
        pendantSplits = new HashSet();
        Iterator<Split> it = allSplits.iterator();
        for (int i=0; i<allSplits.size(); i++) {
            Split p = it.next();
            if (p.isTerminal()==null) allNonTrivialSplits.add(p);
            else pendantSplits.add(p);
        }

        // Done!
        if (verb) {
            System.out.println("There are "+allSplits.size()+" different splits from "+numTrees+" trees on "+numTaxa+" taxa.");
            System.out.println("Data set contains "+this.countTreeTopologies()+" different topologies.");
        }
        
        /* Check resolution */
        int numUnresolved = 0;
        int s = 2*numTaxa-3;
        for (int i=0; i<numTrees; i++) {
            TreeAsSplits t = theTrees.get(i);
            if (t.getNumSplits()<s) numUnresolved++;
        }
        if ((verb)&&(numUnresolved>0)) System.out.println(numUnresolved+" trees in the data set were unresolved.");
    }
    
    
    /** Build from a file containing nexus list */
    private void buildFromNexusFile(File theFile, boolean verb) throws IOException {
        if (verb) System.out.println("Reading nexus file. \nWarning: nexus file must contain a list of trees. Fully general nexus format is not handled.\n");

        // Read until line is "begin trees;"
        int lineCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String line;
        boolean readingTrees = false;
        while (!readingTrees) {
            line = br.readLine();
            if (line==null) throw new IOException("File ended before there was a list of trees.");
            line = line.trim();
            if (line.startsWith("begin trees;")||line.startsWith("BEGIN TREES;")||line.startsWith("Begin trees;")||line.startsWith("Begin Trees;")) readingTrees = true;
            lineCount++;
        }

        // Read until either "translate" or trees start
        boolean translation = false;
        HashMap<String,String> translationInfo = new HashMap();
        boolean validLine = false;
        while (!validLine) {
            line = br.readLine();
            if (line==null) throw new IOException("File ended before there was a list of trees.");
            line = removeCommentFromNexusString(line);
            line = line.trim();
            if (!line.equals("")) validLine = true;
            if (line.startsWith("translate")||line.startsWith("TRANSLATE")||line.startsWith("Translate")) {
                line = line.substring(9);
                translation = true;
            }
            lineCount++;
        }

        if (translation) {
            // Read in translation information
            // Get all info in a single string
            String s;
            String infoString = "";
            while (!(infoString.endsWith(";"))) {
                s=br.readLine();
                infoString += removeCommentFromNexusString(s);
                infoString = infoString.trim();
                lineCount++;
            }

            // Parse the translation string
            // Remove trailing ";"
            infoString = infoString.substring(0, infoString.length()-1);
            // Loop thru' comma separated entries
            while (!infoString.equals("")) {
                int index = infoString.indexOf(',');
                String chunk;
                if (index<0) {
                    chunk = new String(infoString);
                    infoString = "";
                }
                else {
                    chunk = new String(infoString.substring(0,index));
                    infoString = infoString.substring(index+1);
                }
                chunk = chunk.trim();

                index = chunk.indexOf(' ');
                if (index<0) index = chunk.indexOf("\t");
                String key, taxon;
                if (index<0) {
                    throw new IOException("Error in nexus format. Translation block invalid.");
                }
                else {
                    key = new String(chunk.substring(0, index));
                    taxon = new String(chunk.substring(index));
                    key = key.trim();
                    taxon = taxon.trim();
                }
                translationInfo.put(key, taxon);
            } // End loop thru translation info comma separated blocks

        } // End translation info
        else {
            lineCount--;
            translationInfo = null;
        }

        // Read in trees from current point
        this.buildFromTreeListFile(theFile, lineCount, verb);

        if (translation) {
            // Back-translate the taxa / splits
            HashSet<String> newTaxonSet = translateTaxa(taxa, translationInfo);
            taxa.clear();
            taxa.addAll(newTaxonSet);

            // Now do store of splits
            HashSet<Split> newSplits = new HashSet();
            Iterator<Split> it = allSplits.iterator();
            for (int i=0; i<allSplits.size(); i++) {
                Split p = it.next();
                Split q = null;
                try {
                    q = p.translate(translationInfo);
                }
                catch (AlgorithmError anErr) {
                    System.out.println("Algorithm error: couldn't back-translate from a nexus file.");
                }
                newSplits.add(q);
            }
            allSplits.clear();
            allSplits.addAll(newSplits);

            // Now do non-triv splits again
            allNonTrivialSplits.clear();
            pendantSplits.clear();
            it = allSplits.iterator();
            for (int i=0; i<allSplits.size(); i++) {
                Split p = it.next();
                if (p.isTerminal()==null) allNonTrivialSplits.add(p);
                else pendantSplits.add(p);
            }

            /* Finally loop thru all trees and convert splits */
            for (int i=0; i<theTrees.size(); i++) {

                TreeAsSplits splitsAndLengths = theTrees.get(i);
                try {
                    splitsAndLengths.translateTaxa(translationInfo);
                }
                catch (AlgorithmException anotherError) {
                    System.out.println("Error: Failed to find a split when reading nexus tree file a second time. This should not be possible."+anotherError.getMessage());
                }

            }
        }

    }
    
    
    /* File utils ----------------------------------------------------------- */
    
    /** Remove comment from a nexus string */
    private static String removeCommentFromNexusString(String theString) {
        int a,b;
        String s = new String(theString);
        boolean containsComment = true;
        while (containsComment) {
            a = s.indexOf('[');
            if (a<0) containsComment = false;
            else {
                b = s.indexOf(']', a);
                if (b<0) return s;
                else {
                    String t =new String( s.substring(0, a)+s.substring(b+1) );
                    s = t;
                }
            }
        }
        s = s.trim();
        return s;
    }


    /** Read a tree string from a buffered reader */
    private static String readSingleTree(BufferedReader br) throws IOException {
        String s;
        String theString = "";
        while (!(theString.endsWith(";"))) {
            if ( (s=br.readLine()) == null) return null;
            else {
                theString += removeCommentFromNexusString(s);
            }
            theString = theString.trim();
        }
        return theString;
    }

    /** Get a name from a string and tidy it up (remove ";" and any outer brackets and trailing number) */
    protected String getTreeNameAndTidyString(String treeString, int n) throws NewickStringUtils.ParsingException {

        // Get name
        String theName = null;
        int posEq = treeString.indexOf("=");
        if (posEq<0) {
            theName = new String("Tree ");
            theName = theName+Integer.toString(n+1);
        }
        else {
            theName = new String(treeString.substring(0,posEq));
            theName = theName.trim();
            int start = treeString.indexOf("(", posEq);
            treeString = treeString.substring(start);
        }
        treeNames.add(theName);
 
        // Tidy string
        
        // Delete final ";"
        int index = treeString.indexOf(";");
        if (index>-1) {
            treeString = treeString.substring(0, index);
        }

        // Remove outer brackets and any trailing number
        int endPos = 1;
        endPos = NewickStringUtils.findMatchingBracket(0, treeString);
        treeString = treeString.substring(1, endPos);

        return treeString;
    }


    /** Get all the taxon names from a Newick string */
    public static HashSet<String> extractTaxonNamesFromString(String theString) throws NewickStringUtils.ParsingException {
        HashSet<String> u = new HashSet();
        String s = new String(theString);
        // Delete final ";"
        int index = s.indexOf(";");
        if (index>-1) {
            s = s.substring(0, index);
        }

        while (s.length()>0) {
            if ((s.startsWith(" "))||(s.startsWith("("))||(s.startsWith(","))||(s.startsWith(")"))) {
                // Remove first character
                s = s.substring(1);
            }
            else {
                // Find index of next colon
                index = s.indexOf(":");
                if (index == -1) {
                    throw new NewickStringUtils.ParsingException("Error parsing dendrogram string to extract taxa. " +
                                                               "Unable to find branch length.","Missing colon");
                }
                String t = new String(s.substring(0, index));
                u.add(t);
                // Delete up to next comma or close bracket
                s = s.substring(index+1); // up to colon first
                index = NewickStringUtils.getOpenBracketOrCommaPos(s);
                if (index<0) {
                    // Reached end of string
                    s = "";
                }
                else {
                    s = s.substring(index);
                }
            }
        }

        return u;
    }


    /** Turn a Newick string into a hashmap from split to branch length.
     Function expects a comma separated list of Newick strings. 
     Calls itself recursively, hence void. */
    protected static void extractSplitsAndLengthsFromString(TreeAsSplits splitsAndLengths, String theString, HashSet<String> allTaxa) throws NewickStringUtils.ParsingException {

        int startPos, endPos, colonPos;

        while (theString.length()>0) {
            startPos = 0;

            // If leading bracket, find matching close
            if (theString.charAt(startPos)=='(') {
                endPos = NewickStringUtils.findMatchingBracket(startPos,theString);

                if (theString.charAt(endPos+1)==':') {
                    endPos = endPos+1;
                    startPos = endPos+1;
                }
                else {
                    // There must be a bootstrap proportion here
                    colonPos = theString.indexOf(':', endPos);
                    if (colonPos<0) throw new NewickStringUtils.ParsingException("Missing colon after bracket and bootstrap proportion in Newick string","Missing colon");
                    endPos = endPos+1;
                    startPos = colonPos+1;
                }

            }
            else {
                // The current element in the list is a single taxon
                endPos = theString.indexOf(":");
                startPos = endPos+1;
                if (endPos<0) {
                    throw new NewickStringUtils.ParsingException("Missing colon after taxon name in Newick string","Missing colon");
                }
            }

            String subTreeString = theString.substring(0,endPos);
            theString = theString.substring(startPos);

            // Get branch length
            int endOfNumberPos = NewickStringUtils.getCloseBracketOrCommaPos(theString);
            if (endOfNumberPos<0) endOfNumberPos = theString.length();
            // Extract the string containing the number
            String numString = theString.substring(0,endOfNumberPos);
            // Remove the number from the start of the string
            theString = theString.substring(endOfNumberPos);
            if (theString.startsWith(",")) theString = theString.substring(1);
            theString = theString.trim();

            // Convert the number string into a branch length
            double t = 0.0;
            try {
                t = Double.parseDouble(numString);
            }
            catch (NumberFormatException anErr) {
                throw new NewickStringUtils.ParsingException("Error parsing newick string. A branch length for a leaf node was formatted incorrectly.","Bad branch length");
            }

            /* OK. Now got a sub-tree string, a branch length and remainder of string
             1. make a split
             2. put into hash map
             3. parse next bit of string */

            // Make the corresponding split
            String treeWithLengthString = subTreeString+":"+numString;
            HashSet<String> setA = extractTaxonNamesFromString(treeWithLengthString);
            HashSet<String> setB = new HashSet();
            setB.addAll(allTaxa);
            setB.removeAll(setA);
            Split p = null;
            try {
                p = new Split(setA, setB);
            }
            catch (AlgorithmError anErr) {
                System.out.println("Algorithm error: couldn't make a split while parsing a tree.");
                throw new NewickStringUtils.ParsingException("Error parsing dendrogram string. Couldn't make split.","Error extracting taxon names");
            }

            if (splitsAndLengths.contains(p)) {
                double x = splitsAndLengths.getSplitLength(p);
                try {
                    splitsAndLengths.setSplitLength(p,x+t);
                }
                catch (AlgorithmError anErr) {
                    System.out.println("AlgorithmError setting split length. This should not be possible.");
                }
            }
            else {
                splitsAndLengths.add(p, t); // No error check! p should already be compatible
            }
            
            if (subTreeString.length()>0) {
                // if subTreeString begins and ends with brackets then keep going
                if ((subTreeString.startsWith("("))&&(subTreeString.endsWith(")"))) {
                    subTreeString = subTreeString.substring(1,subTreeString.length()-1);
                    extractSplitsAndLengthsFromString(splitsAndLengths, subTreeString, allTaxa);
                }
            }
        } // End while loop, parsing sections of the string

        // OK done -- return the HashMap
    }

    /* Conversion utils ----------------------------------------------------- */

    /*
    public static Split findEquivalentSplitInCollection(Split p, Collection<Split> h) {
        Iterator<Split> it = h.iterator();
        for (int i=0; i<h.size(); i++) {
            Split q = it.next();
            if (p.equals(q)) return q;
        }
        return null;
    }
    
    public Split findEquivalentSplitInData(Split p) {
        return findEquivalentSplitInCollection(p, allSplits);
    }
    */
    
    /** Translate a hashset of taxon labels */
    private static HashSet<String> translateTaxa(HashSet<String> x, HashMap<String,String> translationInfo) {
        HashSet<String> y = new HashSet();
        Iterator<String> it = x.iterator();
        for (int i=0; i<x.size(); i++) {
            String s =  it.next();
            String t = translationInfo.get(s);
            y.add(t);
        }
        return y;
    }


    /** Output */
    private void outputToNexusFile(File theFile) throws java.io.IOException {
        ArrayList<String> orderedTaxa = new ArrayList();
        orderedTaxa.addAll(taxa);
        Collections.sort(orderedTaxa);

        BufferedWriter writer = new BufferedWriter(new FileWriter(theFile));
        writer.write("#NEXUS\n");
        writer.write("\n");

        writer.write("begin taxa;\n");
        writer.write("    dimensions ntax="+numTaxa+"\n");
        writer.write("    taxlabels\n");
        for (int i=0; i<numTaxa; i++) {
            writer.write("        "+orderedTaxa.get(i)+"\n");
        }
        writer.write("        ;\n");
        writer.write("end;\n");

        HashMap<String,String> trans = new HashMap();
        writer.write("begin trees;\n");
        writer.write("    translate\n");
        for (int i=0; i<numTaxa; i++) {
            writer.write("        "+(i+1)+"\t"+orderedTaxa.get(i));
            if (i<(numTaxa-1)) writer.write(",\n");
            else writer.write("\n");
            trans.put(orderedTaxa.get(i), new String(Integer.toString(i+1)));
        }
        writer.write("        ;\n");
        try {
            for (int i=0; i<numTrees; i++) {
                TreeAsSplits t = theTrees.get(i).efficientClone();
                t.translateTaxa(trans);
                writer.write("    "+treeNames.get(i)+" = "+t.toString()+"\n");
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error translating taxa when outputting to nexus file.");
        }
        writer.write("end;\n");
        writer.close();
    }


}
