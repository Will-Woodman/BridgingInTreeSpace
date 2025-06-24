/**
    Alignment
    Class for representing sequence alignments

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

package sequencemodels;

/*

 */

import java.util.*;
import treebase.*;
import java.io.*;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.DoubleUniform;
import simulation.CategoricalDistribution;

public class Alignment implements java.io.Serializable {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    protected int[][] data; // [taxon][site] set of alignment data
    protected int numSites, numTaxa;
    protected ArrayList taxa;
    protected Alphabet theAlphabet;
    
    protected int numSitePatterns;
    protected int[] numSitesPerPattern;
    protected ArrayList<ArrayList<Integer>> sitesPerPattern;
    protected boolean isCondensed;

    /** Blank constructor */
    public Alignment(int nt, int nc, Alphabet a) {
        data = new int[nt][nc];
        numSites = nc;
        numTaxa = nt;
        theAlphabet = a;
        taxa = new ArrayList(numTaxa);
        numSitePatterns = numSites;
        numSitesPerPattern = new int[numSites];
        Arrays.fill(numSitesPerPattern, 1);
        sitesPerPattern = new ArrayList<ArrayList<Integer>>();
        for(int i=0; i<numSites; i++) {
            ArrayList<Integer> tmp = new ArrayList<Integer>();
            tmp.add(i);
            sitesPerPattern.add(tmp);
        }
        isCondensed = false;
    }
    
    public void condenseAlignment() {
        // Identify distinct patterns
        ArrayList<int[]> condensedDataTmp = new ArrayList<int[]>();
        condensedDataTmp.add(getColumn(0));
        sitesPerPattern = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> theSites = new ArrayList<Integer>();
        theSites.add(0);
        sitesPerPattern.add(theSites);
        for(int i=1; i<getNumSites(); i++) {
            int[] col = getColumn(i);
            for(int j=0; j<condensedDataTmp.size(); j++) {
                if(Arrays.equals(col, condensedDataTmp.get(j))) {
                    sitesPerPattern.get(j).add(i);
                    break;
                }
                if(j==condensedDataTmp.size()-1) {
                    condensedDataTmp.add(col); 
                    theSites = new ArrayList<Integer>();
                    theSites.add(i);
                    sitesPerPattern.add(theSites);
                    break;
                }
            }
        }
        numSitePatterns = condensedDataTmp.size();
        numSitesPerPattern = new int[numSitePatterns];
        for(int i=0; i<numSitePatterns; i++) numSitesPerPattern[i] = sitesPerPattern.get(i).size();
        data = new int[numTaxa][numSitePatterns];
        for(int i=0; i<numSitePatterns; i++) {
            for(int j=0; j<numTaxa; j++) {
                data[j][i] = condensedDataTmp.get(i)[j];
            }
        }
        isCondensed = true;
    }

    /* UTILS */

    public int getNumTaxa() {
        return numTaxa;
    }
    public int getNumSites() {
        return numSites;
    }
    public Alphabet getAlphabet() {
        return theAlphabet;
    }
    public boolean matchesAlphabet(Alignment a) {
        return (theAlphabet==a.getAlphabet());
    }
    public ArrayList getTaxa() {
        ArrayList a = new ArrayList();
        a.addAll(taxa);
        return a;
    }
    public boolean matchesTaxonSet(Alignment a) {
        ArrayList t = a.getTaxa();
        return ((t.containsAll(taxa))&&(taxa.containsAll(t)));
    }

    public int[] getColumn(int k) {
        int[] res = new int[numTaxa];
        for (int i=0; i<numTaxa; i++) res[i] = data[i][k];
        return res;
    }

    /** Set a single letter */
    public void setLetter(String name, int siteNum, int letter) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        int rowNum = taxa.indexOf(name);
        if (rowNum<0) throw new AlgorithmException("Taxon name not found.");
        data[rowNum][siteNum] = letter;
    }

    /** Set a specific column */
    public void setColumn(int colNum, int[] d) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        for (int i=0; i<numTaxa; i++) {
            data[i][colNum] = d[i];
        }
    }

    /** Set a specific row */
    public void setRow(int rowNum, int[] d) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        for (int i=0; i<numSites; i++) {
            data[rowNum][i] = d[i];
        }
    }
    public void setRow(String name, int[] d) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        int rowNum = taxa.indexOf(name);
        if (rowNum<0) throw new AlgorithmException("Taxon name not found.");
        setRow(rowNum,d);
    }
    public void setRow(int rowNum, String s) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        int[] d = theAlphabet.translateLettersToIndex(s);
        if (d.length!=numSites) throw new AlgorithmException("Wrong number of sites setting alignment row.");
        setRow(rowNum, d);
    }
    public void setRow(int rowNum, String name, String s) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        if (rowNum>=taxa.size()) {
            taxa.add(rowNum, name);
        }
        else {
            taxa.set(rowNum, name);
        }

        setRow(rowNum, s);
    }
    public void setRow(String name, String s) throws AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot manipulate an alignment once it has been "
                + "condensed.");
        int rowNum = taxa.indexOf(name);
        if (rowNum<0) throw new AlgorithmException("Taxon name not found.");
        setRow(rowNum, s);
    }

    public void setTaxa(ArrayList t) throws AlgorithmError {
        if (t.size()!=numTaxa) throw new AlgorithmError("Wrong number of taxa.");
        taxa = new ArrayList();
        taxa.addAll(t);
    }

    /** Get taxon name */
    public String getTaxon(int i) {
        return (String) taxa.get(i);
    }

    /** Set taxa names */
    public void setTaxonNames(ArrayList names) throws AlgorithmError {
        if (names.size()!=numTaxa) throw new AlgorithmError("Taxon number mismatch when altering an alignment.");
        taxa.clear();
        taxa.addAll(names);
    }

    /** Get a row as a string */
    public String getRowString(int k) {
        try {
            String s = theAlphabet.translateIndicesToLetters(data[k]);
            return s;
        }
        catch (AlgorithmException e) {
            javax.swing.JOptionPane.showMessageDialog(null,"Error translating data to file "+e.toString(),"Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }
        return null;
    }
    public String getRowStringWithName(int k) {
        String s = getRowString(k);
        String t = new String((String)taxa.get(k));
        t += " "+s;
        return t;
    }

    /** Get row as a sequence */
    public int[] getRow(int k) {
        return data[k];
    }
    public int[] getRow(String name) throws AlgorithmError {
        int i = taxa.indexOf(name);
        if (i<0) throw new AlgorithmError("Taxon on tree was missing from alignment. ");
        return getRow(i);
    }
    /** Get a particular letter */
    public int getLetter(String name, int k) throws AlgorithmError {
        int i = taxa.indexOf(name);
        if (i<0) throw new AlgorithmError("Taxon on tree was missing from alignment. ");
        return data[i][k];
    }
    
    public int getNumSitePatterns() {
        return numSitePatterns;
    }
    public int getNumSitesPerPattern(int k) {
        return numSitesPerPattern[k];
    }
    public int[] getSitesPerPattern(int k) {
        int[] sites = new int[sitesPerPattern.get(k).size()];
        for(int i=0; i<sitesPerPattern.get(k).size(); i++) sites[i] = sitesPerPattern.get(k).get(i);
        return sites;
    }
    public boolean getIsCondensed() {
        return isCondensed;
    }

    /* Concatenation and resampling ----------------------------------------- */

    /** If alignment has been condensed, returns column number in condensed alignment
        corresponding to column k in the original */
    private int getColumnForSite(int k) {
        if(!isCondensed) return k;
        else {
            int l = -1;
            for(int i=0; i<numSitePatterns; i++) {
                for(int j=0; j<sitesPerPattern.get(i).size(); j++) {
                    if(k==sitesPerPattern.get(i).get(j)) {
                        l = i;
                        break;
                    }
                }
            }
            return l;
        }
    }
    
    /** Truncate an alignment: first site and last site must lie between 1 and a.getNumSites() */
    public static Alignment truncate(Alignment a, int firstSite, int lastSite) throws AlgorithmException {
        if (firstSite<1 || firstSite> a.getNumSites() || lastSite<1 || lastSite> a.getNumSites() || lastSite < firstSite) 
            throw new AlgorithmException("Invalid truncation");
        Alignment x = new Alignment(a.getNumTaxa(), lastSite - firstSite + 1, a.getAlphabet());
        x.setTaxa(a.getTaxa());
        
        int j = 0;
        for (int i=firstSite-1; i<lastSite; i++) {
            int col = a.getColumnForSite(i);
            x.setColumn(j, a.getColumn(col));
            j++;
        }

        return x;
    }
    
    /** Concatenate two alignments */
    public static Alignment concat(Alignment a, Alignment b) throws AlgorithmException {
        if (!a.matchesAlphabet(b)) throw new AlgorithmException("Alphabet mismatch between alignments");
        if (!a.matchesTaxonSet(b)) throw new AlgorithmException("Taxon set mismatch between alignments");

        Alignment x = new Alignment(a.getNumTaxa(), a.getNumSites()+b.getNumSites(), a.getAlphabet());

        for (int i=0; i<x.getNumSites(); i++) {
            if (i<a.getNumSites()) {
                int col = a.getColumnForSite(i);
                x.setColumn(i, a.getColumn(col));
            }
            else {
                int col = b.getColumnForSite(i-a.getNumSites());
                x.setColumn(i, b.getColumn(col));
            }
        }

        return x;
    }

    /** Bootstrap sampler */
    public Alignment bootstrap(DoubleMersenneTwister theRandomEngine) throws AlgorithmException {
        Alignment x = new Alignment(numTaxa, numSites, theAlphabet);
        try {
            x.setTaxa(taxa);
        }
        catch (AlgorithmError anErr) {} // Should be impossible

        // Create random number generators etc
        DoubleUniform unif = new DoubleUniform(theRandomEngine);

        // Loop thru sites
        double[] prob = new double[numSitesPerPattern.length];
        for(int i=0; i<prob.length; i++) prob[i] = ((double)numSitesPerPattern[i])/prob.length;
        simulation.CategoricalDistribution catSampler = new simulation.CategoricalDistribution(prob);
        int k;
        for (int i=0; i<numSites; i++) {
            // Set k to be random column number
            k = catSampler.sample();
            x.setColumn(i, getColumn(k));
        }
        return x;
    }

    /* Input from file ------------------------------------------------------- */

    /** Read from  file */
    public static Alignment readFromFile(File theFile, Alphabet a) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String firstLine = br.readLine();
        br.close();
        firstLine = firstLine.trim();
        Alignment x;
        if (firstLine.equalsIgnoreCase("#nexus")) {
            x = readFromNexusFile(theFile);
        }
        else {
            x = readFromPhylipFile(theFile, a);
        }
        return x;
    }

    /** Read from phylip file */
    public static Alignment readFromPhylipFile(File theFile, Alphabet a) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(theFile));

        // Read in first line and extract the number of taxa and number of characters
        String firstLine = br.readLine();
        firstLine = firstLine.trim();
        int spaceInd = firstLine.indexOf(" ");
        if (spaceInd<0) throw new IOException("First line of phylip file cannot be parsed into two integers.");
        int nt = Integer.parseInt(firstLine.substring(0, spaceInd));
        int ns = Integer.parseInt(firstLine.substring(spaceInd+1));

        Alignment x = new Alignment(nt, ns, a);

        // Loop through one line per taxon
        for (int i=0; i<nt; i++) {
            String dataLine = br.readLine();
            String name, d;

            /* Extract taxon name:
             Either text upto to a space or first 10 charcters if there is no space */
            spaceInd = dataLine.indexOf(" ");
            spaceInd = (spaceInd<0) ? dataLine.indexOf("\t") : spaceInd; // Might be a tab instead of a space?
            spaceInd = (spaceInd<0) ? 9 : spaceInd;
            d = dataLine.substring(spaceInd+1);
            name = dataLine.substring(0, spaceInd);
            d = d.trim();

            try {
                x.setRow(i, name, d);
            }
            catch (AlgorithmException e) {
                throw new IOException("Alignment translation error on row "+(i+1)+" of file. Expecting data of type "+a.getClass().toString());
            }
        }

        br.close();
        return x;
    }


    /** Read from Nexus file */
    public static Alignment readFromNexusFile(File theFile) throws IOException {
        Alphabet a = null;

        // Read until line is "begin data;"
        int lineCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String line;
        boolean readingHeader = true;
        while (readingHeader) {
            line = br.readLine();
            if (line==null) throw new IOException("File ended before there was any sequence data.");
            line = line.trim();
            if ((line.startsWith("begin data;"))||(line.startsWith("Begin data;"))||(line.startsWith("Begin Data;"))||(line.startsWith("BEGIN DATA;"))) readingHeader = false;
            lineCount++;
        }
        boolean readingSequences = false;
        int nt=0, nc=0;
        while (!readingSequences) {
            line = br.readLine();
            if (line==null) throw new IOException("File ended before there was any sequence data.");
            line = line.trim();
            if ((line.startsWith("matrix"))||(line.startsWith("Matrix"))||(line.startsWith("MATRIX"))) readingSequences = true;
            lineCount++;
            // Parse number of taxa and characters
            int ntInd = line.indexOf("ntax=");
            if (ntInd<0) ntInd = line.indexOf("NTAX=");
            if (ntInd>0) {
                int ntEnd = line.indexOf(" ", ntInd+6);
                ntEnd = (ntEnd<0) ? line.indexOf(";", ntInd+6) : ntEnd;
                ntEnd = (ntEnd<0) ? line.length() : ntEnd;
                nt = Integer.parseInt(line.substring(ntInd+5,ntEnd));
            }
            int ncInd = line.indexOf("nchar=");
            if (ncInd<0) ncInd = line.indexOf("NCHAR=");
            if (ncInd>0) {
                int ncEnd = line.indexOf(" ", ncInd+7);
                ncEnd = (ncEnd<0) ? line.indexOf(";", ncInd+7) : ncEnd;
                ncEnd = (ncEnd<0) ? line.length() : ncEnd;
                nc = Integer.parseInt(line.substring(ncInd+6,ncEnd));
            }
            int abInd = line.indexOf("datatype=");
            if (abInd<0) abInd = line.indexOf("DATATYPE=");
            if (abInd>0) {
                int abEnd = line.indexOf(" ", abInd+10);
                abEnd = (abEnd<0) ? line.indexOf(";", abInd+10) : abEnd;
                abEnd = (abEnd<0) ? line.length() : abEnd;
                String ab = line.substring(abInd+9,abEnd).trim();
                if (ab.startsWith("dna")||ab.startsWith("DNA")) a = DNAAlphabet.getInstance();
                if (ab.startsWith("protein")||ab.startsWith("PROTEIN")) a = AminoAcidAlphabet.getInstance();
                if (ab.startsWith("dayhoff6")||ab.startsWith("DAYHOFF6")) a = DayhoffSixAlphabet.getInstance();
            }
        }

        if (a==null) throw new IOException("Unable to identify sequence type from nexus file.");
        Alignment x = new Alignment(nt, nc, a);
        
        // Create a hashmap from species names to sequence data
        HashMap sequences = new HashMap();
        ArrayList ot = new ArrayList();
        boolean done = false;
        int ind;
        while (!done) {
            line = br.readLine();
            if (line == null) break;
            line = line.trim();
            done = ((line.startsWith(";"))||(line.startsWith("end;"))||(line.startsWith("END;")));
            if (done) break;
            ind = line.indexOf(" ");
            if ((line.length()>2)&&(ind>0)) { // Check line isn't just blank!
                String name = line.substring(0, ind);
                line = line.substring(ind).trim();
                if (sequences.containsKey(name)) {
                    String seq = (String) sequences.get(name);
                    seq = seq+line;
                    sequences.put(name, seq);
                }
                else {
                    String seq = new String(line);
                    sequences.put(name, seq);
                    ot.add(name);
                }
            }
        }

        if (sequences.keySet().size()!=x.getNumTaxa()) throw new IOException("Number of taxa didn't match dimension info in nexus file.");
        
        try {
            for (int i=0; i<ot.size(); i++) {
                String name = (String) ot.get(i);
                String seq = (String) sequences.get(name);
                seq = seq.replaceAll(" ", ""); // Kill off any spaces
                x.setRow(i, name, seq);
            }
        }
        catch (AlgorithmException e) {
            throw new IOException("Error translating sequences when reading nexus file. "+e.getMessage());
        }

        br.close();
        return x;
    }

    /* Output to file ------------------------------------------------------- */

    public void writeToPhylipFile(File theFile) throws IOException, AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot write condensed alignment to file.");
        BufferedWriter bw = new BufferedWriter(new FileWriter(theFile));
        bw.write(numTaxa+" "+numSites+"\n");
        int k;
        for (int i=0; i<numTaxa; i++) {
            String name = getTaxon(i);
            k = 40-name.length();
            for (int j=0; j<k; j++) {
                // Add on a space
                name = name + " ";
            }
            bw.write(name.substring(0, 40)+getRowString(i)+"\n");
        }
        bw.close();
    }

    public void writeToNexusFile(File theFile) throws IOException, AlgorithmException {
        if(isCondensed) throw new AlgorithmException("Cannot write condensed alignment to file.");
        BufferedWriter bw = new BufferedWriter(new FileWriter(theFile));
        bw.write("#NEXUS\n");
        bw.write("\n");
        bw.write("begin data;\n");
        bw.write("dimensions ntax="+numTaxa+" nchar="+numSites+";\n");
        String typeStr = "unknown";
        if (DNAAlphabet.class.isInstance(theAlphabet)) typeStr = "dna";
        if (AminoAcidAlphabet.class.isInstance(theAlphabet)) typeStr = "protein";
        if (DayhoffSixAlphabet.class.isInstance(theAlphabet)) typeStr = "dayhoff6";
        bw.write("format interleave datatype="+typeStr+" gap="+theAlphabet.gap+" symbols=\""+theAlphabet.getAllLetters()+"\";\n");
        bw.write("\nmatrix\n");

        // Write out data
        int site = 0, blockCount, currentSite;
        try {
            while (site<numSites) {

                // Loop thru taxa
                int lineCount = 0;
                for (int i=0; i<numTaxa; i++) {
                    lineCount = 0;

                    // Write out taxon name
                    String name = getTaxon(i);
                    int numSpaces = 40-name.length();
                    for (int l=0; l<numSpaces; l++) name += " ";
                    bw.write(name.substring(0, 40));

                    // Write out blocks of data 10 characters long
                    blockCount = 0;
                    currentSite = site;
                    while ((lineCount<50)&&(currentSite<numSites)) {
                        bw.write(theAlphabet.translateIndexToLetter(data[i][currentSite]));
                        currentSite++;
                        lineCount++;
                        blockCount++;
                        if (blockCount==10) {
                            blockCount = 0;
                            bw.write(" ");
                        }
                    }
                    bw.write("\n");
                }
                site += lineCount;
                bw.write("\n");
            }
        }
        catch (AlgorithmException e) {
            throw new IOException("Failed to translate during nexus file output. "+e);
        }


        bw.write(";\n");
        bw.write("end;\n");

        bw.close();
    }

    
    /* Static methods for simulation ----------------------------------------------------------- */

    /** Generate letters from the stationary distribution */
    public static int[] simulateFromStationaryDistribution(int numSites, SubstitutionModel substModel) {

        CategoricalDistribution x = new CategoricalDistribution(substModel.getStationaryDistrib());
        // Loop thru sites
        int[] res = new int[numSites];
        for (int i=0; i<numSites; i++) {
            res[i] = x.sample();
        }
        return res;
    }

    /** Evolve a sequence of letters along a branch */
    public static int[] evolveSequence(int[] startingSeq, double branchLength, SubstitutionModel substModel) {
        double[][] P = substModel.getTransitionMatrix(branchLength);
        // Take cumulative sum across rows.
        int n = substModel.getAlphabet().size();
        CategoricalDistribution[] samplers = new CategoricalDistribution[n];

        for (int i=0; i<n; i++) {
            samplers[i] = new CategoricalDistribution(P[i]);
        }

        int[] seq = new int[startingSeq.length];
        // Loop through sites
        for (int i=0; i<startingSeq.length; i++) {
            seq[i] = samplers[startingSeq[i]].sample();
        }
        return seq;
    }

    /** Evolve a sequence of letters along a branch 
     rateDistib[0] are the relative rates, rateDistib[1] are the corresponding probabilities 
     A rate of 1E-6 or less is treated as zero */
    public static int[] evolveSequenceWithDiscreteRateHeterogeneity(int[] startingSeq, double[][] rateDistrib, double branchLength, SubstitutionModel substModel) {
        
        int numCats = rateDistrib[0].length;
        int[] numSitesPerCategory = new int[numCats];
        CategoricalDistribution cd = new CategoricalDistribution(rateDistrib[1]);
        for (int i=0; i<startingSeq.length; i++) {
            int k= cd.sample();
            numSitesPerCategory[k]++;
        }
        
        int[] seq = new int[startingSeq.length];
        
        int currentSite = 0;
        for (int j=0; j<numCats; j++) {
            
            if (rateDistrib[0][j]<1E-6) {
        
                double[][] P = substModel.getTransitionMatrix(branchLength*rateDistrib[0][j]); // Multiply edge length by relative rate for this category

                int n = substModel.getAlphabet().size();
                CategoricalDistribution[] samplers = new CategoricalDistribution[n];

                for (int i=0; i<n; i++) {
                    samplers[i] = new CategoricalDistribution(P[i]);
                }

                // Loop through sites
                for (int i=0; i<numSitesPerCategory[j]; i++) {
                    seq[currentSite] = samplers[startingSeq[currentSite]].sample();
                    currentSite++;
                }
            }
            else {
                // Loop through sites
                for (int i=0; i<numSitesPerCategory[j]; i++) {
                    seq[currentSite] = startingSeq[currentSite];
                    currentSite++;
                }
            }
            
        }
        return seq;
    }

    /** Simulate alignments when there is no rate heterogeneity. */

    public static Alignment simulateAlignment(int numSites, Tree theTree, SubstitutionModel substModel) {
        Graph.Vertex v = StaticLikelihoodCalculator.findBaseVertex(theTree);
        return simulateAlignment(numSites, theTree, substModel, v);
    }
    public static Alignment simulateAlignment(int numSites, Tree theTree, SubstitutionModel substModel, Graph.Vertex v) {
        cern.jet.random.tdouble.DoubleUniform unif = new cern.jet.random.tdouble.DoubleUniform(simulation.Random.getEngine());
        int[] startingSeq = simulateFromStationaryDistribution(numSites, substModel);
        Alignment theAlignment = new Alignment(theTree.numTaxa(), numSites, substModel.getAlphabet());
        ArrayList taxa = new ArrayList(theTree.getTaxa());
        Collections.sort(taxa);
        try {
            theAlignment.setTaxonNames(taxa);
        }
        catch (AlgorithmError e) {
            System.out.println("Algorithm error setting taxon names when simulating an alignment. This should not be possible. ");
        }
        recursiveSimulateAlignment(v, null, startingSeq, theAlignment, substModel);
        return theAlignment;
    }

    /** Private method used to traverse down tree generating sequences */
    private static void recursiveSimulateAlignment(Graph.Vertex v, Graph.Vertex fromV, int[] startingSeq, Alignment theAlignment, SubstitutionModel substModel) {
        int[] seq = null;
        if (fromV!=null) {
            // Evolve sequence down branch
            seq = evolveSequence(startingSeq, v.getNeighbourDistance(fromV), substModel);
        }
        else {
            seq = startingSeq;
        }
        if (v.degree()==1) {
            // Store in alignment
            try {
                theAlignment.setRow(v.label, seq);
            }
            catch (AlgorithmException e) {
                System.out.println("Error simulating alignment for taxon "+v.label+". This should not be possible. "+e.getMessage());
            }
        }

        // Now loop thru' children
        HashSet children = v.getNeighbours();
        if (fromV!=null) children.remove(fromV);
        Iterator it = children.iterator();
        for (int i=0; i<children.size(); i++) {
            Graph.Vertex w = (Graph.Vertex) it.next();
            recursiveSimulateAlignment(w, v, seq, theAlignment, substModel);
        }
    }

    
    /** Simulate alignments with rate heterogeneity
     rateDistib[0] are the relative rates, rateDistib[1] are the corresponding probabilities*/
    public static Alignment simulateAlignment(int numSites, Tree theTree, double[][] rateDistrib, SubstitutionModel substModel){
        Graph.Vertex v = StaticLikelihoodCalculator.findBaseVertex(theTree);
        return simulateAlignment(numSites, v, theTree, rateDistrib, substModel);
    }
    public static Alignment simulateAlignment(int numSites, Graph.Vertex v, Tree theTree, double[][] rateDistrib, SubstitutionModel substModel) {
        cern.jet.random.tdouble.DoubleUniform unif = new cern.jet.random.tdouble.DoubleUniform(simulation.Random.getEngine());
        int[] startingSeq = simulateFromStationaryDistribution(numSites, substModel);
        Alignment theAlignment = new Alignment(theTree.numTaxa(), numSites, substModel.getAlphabet());
        ArrayList taxa = new ArrayList(theTree.getTaxa());
        Collections.sort(taxa);
        try {
            theAlignment.setTaxonNames(taxa);
        }
        catch (AlgorithmError e) {
            System.out.println("Algorithm error setting taxon names when simulating an alignment. This should not be possible. ");
        }
        recursiveSimulateAlignment(v, null, startingSeq, theAlignment, substModel, rateDistrib);
        return theAlignment;
    }
    
    
    /** Private method used to traverse down tree generating sequences */
    private static void recursiveSimulateAlignment(Graph.Vertex v, Graph.Vertex fromV, int[] startingSeq, Alignment theAlignment, SubstitutionModel substModel, double[][] rateDistrib) {
        int[] seq = null;
        if (fromV!=null) {
            // Evolve sequence down branch
            seq = evolveSequenceWithDiscreteRateHeterogeneity(startingSeq, rateDistrib, v.getNeighbourDistance(fromV), substModel);
        }
        else {
            seq = startingSeq;
        }
        if (v.degree()==1) {
            // Store in alignment
            try {
                theAlignment.setRow(v.label, seq);
            }
            catch (AlgorithmException e) {
                System.out.println("Error simulating alignment for taxon "+v.label+". This should not be possible. "+e.getMessage());
            }
        }

        // Now loop thru' children
        HashSet children = v.getNeighbours();
        if (fromV!=null) children.remove(fromV);
        Iterator it = children.iterator();
        for (int i=0; i<children.size(); i++) {
            Graph.Vertex w = (Graph.Vertex) it.next();
            recursiveSimulateAlignment(w, v, seq, theAlignment, substModel, rateDistrib);
        }
    }


    
    /*------------------------------------------------------------------------*/
    
    /* Legacy code: Simulate alignments with continuous gamma rate heterogeneity */
    /** Method used to traverse down tree generating sequences */
    public static void recursiveSimulateSite(int siteNum, int startingLetter, double relRate, Graph.Vertex v, Graph.Vertex fromV, Alignment theAlignment, SubstitutionModel substModel) {
        int nextLetter;
        if (fromV!=null) {
            // Evolve sequence down branch
            double[][] P = substModel.getTransitionMatrix(v.getNeighbourDistance(fromV)*relRate);
            CategoricalDistribution x = new CategoricalDistribution(P[startingLetter]);
            nextLetter = x.sample();
        }
        else {
            nextLetter = startingLetter;
        }
        if (v.degree()==1) {
            // Store in alignment
            try {
                theAlignment.setLetter(v.label, siteNum, nextLetter);
            }
            catch (AlgorithmException e) {
                System.out.println("Error simulating alignment for taxon "+v.label+". This should not be possible. "+e.getMessage());
            }
        }

        // Now loop thru' children
        HashSet children = v.getNeighbours();
        if (fromV!=null) children.remove(fromV);
        Iterator it = children.iterator();
        for (int i=0; i<children.size(); i++) {
            Graph.Vertex w = (Graph.Vertex) it.next();
            recursiveSimulateSite(siteNum, nextLetter, relRate, w, v, theAlignment, substModel);
        }
    }

    public static Alignment simulateAlignment(Tree theTree, double[] rates, SubstitutionModel substModel){
        Graph.Vertex v = StaticLikelihoodCalculator.findBaseVertex(theTree);
        return simulateAlignment(v, theTree, rates, substModel);
    }
    public static Alignment simulateAlignment(Graph.Vertex v, Tree theTree, double[] rates, SubstitutionModel substModel) {
        Alignment theAlignment = new Alignment(theTree.numTaxa(), rates.length, substModel.getAlphabet());
        ArrayList orderedTaxa = new ArrayList(theTree.getTaxa());
        Collections.sort(orderedTaxa);
        try {
            theAlignment.setTaxonNames(orderedTaxa);
        }
        catch (AlgorithmError e) {
            System.out.println("Algorithm error setting taxon names when simulating an alignment. This should not be possible. ");
        }

        // Loop across sites
        for (int i=0; i<rates.length; i++) {
            simulateAlignmentColumn(v, theAlignment, substModel, i, rates[i]);
        }

        return theAlignment;
    }

    public static void simulateAlignmentColumn(Graph.Vertex v, Alignment theAlignment, SubstitutionModel substModel, int colNum, double rate) {
        int initialLetter = simulateFromStationaryDistribution(1, substModel)[0];
        recursiveSimulateSite(colNum, initialLetter, rate, v, null, theAlignment, substModel);
    }

    
}
