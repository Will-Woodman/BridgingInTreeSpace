/**
    NewickStringUtils
    Static methods for parsing Newick Strings etc

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

/**
 * Static methods for parsing Newick Strings etc
 */

import java.util.HashMap;

public class NewickStringUtils {
    
    public static String removeTrailingBrackets(String theString) {
        // Remove suffix "(" and tail ");" 
        theString = theString.trim();
        if (theString.startsWith("(")) {
            theString = theString.substring(1);
        }
        if (theString.endsWith(";")) {
            theString = theString.substring(0,theString.length()-1);
        }
        theString = theString.trim();
        if (theString.endsWith(")")) {
            theString = theString.substring(0,theString.length()-1);
        }
        return theString;
    }

    /** Splits a Newick string (with opening and closing brackets removed)
     into comma separated substrings. The hashmap returned points from each
     substring to the corresponding branch length. */
    public static HashMap<String,Double> parseToSubstrings(String theString) throws ParsingException {
        HashMap<String,Double> blocksAndLengths = new HashMap<String,Double>();

        int endPos, endOfNumberPos;
        String numString;
        double t;

        while (theString.length()>0) {

            // If leading bracket, find matching close
            if (theString.charAt(0)=='(') {
                endPos = findMatchingBracket(0,theString);
                if (theString.charAt(endPos+1)==':') {
                    endPos = endPos+1;
                }
                else {
                    throw new ParsingException("Close bracket not followed by colon in Newick string","Missing colon");
                }
            }
            else {
                // The current element in the list is a single taxon
                endPos = theString.indexOf(":");
                if (endPos<0) {
                    if (blocksAndLengths.size()==0) {
                        // The entire string is just a taxon label
                        return blocksAndLengths; // An empty map
                    }
                    else {
                        // Error!
                        throw new ParsingException("Missing colon after taxon name in Newick string","Missing colon");
                    }
                }
            }

            String subTreeString = theString.substring(0,endPos);
            theString = theString.substring(endPos+1);

            // Get branch length
            endOfNumberPos = getCloseBracketOrCommaPos(theString);
            if (endOfNumberPos<0) endOfNumberPos = theString.length();
            // Extract the string containing the number
            numString = theString.substring(0,endOfNumberPos);
            // Remove the number from the start of the string
            theString = theString.substring(endOfNumberPos);
            if (theString.startsWith(",")) theString = theString.substring(1);
            theString = theString.trim();

            // Convert the number string into a branch length
            t = 0.0;
            try {
                t = Double.parseDouble(numString);
            }
            catch (NumberFormatException anErr) {
               throw new ParsingException("Error parsing dendrogram string. A branch length for a leaf node was formatted incorrectly.","Bad branch length");
            }

            /* OK. Now got a sub-tree string, a branch length and remainder of string */
            
            // Trim brackets off the sub-tree
            subTreeString = removeTrailingBrackets(subTreeString);
            // Add to the hashmap
            blocksAndLengths.put(subTreeString, new Double(t));


        } // End while loop thru' blocks in the Newick string
        return blocksAndLengths;
    }


    /** Find the matching close bracket in a string */
    public static int findMatchingBracket(int startPos, String theString) throws ParsingException {
        // Check string starts "("
        if (theString.charAt(startPos)!='(') {
            throw new ParsingException("Newick string parser was expecting an open bracket.","Error matching brackets");
        }

        int score = 1;
        int currentPos = startPos+1;
        while (currentPos<theString.length()){
            String s = theString.substring(currentPos);
            int indexOpen = s.indexOf("(");
            int indexClose = s.indexOf(")");
            if ((indexOpen==-1)&&(indexClose==-1))
            {
                throw new ParsingException("Error parsing Newick string. " +
                                                      "Couldn't find a matching bracket.","Bad DND format.");
            }
            // Which came first? "(" or ")"?
            int index;
            if (indexOpen==-1) index = indexClose;
            else if (indexClose == -1) index = indexOpen;
            else index = (indexOpen<indexClose ? indexOpen : indexClose);

            if (s.charAt(index)=='(') {
                score++;
            }
            else {
                score--;
            }

            if (score<0) throw new ParsingException("Error parsing Newick string. " +
                                                      "Unbalanced close brackets","Bad DND format.");

            currentPos = currentPos+index+1;
            if (score==0) return (currentPos-1);
        }
        throw new ParsingException("Error parsing Newick string. Unbalanced close brackets","Bad DND format.");
    }



    /** Get the position of the "," or "(" following a branch length */
    public static int getOpenBracketOrCommaPos(String theString) {
        int indexComma = theString.indexOf(",");
        int indexBracket = theString.indexOf("(");
        if ((indexComma==-1)&&(indexBracket==-1))
        {
            return -1;
        }
        // Which came first? "," or "("?
        // NB: you might not find both of these characters
        int index;
        if (indexComma==-1) index = indexBracket;
        else if (indexBracket == -1) index = indexComma;
        else index = (indexComma<indexBracket ? indexComma : indexBracket);
        return index;
    }

    /** Get the position of the "," or ")" following a branch length */
    public static int getCloseBracketOrCommaPos(String theString) throws ParsingException {
        int indexComma = theString.indexOf(",");
        int indexBracket = theString.indexOf(")");
        if ((indexComma==-1)&&(indexBracket==-1))
        {
            return -1;
        }
        // Which came first? "," or ")"?
        // NB: you might not find both of these characters
        int index;
        if (indexComma==-1) index = indexBracket;
        else if (indexBracket == -1) index = indexComma;
        else index = (indexComma<indexBracket ? indexComma : indexBracket);
        return index;
    }


    /** A class representing a parsing exception */
    public static class ParsingException extends AlgorithmException {

        private String reason;

        public ParsingException() {
            super("An unspecified ParsingException has been thrown...");
            reason = new String();
        }

        public ParsingException(String msg, String r) {
            super(msg);
            reason = new String(r);
        }

        public String getReason() {
            return reason;
        }
    }

    /** Find out whether two Newick strings represent the same tree topology */
    public static boolean compareNewickStringTopologies(String s1, String s2) throws AlgorithmException {
        Tree t1 = new Tree(s1);
        Tree t2 = new Tree(s2);
        java.util.HashSet<Split> splits1 = t1.getSplits();
        java.util.HashSet<Split> splits2 = t2.getSplits();
        if(splits1.size()!=splits2.size()) return false;
        java.util.Iterator<Split> it = splits1.iterator();
        for(int i=0; i<splits1.size(); i++) {
            Split sp1 = it.next();
            java.util.Iterator<Split> it2 = splits2.iterator();
            for(int j=0; j<splits2.size(); j++) {
                Split sp2 = it2.next();
                if(sp1.equals(sp2)) {
                    splits2.remove(sp2);
                    break;
                }
                if(j==splits2.size()-1) return false;
            }
        }
        return true;
    }

        /** Find out whether two Newick strings represent the same tree (ie trees with the same sets
        of splits AND, for a degree two root, trees where the root split is the same) */
    public static boolean compareRootedNewickStringTopologies(String s1, String s2) throws AlgorithmException {
        RootedTree t1 = new RootedTree(s1);
        RootedTree t2 = new RootedTree(s2);
        if(t1.getRoot().degree()!=t2.getRoot().degree()) return false;
        if(t1.getRoot().degree()==2) {
            java.util.Iterator<Graph.Vertex> itV = t1.getRoot().neighboursIterator();
            Split rootSplit1 = t1.getSplit(t1.getRoot().getEdge(itV.next()));
            itV = t2.getRoot().neighboursIterator();
            Split rootSplit2 = t2.getSplit(t2.getRoot().getEdge(itV.next()));
            if(!rootSplit1.equals(rootSplit2)) return false;
        }
        java.util.HashSet<Split> splits1 = t1.getSplits();
        java.util.HashSet<Split> splits2 = t2.getSplits();
        if(splits1.size()!=splits2.size()) return false;
        java.util.Iterator<Split> it = splits1.iterator();
        for(int i=0; i<splits1.size(); i++) {
            Split sp1 = it.next();
            java.util.Iterator<Split> it2 = splits2.iterator();
            for(int j=0; j<splits2.size(); j++) {
                Split sp2 = it2.next();
                if(sp1.equals(sp2)) {
                    splits2.remove(sp2);
                    break;
                }
                if(j==splits2.size()-1) return false;
            }
        }
        return true;
    }


}
