
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 *  Code for different accuracy measures was written by Robert Bradley.
 */


package mad;

import java.util.Scanner;
import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.net.*;

/**
 * Class Name:	Alignments
 *
 * Description: Acts as a wrapper for the DAG to communicate with viewer as if it held a list 
 * 				of the different alignments.  Also parses and stores data from the log file. 
 **/

public class Alignments{
	
	AlignDAG aDag;
	Alignment currAlign;
	Alignment altAlign = null;
	ProbabilityMatrices probMatrices;
	int index;
	int size;
	int[] seqLengths;
	String[] sequences;

	public Double[] alignAccs;
	public Double[] alignSPSs;
	public Double[] alignPPVs;
	public Double[] alignCerts;
	public Double[] alignConss;

        public Double[] alignGFs;

	LinkedList<String> keys;
	int[][] merges; // Merge is indexed so that the current merge will take you to the next alignment.
	boolean suppressGaps;
        int colorScheme;
	
	public Alignments(String seqPath, String guiPath, String probPath){
		suppressGaps = true;
		colorScheme = 0;
		loadData(seqPath, guiPath, probPath);
		makeAccArrays();
	}
	
	public Alignments(String seqPath, String guiPath, String probPath, String altAlignPath){
		this(seqPath, guiPath, probPath);
		loadAltAlign(altAlignPath);
	}

	public int size(){
		return this.size;
	}
	
	private void loadAltAlign(String altAlignPath){
		String[] altSeqs = new String[keys.size()];
		String line;
		BufferedReader br;
		
		try{
			if (altAlignPath.startsWith("http://")){
				br = new BufferedReader( new InputStreamReader( (new URL(altAlignPath)).openStream() ));	
			}
			else{
				br = new BufferedReader( new FileReader(altAlignPath) );	
			}
			line = br.readLine();
			while (line != null && line.length() > 0){
				String key = line.substring(1);
				int i = keys.indexOf(key);
				
				if (i == -1){
					System.err.println("FATAL ERROR: Alternate alignment key '"+ key +"' does not match FSA alignment.");
					System.exit(-1);
				}
				
				StringBuffer seq = new StringBuffer("");
				while ((line = br.readLine()) != null && line.length() > 0 && line.charAt(0) != '>')
					seq.append(line);
				altSeqs[i] = seq.toString();
			}
		} catch (IOException ioe){
			ioe.printStackTrace();
			System.exit(-1);
		}
				
		AlignDAG altDag = new AlignDAG(keys, sequences, altSeqs, probMatrices);
		this.altAlign = new Alignment(keys, sequences, altDag, suppressGaps, colorScheme);
	}

	private String firstWord(String s) {
		int space_pos = s.indexOf(' ');
		if (space_pos == -1) {
			return s;
		} else {
			return s.substring(0, space_pos);
		}
	}
	
	private void loadData(String seqPath, String guiPath, String probPath){

		keys = new LinkedList<String>();
		LinkedList<String> seqs = new LinkedList<String>();
		LinkedList<int[]> merges = new LinkedList<int[]>();
		HashMap<String, Integer> initDAG = new HashMap<String, Integer>();
		String line;
		BufferedReader br;
		Scanner scan;
		
		
		String x;
		String y;

		int i=0;
		
		try{
			if (seqPath.startsWith("http://")){
				System.out.println("WEB");
				br = new BufferedReader( new InputStreamReader( (new URL(seqPath)).openStream() ));	
			}
			else{
				br = new BufferedReader( new FileReader(seqPath) );
			}
			String key = null;
			StringBuffer seq = new StringBuffer("");
			while (true) {
				line = br.readLine();

				// Check for end of previous sequence
				if ((line == null || line.startsWith(">")) && key != null) {
					// Only use the first word of title line as key
					keys.add(firstWord(key.trim()));
					String sequence = seq.toString();
					sequence = sequence.replaceAll("\\s", "");   // remove whitespace
					sequence = sequence.replaceAll("[-_.]", ""); // remove gaps
					seqs.add(sequence);
				}

				// Process next line
				if (line == null) {
					// end of file
					break;
				} else if (line.startsWith(">")) {
					// start a new sequence
					key = line.substring(1);
					seq.setLength(0);
				} else {
					// add to current sequence
					seq.append(line);
				}
			}
			
		} catch (IOException ioe){
			ioe.printStackTrace();
			System.exit(-1);
		} 
		
				
		try{
			if (guiPath.startsWith("http://")){
				br = new BufferedReader( new InputStreamReader( (new URL(guiPath)).openStream() ));	
			}
			else{
				br = new BufferedReader( new FileReader(guiPath) );
			}
			
			// Skip blank lines and comments.		
			while ((line = br.readLine()).length() == 0 || line.charAt(0) == ';')
				continue;
			
			// Read nodes.
			while(line.length() > 0 ){
				if (line.charAt(0) == ';'){ // Skip comments
					line = br.readLine();
					continue;
				}
				
				scan = new Scanner(line);
				//				scan.findInLine("(\\d+):\\s[(](\\d+),\\s(\\d+)[)]\\s=>\\s(\\d+.\\d+|\\d+)");
				scan.findInLine("(\\d+):\\s[(](\\d+),\\s(\\d+)[)]");
		
				
     			MatchResult result = scan.match();
     			i = new Integer(result.group(1));
				x = result.group(2);
				y = result.group(3);
				initDAG.put(x + "." + y, i); // preserve 0-based indexing
				line = br.readLine();

			}
			
	
			//Skip blank lines and comments
			while ((line = br.readLine()).length() == 0 || line.charAt(0) == ';')
				continue;
			
			// Read merges
			while (line != null){
				scan = new Scanner(line);
				scan.findInLine("(\\d+)\\s->\\s(\\d+)");
				MatchResult result = scan.match();
				int[] merge = {new Integer(result.group(1)), new Integer(result.group(2))};
				merges.add(merge);
				line = br.readLine();
			}
				
		} catch (IOException ioe){
			ioe.printStackTrace();
			System.exit(-1);
		}
				
		// Read probabilities and set up matrix
		
		seqLengths = new int[seqs.size()];
		Iterator<String> iter = seqs.iterator();
		i=0;
		while(iter.hasNext())
			seqLengths[i++] = iter.next().length();
		
		probMatrices = new ProbabilityMatrices(seqLengths);
		
		int seq1, seq2;
		int i1, i2;
		double prob;
		
		
		try{
			if (probPath.startsWith("http://")){
				br = new BufferedReader( new InputStreamReader( (new URL(probPath)).openStream() ));	
			}
			else{
				br = new BufferedReader( new FileReader(probPath) );
			}
			while ((line = br.readLine()) != null){
				if ( line.length() == 0 || line.charAt(0) == ';') // Skip comments and blanks
					continue;
				scan = new Scanner(line);
				// cover case of scientific notation
				scan.findInLine("[(](\\d+|-\\d+),\\s(\\d+|-\\d+)[)]\\s~\\s[(](\\d+|-\\d+),\\s(\\d+|-\\d+)[)]\\s=>\\s(\\d.\\d+|\\d+)\\*?[e|E]?(-\\d+)?");
				MatchResult result = scan.match();
			
				seq1 = new Integer(result.group(1)).intValue();
				i1 = new Integer(result.group(2)).intValue();
				seq2 = new Integer(result.group(3)).intValue();
				i2 = new Integer(result.group(4)).intValue();
				prob = new Double(result.group(5)).doubleValue();
				if (result.group(6) != null) {
				    Integer exponent = new Integer((String) result.group(6));
				    prob *= Math.pow(10, (double) exponent);
				}
				probMatrices.addElement(seq1, seq2, i1, i2, prob); // preserve 0-based indexing
				
				//double testProb = probMatrices.getElement(seq1, seq2, i1, i2);
				//if (testProb != prob){
				//	System.err.println("(" + seq1 + ", " + i1 + ") ~ (" + seq2 + ", " + i2 + ") => " + prob + "/" + testProb);
				//	System.exit(-1);
				//}
			}
		} catch (IOException ioe){
			ioe.printStackTrace();
			System.exit(-1);
		}
		
		this.index = 0;
		this.merges = merges.toArray(new int[0][0]);
		this.size = this.merges.length+1;
		this.sequences = seqs.toArray(new String[0]);
		this.aDag = new AlignDAG(keys, sequences, initDAG, probMatrices);
	}
	

    /// Calculate and store Accuracy, SPS, PPV, Certainty and Consistency for all alignments.
    private void makeAccArrays(){
	alignAccs = new Double[merges.length + 1];
	alignSPSs = new Double[merges.length + 1];
	alignPPVs = new Double[merges.length + 1];
	alignCerts = new Double[merges.length + 1];
	alignConss = new Double[merges.length + 1];
	alignGFs = new Double[merges.length + 1];

	alignAccs[0] = (Double) aDag.getAlignAccNorm();
	alignSPSs[0] = (Double) aDag.getAlignSPSNorm();
	alignPPVs[0] = (Double) aDag.getAlignPPVNorm();
	alignCerts[0] = (Double) aDag.getAlignCertNorm();
	alignConss[0] = (Double) aDag.getAlignConsNorm();
	alignGFs[0] = (Double) aDag.getGapFactor();
	
	for (int i = 0; i < merges.length; i ++){
	    next();
	    alignAccs[i+1] = (Double) aDag.getAlignAccNorm();
	    alignSPSs[i+1] = (Double) aDag.getAlignSPSNorm();
	    alignPPVs[i+1] = (Double) aDag.getAlignPPVNorm();
	    alignCerts[i+1] = (Double) aDag.getAlignCertNorm();
	    alignConss[i+1] = (Double) aDag.getAlignConsNorm();
	    alignGFs[i+1] = (Double) aDag.getGapFactor();
	    //	    System.err.println ("i = " + i + "; acc = " + alignAccs[i+1] + "; sps = " + alignSPSs[i+1] + "; ppv = " + alignPPVs[i+1] + "; cert = " + alignCerts[i+1]);
	}
    }
	
	public Alignment get(int index){
				
		if (index == -1){ // -1 should only be request when an alternate alignment is available
						  // otherwise the button is disabled
			return this.altAlign;
		}
		
		if (this.index == index && currAlign != null)
			return this.currAlign;

		while (this.index < index)
			next();
			
		while (this.index > index)
			previous();

		currAlign = new Alignment(keys, sequences, aDag, suppressGaps, colorScheme);
		System.gc();
		return currAlign;
	}
	
	// Merges with the current index and then increments it.
	private void next(){
		aDag.merge(merges[index][0], merges[index][1]);
		index++;
	}
	
	// Decrements the index and uses the last merge as a demerge.
        // Tell AlignDAG::demerge about the next-to-last merge as well to lookup the previous gap factor.
	private void previous(){
		index --;
		if (index > 0)
		    aDag.demerge(merges[index][0], merges[index-1][0]);
		else
		    aDag.demerge(merges[index][0], -1);
	}
	
	public void setGapSuppress(boolean val){
		suppressGaps = val;
		currAlign = null;
		if (altAlign != null && altAlign.gapsSuppressed != suppressGaps)
		    altAlign.buildSeqsAndCols(suppressGaps, colorScheme);
	}

    /// Set color scheme.
    public void setColorScheme(int scheme){

	aDag.setColorScheme (scheme);

	this.colorScheme = scheme;
	currAlign = null;
	if (altAlign != null && altAlign.colorScheme != scheme)
	    altAlign.buildSeqsAndCols(suppressGaps, scheme);
    }

}
