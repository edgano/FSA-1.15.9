
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 *  Code for different accuracy measures was written by Robert Bradley.
 */

package mad;

import java.util.*;

public class Alignment {

    Node[] L; // nodes in topologically sorted order
    String[] seqs;
    String[] cols;
    List<String> keys;
    AlignDAG aDag;
    String[] fullSequences;
    boolean gapsSuppressed;
    int colorScheme;

    public Alignment(List<String> keys, String[] fullSequences, AlignDAG aDag, boolean suppressGaps, int colorScheme) {
	this.keys = keys;
	this.fullSequences = fullSequences;
	this.aDag = aDag;
	this.colorScheme = colorScheme;
	buildSeqsAndCols(suppressGaps, colorScheme);
    }
	

    public void buildSeqsAndCols(boolean suppressGaps, int scheme){

	// update color scheme
	this.colorScheme = scheme;

       	int N = keys.size();
       	
    	// Initialize empty StringBuffer array and sanity check
    	StringBuffer[] cols = new StringBuffer[N];
    	StringBuffer[] seqs = new StringBuffer[N];
	int[] sanityCheck = new int[N];
    	for (int i = 0; i < N; i++){
	    seqs[i] = new StringBuffer();
	    cols[i] = new StringBuffer();
	    sanityCheck[i] = 0;
    	}
    	
    	StringBuffer color = new StringBuffer();
    	char[] colors;
    	int[] seqIndices;
    	Node v;
    	
	// update color scheme
	aDag.setColorScheme (colorScheme);
    	L = aDag.topoSort();
    	
    	//Gap suppression variables
    	int [] gapCounts = new int[1];
    	int gapCountDown = 0;
    	int lastGapCount = 0;
    	int[] colorSum = new int[N];
	if (suppressGaps)
	    gapCounts = getGapCounts(L);
    	
    	
    	// Hide consecutive graps of length greater than 3.
       	for(int i = 0; i < L.length; i++){
        	          
            v = L[i];
            
            int[] sI = v.getSeqIndices();
            for (int z = 0; z<N; z++){
				if (sI[z] != -1){
					if (sI[z] != sanityCheck[z])
						System.err.println("ERROR: " + z + " " + sI[z] + " " + sanityCheck[z]);
					sanityCheck[z] ++;
				}		
            }

			//Handle gap suppression...
            if (!suppressGaps){}
            
            else if (gapCountDown-- == 1){
            	colorSum = sumArrays(colorSum, v.getColors());
            	String[][] gapFill = getGapFill(lastGapCount, colorSum, v.seqI);
            	for (int j = 0; j < N; j++){
       				seqs[j].append(gapFill[j][0]);
            		cols[j].append(gapFill[j][1]);
            	}
            	continue;
            }
    
            else if (gapCountDown+1 > 0){ // Adding 1 because we already decremented
            	colorSum = sumArrays(colorSum, v.getColors());
            	continue;
            }  
    
    		else if ((lastGapCount = gapCounts[i]) > 3){
            	gapCountDown = lastGapCount - 1;
            	colorSum = new int[N];
            	for (int x = 0; x < N; x++)
            		colorSum[x] = 0;
            	colorSum = sumArrays(colorSum, v.getColors());
            	continue;
            }
            //...done handling gap suppression.
            
                     
            colors = v.getColors();
            seqIndices = v.getSeqIndices();
            
            for (int j = 0; j < N; j++){
            	if (seqIndices[j] == -1)
            		seqs[j].append('-');
            	else
            		seqs[j].append(fullSequences[j].charAt(seqIndices[j])); 
            	cols[j].append(colors[j]);
            }
        }
        
        // Convert StringBuffers to Strings
        this.cols = new String[N];
        this.seqs = new String[N];

        for (int i = 0; i < N; i++){
        	this.cols[i] = cols[i].toString();
        	this.seqs[i] = seqs[i].toString();
        }
          
        this.gapsSuppressed = suppressGaps;

    }
    
    // Creates an array specifying how many consecutive gaps follow at each starting point (i).
   	public int[] getGapCounts(Node[] L){
   	
    	int[] gapCounts = new int[L.length];
    	int lastCount = 0;
    	int lastSeqI = -1;
    	Node u;
    	
    	for (int i = L.length-1; i>=0; i--){
    		u = L[i];
    		
    		if (u.isMerge){
    			lastCount = 0;
    			gapCounts[i] = 0;
    		}
    		else{
    			if (lastSeqI != u.seqI)
    				lastCount = 0;
    			gapCounts[i] = ++lastCount;
    			lastSeqI = u.seqI;
    		}
    	}
    	
    	return gapCounts;
    }
    
    private int[] sumArrays(int[] a1, char[] a2){
    	for (int i = 0; i < a1.length; i ++){
    		int c = (int)a2[i];
    		if (c > 256) // Remove highlighting
    			c -= 256;
    		a1[i] += c;
    		
    	}
    	return a1;
    }
    
    // Returns a sequences to suppress a gap of size nGaps.    
    private String[][] getGapFill(int nGaps, int[] colorSum, int seqI){
    	//String nonGapSeq = "\u25B6" + nGaps + "\u25C0";
    	String nonGapSeq = ">" + nGaps + "<";
    	StringBuffer col;
    	StringBuffer seq;
    	
    	int N = keys.size();
    	
    	String[][] results = new String[N][2];
    	
    	char avgColor;
    	
    	for (int i = 0; i < N; i ++){
    		seq = new StringBuffer("");
    		col = new StringBuffer("");
    		avgColor = (char)Math.rint(colorSum[i]/nGaps);
    		for (int j = 0; j < nonGapSeq.length(); j++){
    			seq.append('-');
    			col.append(avgColor);
       		}
    		results[i][0] = seq.toString();
    		results[i][1] = col.toString();
    	}
    	
    	results[seqI][0] = nonGapSeq;
    
    	return results;
    }


	public List<String> getOrderedKeys() {
		return keys;
	}

	public String[] getSequences() {
		return seqs;
	}

	public String[] getColors() {
		return cols;
	}

    /// Get normalized accuracy for alignment.
    /*
     * \return accuracy, or -1 if undefined
     */
    public double getAccNorm() {
	return aDag.getAlignAccNorm();
    }

    /// Get normalized SPS for alignment.
    /*
     * \return SPS, or -1 if undefined
     */
    public double getSPSNorm() {
	return aDag.getAlignSPSNorm();
    }

    /// Get normalized PPV for alignment.
    /*
     * \return PPV, or -1 if undefined
     */
    public double getPPVNorm() {
	return aDag.getAlignPPVNorm();
    }

    /// Get normalized certainty for alignment.
    /*
     * \return certainty, or -1 if undefined
     */
    public double getCertNorm() {
	return aDag.getAlignCertNorm();
    }

    /// Get normalized consistency for alignment.
    /*
     * \return consistency, or -1 if undefined
     */
    public double getConsNorm() {
	return aDag.getAlignConsNorm();
    }

    /// Get the change in normalized accuracy.
    public double getDeltaAccNorm(){
	return aDag.getDeltaAccNorm();
    }

    /// Get the change in normalized SPS.
    public double getDeltaSPSNorm(){
	return aDag.getDeltaSPSNorm();
    }

    /// Get the change in normalized PPV.
    public double getDeltaPPVNorm(){
	return aDag.getDeltaPPVNorm();
    }

    /// Get the change in normalized certainty.
    public double getDeltaCertNorm(){
	return aDag.getDeltaCertNorm();
    }

    /// Get the change in normalized consistency.
    public double getDeltaConsNorm(){
	return aDag.getDeltaConsNorm();
    }

    /// Get the implicit gap factor.
    /*
     * \return -1 if Infinity, implicit current gap factor otherwise
     */
    public double getGapFactor(){
	return aDag.getGapFactor();
    }

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Weight      = ");
		sb.append((new Double(getAccNorm())).toString());
		sb.append("\n");
		int i = 0;
		for ( String key : keys ) {
			sb.append(key);
			sb.append("\t");
			sb.append(seqs[i]);
			sb.append("\n");
			i++;
		}
		sb.append("\n");
		sb.append("\n");

		return sb.toString();
	}

	public String toMultiFasta() {
		if (gapsSuppressed) // This is so that the sequences will be saved without gaps suppression
		    buildSeqsAndCols(false, 0); // throwaway color scheme
		
		StringBuffer sb = new StringBuffer();
		int i = 0;
		
		for ( String key : keys ) {
			sb.append(">");
			sb.append(key);
			sb.append("\n");
			sb.append(seqs[i]);
			sb.append("\n");
			i++;
		}
		
		if (gapsSuppressed)
		    buildSeqsAndCols(true, 0); // throwaway color scheme
		
		return sb.toString();
	}
}
