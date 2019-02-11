
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 *  Code for different accuracy measures was written by Robert Bradley.
 */


// to do:


package mad;

import java.util.*;

//import it.unimi.dsi.fastutil.objects.*;

public class Node{
    
    private int[] seqVals;// Indices of position in each sequence contained in node
    private char[] colors;
    public int seqI; 		// If an initial node, what sequence I contain (useful for gap fills).
    public int seqVal;
    public int numSeqs;
    
    public double acc;
    public double sps;
    public double ppv;
    public double cert;
    public double cons;

    public double accDenom;
    public double spsDenom;
    public double ppvDenom;
    public double certDenom;
    public double consDenom;

    public double weightTgf;

    public Node fwdEdge;
    public boolean highlighted = false;

    public int mergeOrder;

    /* 
     * 0 = Accuracy
     * 1 = Sensitivity
     * 2 = Specificity
     * 3 = Certainty
     * 4 = Consistency
     */ 
    public int colorScheme;

    // Used for a merged node
    public boolean isMerge;
    public Node n1;
    public Node n2;
    
    // Used for topological sort (DFS)
    public int finishTime; 	
    public boolean visited;   
    
    public ProbabilityMatrices pm; // allows access to probabilities to calculate accurcy
    
    public Node(int seqI, int seqVal, ProbabilityMatrices pm){
        this.seqVal = seqVal;
        this.seqI = seqI;
        this.pm = pm;
        this.isMerge = false;

	this.acc = -1;
	this.sps = -1;
	this.ppv = -1;
	this.cert = -1;
	this.cons = -1;

	this.accDenom = -1;
	this.spsDenom = -1;
	this.ppvDenom = -1;
	this.certDenom = -1;
	this.consDenom = -1;

	this.weightTgf = -1;
	this.mergeOrder = -1;

     }
         
    public Node(int[] seqVals, int seqI, ProbabilityMatrices pm){
        this.seqVals = seqVals;
        this.seqI = seqI;
        this.pm = pm;
        this.isMerge = false;

	this.acc = -1;
	this.sps = -1;
	this.ppv = -1;
	this.cert = -1;
	this.cons = -1;

	this.accDenom = -1;
	this.spsDenom = -1;
	this.ppvDenom = -1;
	this.certDenom = -1;
	this.consDenom = -1;

	this.weightTgf = -1;
	this.mergeOrder = -1;

     }
    
    public Node(Node n1, Node n2, int mergeOrder){

        this.n1 = n1;
        this.n2 = n2;
        this.pm = n1.pm;
        this.isMerge = true;

	this.acc = -1;
	this.sps = -1;
	this.ppv = -1;
	this.cert = -1;
	this.cons = -1;

	this.accDenom = -1;
	this.spsDenom = -1;
	this.ppvDenom = -1;
	this.certDenom = -1;
	this.consDenom = -1;

	this.weightTgf = calcWeightTgf (n1, n2);
	this.mergeOrder = mergeOrder;

        mergeNodes();
    }
    
  
    private void mergeNodes(){
// 		int[] s1 = n1.seqIndices;
//     	int[] s2 = n2.seqIndices;
//     		
//     	seqIndices = new int[s1.length];
//     	
//     	//Merge columns.
//     	for (int i = 0; i < s1.length; i++) {
//     		if (s1[i] != -1 && s2[i] != -1)
//     			System.out.println("Sanity Check Failed: Sequence has been overwritten.");
//     			
//     		if (s1[i] != -1)
//     			seqIndices[i] = s1[i];
//     		else if (s2[i] != -1)
//     			seqIndices[i] = s2[i];
//     		else
//     			seqIndices[i] = -1;
//     	}
    	
  
  		//Merge highlighting
  		this.highlighted = (n1.highlighted || n2.highlighted);
  		
    	//Merge edges
    	//this.fEdges = new HashSet<Node>();
    	//this.fEdges.addAll(n1.getFwdEdges());
    	//this.fEdges.addAll(n2.getFwdEdges());   	
    }
    
    public void setHighlight(boolean val){
    	if (val == highlighted)
    		return;
    	
    	if (isMerge){
    		n1.setHighlight(val);
    		n2.setHighlight(val);
    	}
    	
    	highlighted = val;
    	calcColorAccuracy();
    }
    
    /// Calculate (updated) tgf weight for the merged column.
    private double calcWeightTgf(Node n1, Node n2){
	
	// initialize
	weightTgf = -1;
	
	int x,y;
	double prob;

	double pMatch = 0.;
	double pGap = 0.;

	int[] n1SeqIndices = n1.getSeqIndices();
	int[] n2SeqIndices = n2.getSeqIndices();

	for(int seq1=0; seq1 < pm.numSeqs; seq1++){

	    x = n1SeqIndices[seq1];

	    for (int seq2 = 0; seq2 < pm.numSeqs; seq2++){

		// if we don't have probability information for the sequence pair seq1, seq2 then skip it
		if ((seq1 == seq2) || !pm.matrixExists(seq1, seq2))
		    continue;
		
		y = n2SeqIndices[seq2];

		// only look at pairs of characters
		if (x == -1 || y == -1)
		    continue;

		// increment pMatch
		pMatch += pm.getElement(seq1,seq2,x,y);
		// increment pGap
		pGap += pm.getElement(seq1,seq2,x,-1);
		pGap += pm.getElement(seq1,seq2,-1,y);

	    }

	}

//	System.err.print ("n1: "); n1.output();
//	System.err.print ("n2: "); n2.output();
//	System.err.println ("  => pMatch = " + pMatch + "; pGap = " + pGap);

	return 2 * pMatch / pGap;

    }

    /// Calculate accuracy metrics for the column and store desired coloration.
    /*
     * Coloration is encoded as an array of chars.
     * \see getColors()
     */
    private void calcColorAccuracy(){
	
	// initialize
	acc = 0;
	sps = 0;
	ppv = 0;
	cert = 0;
	cons = 0;

	accDenom = 0;
	spsDenom = 0;
	ppvDenom = 0;
	certDenom = 0;
	consDenom = 0;

	// coloration for each character or gap in column
	// color scheme chosen based on this.colorScheme
	Double[] colorAccs = new Double[pm.numSeqs];

	double num;    // per-character numerator
	double denom;  // per-character denominator (must be a double b/c takes non-integer values for SPS, Certainty and Consistency)

	int x, y;      // sequence positions
	double prob;   // probability (match or gap)

	boolean hasMatch;  // are there at least two characters in the column?

	int[] seqIndices = getSeqIndices(); // map from sequences to sequence positions in column
	for(int seq1=0; seq1 < pm.numSeqs; seq1++){

	    // initialize values
	    num = 0;
	    denom = 0;
	    colorAccs[seq1] = (Double) 0.0;
	    hasMatch = false;
	
	    x = seqIndices[seq1];

	    for (int seq2 = 0; seq2 < pm.numSeqs; seq2++){

		// if we don't have probability information for the sequence pair seq1, seq2 then skip it
		if ((seq1 == seq2) || !pm.matrixExists(seq1, seq2))
		    continue;
		
		y = seqIndices[seq2];

		int multiplier = 0;
		if (x != -1 && y != -1) {      // neither are gaps
		    multiplier = 2;
		    hasMatch = true;
		}
		else if (x != -1 || y != -1) { // exactly 1 is a gap
		    multiplier = 1;
		}
		else {
		    continue;
		}
		
		prob = pm.getElement(seq1,seq2,x,y);

		// everything contributes to accuracy
		acc += multiplier * prob;
		accDenom += multiplier;
		if (colorScheme == 0) { // Accuracy
		    num += multiplier * prob;
		    denom += multiplier;
		}

		// certainty is colored for both character and gaps,
		// but we only add to the full certainty calculation if 
		// the currect position isn't a gap (otherwise we overcount the gap information)
		if (colorScheme == 3) { // Certainty
		    num += multiplier * prob;
		    if (x != -1) // if x is not a gap, then use altmax_y P(x ~ y)
			denom += multiplier * pm.getAltMaxProb(seq1,seq2,x,y);
		    else         // if x is a gap, then use altmax_x P(x ~ y)
			denom += multiplier * pm.getAltMaxProb(seq2,seq1,y,x);
		}

		// and for consistency
		cons += multiplier * prob;
		if (x != -1) // if x is not a gap, then use max_y P(x ~ y)
		    consDenom += multiplier * pm.getMaxProb(seq1,seq2,x);
		else         // if x is a gap, then use max_x P(x ~ y)
		    consDenom += multiplier * pm.getMaxProb(seq2,seq1,y);
		if (colorScheme == 4) { // Consistency
		    num += multiplier * prob;
		    if (x != -1) // if x is not a gap, then use max_y P(x ~ y)
			denom += multiplier * pm.getMaxProb(seq1,seq2,x);
		    else         // if x is a gap, then use max_x P(x ~ y)
			denom += multiplier * pm.getMaxProb(seq2,seq1,y);
		}

		// in contrast, SPS and PPV aren't defined for gaps
		if (x != -1) {

		    // the numerator for SPS increases only if we're looking at aligned characters (y not a gap);
		    // the denominator increases regardless
		    //   (this is guaranteed by the (multiplier - 1) factor)
		    sps += (multiplier - 1) * prob;
		    spsDenom += pm.getSumProb(seq1,seq2,x);
		    if (colorScheme == 1) { // SPS
			num += (multiplier - 1) * prob;
			denom += pm.getSumProb(seq1,seq2,x);
		    }

		    // only pairs of aligned characters contribute to PPV
		    // (both x and y ungapped)
		    if (y != -1) {
			ppv += multiplier * prob;
			ppvDenom += multiplier;
			if (colorScheme == 2) { // PPV
			    num += multiplier * prob;
			    denom += multiplier;
			}
		    }

		    // certainty
		    cert += multiplier * prob;
		    if (x != -1) // if x is not a gap, then use altmax_y P(x ~ y)
			certDenom += multiplier * pm.getAltMaxProb(seq1,seq2,x,y);
		    else         // if x is a gap, then use altmax_x P(x ~ y)
			certDenom += multiplier * pm.getAltMaxProb(seq2,seq1,y,x);

		}

	    }

	    // store value for this character or gap if it exists
	    // (if the denominator is defined)
	    if (denom > 0.0001) {

		// sps and ppv are only defined if there is at least 1 match
		if (colorScheme == 1 || colorScheme == 2) {
		    if (hasMatch)
			colorAccs[seq1] = (Double) SparseMatrix.saneProb (num / denom);
		    else
			colorAccs[seq1] = (Double) (-1.); // store (temporary) negative value
		}

		// certainty takes values in interval [0, 1/0.01], so it needs to be mapped to [0,1]
		else if (colorScheme == 3)
		    colorAccs[seq1] = (Double) SparseMatrix.saneProb (Node.normalizedCertainty (num / denom));

		// otherwise just store the normalized value
		else
		    colorAccs[seq1] = (Double) SparseMatrix.saneProb (num / denom);

	    }
	    // if undefined, then store (temporary) negative value
	    else {
		colorAccs[seq1] = (Double) (-1.);
	    }

	}//Outer loop

	// Generate color code array
	colors = new char[seqIndices.length];
	for (int i = 0; i < colors.length; i++){

	    double c;

	    // if highlighted
	    if (highlighted) {
		c = 65535.;
	    }
	    // if undefined
	    else if (colorAccs[i].doubleValue() < 0.) {
		c = 256.;
	    }
	    // else encode as integer in interval 0...255
	    else {
		c = Math.rint(colorAccs[i].doubleValue()*255);
	    }

	    // now store the value as a char
	    // remember that 'char' acts as an unsigned short
	    colors[i] = (char)c;
	}

    }
  
    /// Get unnormalized accuracy for column.
    /*
     * Re-calculate if necessary.
     */
    public double getAcc(){
    	if (acc < 0)
	    calcColorAccuracy();
	return acc;
    }

    /// Get normalization for accuracy for column.
    /*
     * Re-calculate if necessary.
     */
    public double getAccDenom(){
    	if (accDenom < 0)
	    calcColorAccuracy();
	return accDenom;
    }

    /// Get unnormalized SPS for column.
    /*
     * Re-calculate if necessary.
     */
    public double getSPS(){
    	if (sps < 0)
	    calcColorAccuracy();
	return sps;
    }

    /// Get normalization for SPS for column.
    /*
     * Re-calculate if necessary.
     */
    public double getSPSDenom(){
    	if (spsDenom < 0)
	    calcColorAccuracy();
	return spsDenom;
    }

    /// Get unnormalized PPV for column.
    /*
     * Re-calculate if necessary.
     */
    public double getPPV(){
    	if (ppv < 0)
	    calcColorAccuracy();
	return ppv;
    }

    /// Get normalization for PPV for column.
    /*
     * Re-calculate if necessary.
     */
    public double getPPVDenom(){
    	if (ppvDenom < 0)
	    calcColorAccuracy();
	return ppvDenom;
    }
  
    /// Get unnormalized certainty for column.
    /*
     * Re-calculate if necessary.
     */
    public double getCert(){
    	if (cert < 0)
	    calcColorAccuracy();
	return cert;
    }

    /// Get normalization for certainty for column.
    /*
     * Re-calculate if necessary.
     */
    public double getCertDenom(){
    	if (certDenom < 0)
	    calcColorAccuracy();
	return certDenom;
    }

    /// Get unnormalized consistency for column.
    /*
     * Re-calculate if necessary.
     */
    public double getCons(){
    	if (cons < 0)
	    calcColorAccuracy();
	return cons;
    }

    /// Get normalization for consistency for column.
    /*
     * Re-calculate if necessary.
     */
    public double getConsDenom(){
    	if (consDenom < 0)
	    calcColorAccuracy();
	return consDenom;
    }
  
    /// Get the tgf weight of the column.
    public double getWeightTgf(){
	if (weightTgf < 0)
	    calcColorAccuracy();
	return weightTgf;
    }

    public int[] getSeqIndices(){
	if (seqVals != null)
	    return seqVals;
	int[] seqIndices = new int[pm.numSeqs];
	Arrays.fill(seqIndices, -1);
	return addSeqIndices(seqIndices);
    }
  	
    // Recursive method that adds all of the children's sequence indices to the array.
    public int[] addSeqIndices(int[] seqIndices){
	if (!isMerge){
	    seqIndices[seqI] = seqVal;
	    return seqIndices;
	}
   		
	seqIndices = n1.addSeqIndices(seqIndices);
	seqIndices = n2.addSeqIndices(seqIndices);
   		
	return seqIndices;
    }


    /// Update color scheme for node.
    /*
     * Updates coloration.
     * It is the responsibility of the parent object
     * to call this function whenever the color scheme is updated globally.
     */
    public void setColorScheme(int scheme){
	this.colorScheme = scheme;
	calcColorAccuracy();
    }

    /// Return an array of chars, each char representing a color.
    /*
     * Chars act as unsigned shorts, taking values from 0 to 2^16.
     * Colors are encoded as:
     * 0...255 coloration for accuracy (=> 256 possible colors in total)
     * 256 indicates highlighting
     * 65,535 indicates undefined
     *  (this value is 2^16 - 1)
     */
    public char[] getColors(){
	if (colors == null)
	    calcColorAccuracy();
	return colors;
    }
   	
    public Set<Node> getFwdEdges(){
	Set<Node> fwdEdges = new HashSet<Node>();
	return addFwdEdges(fwdEdges);
    }	
   	
    // Recursive method that adds all of the children's edges to the set.
    public Set<Node> addFwdEdges(Set<Node> fwdEdges){
   		
	if (isMerge && n1 != null){ // Alternate alignments have merges that w/o any child nodes
	    n1.addFwdEdges(fwdEdges);
	    n2.addFwdEdges(fwdEdges);
	}
	else if (fwdEdge != null){
	    fwdEdges.add(fwdEdge);
	}
   		
	return fwdEdges;
    }
    
    public void output(){
	StringBuffer s = new StringBuffer("");
	s.append(finishTime + ": ");
	int[] seqIndices = getSeqIndices();
	for (int i = 0; i < seqIndices.length; i++)
	    s.append("(" + i + ", " + seqIndices[i] + ") ~ ");
    	
    	System.out.println(s);
    }
    
    /// Normalize a certainty value to interval [0,1].
    /*
     * @param c certainty in range [0, 1 / 0.01]
     * Uses a logarithmic transform.
     */
    static double normalizedCertainty (double c) {

	// if there was a better alternative, then return 0
	if (c < 1.0)
	    return 0;
	// if it was 10 times better than the best alternative, then return 1
	else if (c >= 5.0)
	    return 1;
	// else perform a logarithmic transform
	else
	    return Math.log (c) / Math.log (5.0);

    }

}
