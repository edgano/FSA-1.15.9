
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Michael Smoot and Adam Roberts.
 */


package mad;


import javax.swing.*;
import javax.swing.text.*;
import java.util.*;
import java.awt.event.*;
import java.awt.Dimension;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.font.FontRenderContext;

import java.io.*;
import java.awt.*;
import java.awt.image.*;
import com.sun.image.codec.jpeg.*;
import javax.media.jai.*;
import com.sun.media.jai.codec.*;

public class AlignmentPanel extends JTextPane implements ComponentListener, MouseListener {

	Alignments aligns;
	String[] orderedKeys; 
	HashMap<Integer, Integer> caretToNode;
	int index;
	int numCharsAcross;
	int nameWidth;
	int glyphWidth;
	int glyphHeight;
	int len; 
	Font defaultFont;
	int height;
	int width;
	int os;
	int currPaintIndex;
	
    boolean colored;

    public AlignmentPanel(Alignments aligns, int initWidth, boolean colored) {
		super();
		this.addMouseListener(this);
		this.aligns = aligns;
		this.width = initWidth;
		this.colored = colored;
		String[] seqs = aligns.get(0).getSequences();
		
		nameWidth = 0;
		for (String key : aligns.keys)
			nameWidth = Math.max(nameWidth, key.length());
		nameWidth += 2;	
		

		orderedKeys = new String[aligns.keys.size()];
		int i = 0;
		for (String key : aligns.keys){
			StringBuffer blanks = new StringBuffer("");
			for (int j = 0; j < nameWidth-key.length(); j++)
				blanks.append(" ");
			orderedKeys[i++] = " " + key + blanks;
		}
			
		
		defaultFont = new Font("Courier",Font.BOLD,12);
		glyphWidth=10;
		glyphHeight=10;
		
		index = 0;
		currPaintIndex = -2;

		// Aaaaargh
		if ( System.getProperty("os.name").matches(".*[Mm][Aa][Cc].*") )
			os = -1;
		else
			os = 1;
			
		this.setEditable(false);

		sizeScreen( new Dimension(initWidth,0) );
	}


	private void sizeScreen(Dimension d) {
		width = (int)(d.getWidth()); 

		// find the width of the sequence name
		FontMetrics fm = getFontMetrics(defaultFont);

		glyphWidth = fm.charWidth('W');
		glyphHeight = fm.getHeight();
		
		
		// figure out how many multiple alignment rows there will be	
		int seqAllowedWidth = width - 20 - 2*glyphWidth - nameWidth*glyphWidth - 2*glyphWidth - 20;
		numCharsAcross = seqAllowedWidth/glyphWidth;
		
		if (numCharsAcross > 0){
			currPaintIndex = -2;
			paint();
		}

		revalidate();
	}
	

	public void setIndex(int i) {
		if ( i >= -1 && i < aligns.size() ){
			index = i;
			paint();
		}
		else
			System.out.println("invalid index: " + i);
	}	

	    // The primitive type 'char' in Java acts as an unsigned int taking
	    // values from 0 to 2^16 (65,536).

    public void paint() {
	
	if (index == currPaintIndex)
	    return;
	
	caretToNode = new HashMap<Integer, Integer>();
	//System.out.println("Starting paint: " + index );
	String[] seqs = aligns.get(index).getSequences();
	String[] cols = aligns.get(index).getColors();

		
	BatchDocument doc = new BatchDocument();
	Style style = (Style) doc.getStyle(StyleContext.DEFAULT_STYLE);
	StyleConstants.setFontFamily(style, "Courier");
	StyleConstants.setBold(style, true);
	StyleConstants.setFontSize(style, defaultFont.getSize());
		
	int pos = 0;
	int caretPos = 0;
		
	while (pos < seqs[0].length()){
	
	    for (int i = 0; i < orderedKeys.length; i++ ) {
		String name = orderedKeys[i];
		String seq = seqs[i];
		String colString = cols[i];
			
		if ( seq == null )
		    continue;

		StyleConstants.setForeground(style, Color.black);
		StyleConstants.setBackground(style, Color.white);
		doc.appendBatchString(name, style);
		caretPos += name.length();
		StyleConstants.setForeground(style, Color.white);

		for(int j = pos; j < pos+numCharsAcross && j < seq.length() ; j++){

		    // get color code
		    // see Node::getColors() for a description of the encoding used
		    int colCode = new Integer(colString.charAt(j));

		    if (colored) {
			// if highlighted
			if (colCode == 65535) {
			    StyleConstants.setBackground(style, Color.black);
			}
			// if undefined
			else if (colCode == 256) {
			    StyleConstants.setBackground(style, Color.lightGray);
			}
			// else regular coloration based on accuracy measure
			else{
			    StyleConstants.setBackground(style, getColorFromCode(colCode));
			}
		    }

		    else {
		StyleConstants.setForeground(style, Color.black);
		StyleConstants.setBackground(style, Color.lightGray);

		//			StyleConstants.setBackground(style, Color.lightGray);
		    }

				
		    doc.appendBatchString(seq.substring(j,j+1), style);
					
		    caretToNode.put(caretPos++, j);
		}
		doc.appendBatchLineFeed(style);
		caretPos++;
	    }
	    pos += numCharsAcross;
	    doc.appendBatchLineFeed(style);
	    caretPos++;
	}
		
	try{ doc.processBatchUpdates(0);}
	catch (BadLocationException badLocationException) {System.err.println("Oops");}
	//System.out.println("Setting document: " + index);
	this.setStyledDocument(doc);
	//System.out.println("Finished paint: " + index + "\n");
		
	currPaintIndex = index;
    }
	
    public void setColored(boolean colored){
	this.colored = colored;
    }

    public Color getColorFromCode(int colCode){
	if (colCode >= 256){ //Highlighted
	    return Color.black;
	}
	int r = (int)Math.max(Math.round(-0.015564*(colCode-127)*(colCode-382)),0);
	int g = (int)Math.max(Math.round(-0.015564*(colCode)*(colCode-255)),0);
	int b = (int)Math.max(Math.round(-0.015564*(colCode+127)*(colCode-127)),0);
		
	return new Color(r,g,b);
    }
    
    /// Save the alignment as TIFF.
    public void saveAsTiff(String filename){

	try {

	    // get the iamge
	    BufferedImage saveImage = new BufferedImage (this.getWidth(), this.getHeight(), BufferedImage.TYPE_INT_RGB);
	    // paint it
	    Graphics gc = saveImage.getGraphics();
	    paint(gc);

	    // Store the image in the TIFF format.
	    TIFFEncodeParam encodeParam = new TIFFEncodeParam();
	    //	    encodeParam.setCompression(TIFFEncodeParam.COMPRESSION_DEFLATE);

//	    TIFFField[] extras = new TIFFField[3];
//	    extras[0] = new TIFFField(282,TIFFField.TIFF_RATIONAL, 1, (Object)new long[][] {{300,(long)1},{(long)0 ,(long)0}});
//	    extras[1] = new TIFFField(283,TIFFField.TIFF_RATIONAL, 1, (Object)new long[][] {{300,(long)1},{(long)0 ,(long)0}}); // x
//	    extras[2] = new TIFFField(296, TIFFField.TIFF_SHORT, 1, (Object) new char[] {2}); //2 for inches 	    //y //set resolution unit to inches
//	    encodeParam.setExtraFields(extras);

//	    final int XRES_TAG = 282;
//	    final int YRES_TAG = 283;
//
//	    int[] resolution = { 300, 1};
//	    TIFFField xRes = new TIFFField(XRES_TAG, TIFFField.TIFF_RATIONAL, 1, new int[][] { resolution });
//	    TIFFField yRes = new TIFFField(YRES_TAG, TIFFField.TIFF_RATIONAL, 1, new int[][] { resolution });
//	    encodeParam.setExtraFields(new TIFFField[] { xRes, yRes});

// above are failed attempts to increase the resolution of the resulting image
// -- RKB 10/15/08

	    RenderedOp op = JAI.create("filestore", saveImage, filename, "tiff", encodeParam);

	} catch (Exception e) {
	    System.out.print (e.toString());
	}

    }

    /// Save the alignment as JPEG.
    public void saveAsJPEG(String filename){
	
        File save_file;
 
         // For encoder
        BufferedImage saveImage;
 
       
	try{
	    FileOutputStream file_out = new FileOutputStream(filename);

	    // get the image
	    saveImage = new BufferedImage(this.getWidth(), this.getHeight(),
					   BufferedImage.TYPE_INT_RGB);
	    // paint it
	    Graphics gc = saveImage.getGraphics();
	    paint(gc);

	    float quality = 0.25f;
	    com.sun.image.codec.jpeg.JPEGEncodeParam encodeParam = JPEGCodec.getDefaultJPEGEncodeParam(saveImage);
	    encodeParam.setQuality(quality,false);

	    JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(file_out, encodeParam);
	    encoder.encode(saveImage);
	    System.out.println("Saved frame: " + filename);
	    file_out.close();
	}
	catch(Exception ex)
	    {System.out.print(ex.toString());}
    }
      


	public void iResized(){
		sizeScreen(this.getSize());
	}
	
	public void mouseClicked(MouseEvent e) {
	
		Alignment align = aligns.currAlign;
		int caretPos = this.getCaret().getMark();
		
		if (align.gapsSuppressed || !caretToNode.containsKey(caretPos))
			return;
			
		int nodeI = caretToNode.get(caretPos);

		Node v = align.L[nodeI];
		v.setHighlight(!v.highlighted);
		aligns.currAlign.buildSeqsAndCols(aligns.currAlign.gapsSuppressed, aligns.currAlign.colorScheme);
		currPaintIndex = -2;
		paint();
    }
    
    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}



	
    public void componentResized(ComponentEvent e) {
		sizeScreen( e.getComponent().getSize() );
    }

    public void componentShown(ComponentEvent e) {};
    public void componentHidden(ComponentEvent e) {};
    public void componentMoved(ComponentEvent e) {};
}
