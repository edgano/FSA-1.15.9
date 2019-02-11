
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Michael Smoot and Adam Roberts.
 *  Code for different accuracy measures was written by Robert Bradley.
 */

// to do:


package mad;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.text.NumberFormatter;
import java.beans.*;
import java.util.Map;
import java.util.HashMap;
import java.util.Vector;
import java.util.Iterator;


import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import java.net.URL;


public class MadPanel extends JPanel 
    implements ChangeListener, PropertyChangeListener, ActionListener, ItemListener {

    AlignmentPanel alignPanel;
    JFormattedTextField indexField;
    JLabel legLabel;
    JFormattedTextField accField;

    JLabel accLabel;
    JLabel dAccLabel;

    JLabel gfLabel;
    JFormattedTextField gfField;

    JSlider alignSlider;
    JSlider zoomSlider;
    Alignments aligns;
    JButton nextButton;
    JButton prevButton;
    JButton startStopButton;
    JButton altAlignButton;
    JCheckBox gapSuppressCheck;
    JFrame parent = null;
    
    JCheckBox coloredCheck;
    JComboBox colorSchemes;
    String[] colorSchemesList;
    String[] colorSchemesLabels;
    
    String fsaPath, altPath;
    
    int maxSliderVal;
    Map<Component,Boolean> ok2update;
    Timer timer;
    boolean frozen = true;
    int delay = 100;

    boolean inAlt = false;
	

    public MadPanel(JFrame parent, Alignments aligns, String fsaPath, String altPath){
    	this(aligns, fsaPath, altPath);
    	this.parent = parent;

    }
    public MadPanel(Alignments aligns, String fsaPath, String altPath) {

	this.aligns = aligns;
	this.fsaPath = fsaPath;
	this.altPath = altPath;
		
	maxSliderVal = aligns.size()-1;

	ok2update = new HashMap<Component,Boolean>();

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

	// set up animation timer
	timer = new Timer(delay,this);
	timer.setCoalesce(true);

        //Create the label.
        JLabel sliderLabel = new JLabel("Alignment: ", JLabel.CENTER);
        sliderLabel.setAlignmentX(Component.CENTER_ALIGNMENT);

        //Create the alignment number field 
        java.text.NumberFormat intFormat = java.text.NumberFormat.getIntegerInstance();
        NumberFormatter intFormatter = new NumberFormatter(intFormat);
        intFormatter.setMinimum(new Integer(0));
        intFormatter.setMaximum(new Integer(maxSliderVal));
        indexField = new JFormattedTextField(intFormatter);
        indexField.setColumns(4); //get some space
        indexField.addPropertyChangeListener( this ); 
	handleEnterKeyStroke( indexField );
	ok2update.put(indexField,false);

        accLabel = new JLabel("Accuracy: ", JLabel.CENTER);
       	dAccLabel = new JLabel("\u0394 Accuracy: ", JLabel.CENTER);
       	gfLabel = new JLabel("Gap Factor: ", JLabel.CENTER);

        //Create the alignment acc field 
        java.text.NumberFormat numFormat = java.text.NumberFormat.getNumberInstance();
        NumberFormatter numFormatter = new NumberFormatter(numFormat);
        numFormatter.setMinimum(new Float(0f));
        numFormatter.setMaximum(new Float(Float.MAX_VALUE));
        accField = new JFormattedTextField(numFormatter);
        accField.setColumns(4); //get some space
        accField.addPropertyChangeListener( this ); 
	handleEnterKeyStroke( accField );
	ok2update.put(accField,false);

        //Create the alignment acc field 
        java.text.NumberFormat numFormat2 = java.text.NumberFormat.getNumberInstance();
        NumberFormatter numFormatter2 = new NumberFormatter(numFormat);
        numFormatter2.setMinimum(new Float(0f));
        numFormatter2.setMaximum(new Float(Float.MAX_VALUE));
        gfField = new JFormattedTextField(numFormatter2);
        gfField.setColumns(4); //get some space
        gfField.addPropertyChangeListener( this ); 
	handleEnterKeyStroke( gfField );
	ok2update.put(gfField,false);

	// create next button
	nextButton = new JButton("Next");
	nextButton.addActionListener(this);

	startStopButton = new JButton("Start Animation");
	startStopButton.addActionListener(this);

	// create prev button
	prevButton = new JButton("Prev");
	prevButton.addActionListener(this);

	// create alternate alignment button
	altAlignButton = new JButton("Alt Alignment");
	altAlignButton.addActionListener(this);
	altAlignButton.setEnabled(aligns.altAlign != null);

        //Create the sliders.
        alignSlider = new JSlider(JSlider.HORIZONTAL, 0, maxSliderVal, 0);
        alignSlider.addChangeListener(this);
        zoomSlider = new JSlider(JSlider.VERTICAL, 0, 20, 10);
        zoomSlider.addChangeListener(this);

        //Turn on labels at major tick marks.
        alignSlider.setMajorTickSpacing(calcTickSpacing(maxSliderVal));
        alignSlider.setPaintTicks(true);
        alignSlider.setPaintLabels(true);

	//Create gap suppression checkbox
	gapSuppressCheck = new JCheckBox("Gap Suppression", true);
	gapSuppressCheck.addItemListener(this);

        //Create the panel that displays the animation.
        alignPanel = new AlignmentPanel(aligns,450,true);
        alignPanel.setBorder(BorderFactory.createCompoundBorder(
								BorderFactory.createLoweredBevelBorder(),
								BorderFactory.createEmptyBorder(10,10,10,10)));

	// set up scrolling
	JScrollPane alignScroll = new JScrollPane(alignPanel);
	alignScroll.setPreferredSize(new Dimension(500,800));
	alignScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
	alignScroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	alignScroll.addComponentListener(alignPanel);

	// Set up zoompanel
	JLabel plusLabel = new JLabel("+", JLabel.CENTER);
	plusLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
	JLabel minusLabel = new JLabel("-", JLabel.CENTER);
	minusLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
	JPanel zoomPanel = new JPanel();
	zoomPanel.setLayout(new BoxLayout(zoomPanel, BoxLayout.Y_AXIS));
	zoomPanel.add(plusLabel);
	zoomPanel.add(zoomSlider);
	zoomPanel.add(minusLabel);
		
	JPanel alignAndZoom = new JPanel(new BorderLayout());
	alignAndZoom.add(zoomPanel,BorderLayout.LINE_START);
	alignAndZoom.add(alignScroll,BorderLayout.CENTER);

	// create colorin checkbox
	coloredCheck = new JCheckBox("Coloring", true);
	coloredCheck.addItemListener(this);

	// set up color schemes	
	colorSchemesList = new String[5];
	colorSchemesList[0] = "Accuracy";
	colorSchemesList[1] = "Sensitivity";
	colorSchemesList[2] = "Specificity";
	colorSchemesList[3] = "Certainty";
	colorSchemesList[4] = "Consistency";

	colorSchemesLabels = new String[5];
	colorSchemesLabels[0] = "Low Accuracy                   High Accuracy";
	colorSchemesLabels[1] = "Low Sensitivity               High Sensitivity";
	colorSchemesLabels[2] = "Low Specificity               High Specificity";
	colorSchemesLabels[3] = "Low Certainty                   High Certainty";
	colorSchemesLabels[4] = "Low Consistency          High Consistency";

	colorSchemes = new JComboBox(colorSchemesList);
	colorSchemes.setSelectedIndex(0);
	colorSchemes.addActionListener(this);

	//Create legend
	URL imgURL = getClass().getResource("legend.jpg");
	ImageIcon legIcon = new ImageIcon(imgURL);
	legLabel = new JLabel("Low Accuracy                   High Accuracy", legIcon, JLabel.CENTER);
	legLabel.setVerticalTextPosition(JLabel.BOTTOM);
	legLabel.setHorizontalTextPosition(JLabel.CENTER);
	legLabel.setAlignmentX(Component.LEFT_ALIGNMENT);

        // Create a subpanel for the label and text field.
        JPanel labelAndTextField = new JPanel(); //use FlowLayout
	labelAndTextField.add(colorSchemes);
	labelAndTextField.add(Box.createRigidArea(new Dimension(20,0)));
	labelAndTextField.add(legLabel);
	labelAndTextField.add(Box.createRigidArea(new Dimension(20,0)));
	labelAndTextField.add(accLabel);
	labelAndTextField.add(accField);
	labelAndTextField.add(Box.createRigidArea(new Dimension(10,0)));
	labelAndTextField.add(dAccLabel);
	labelAndTextField.add(Box.createRigidArea(new Dimension(20,0)));
	labelAndTextField.add(gfLabel);
	labelAndTextField.add(gfField);
	labelAndTextField.add(Box.createRigidArea(new Dimension(20,0)));
        labelAndTextField.add(sliderLabel);
        labelAndTextField.add(indexField);

	JPanel controlPane = new JPanel();
	controlPane.add(prevButton);
	controlPane.add(startStopButton);
	controlPane.add(nextButton);
	controlPane.add(altAlignButton);
	controlPane.add(coloredCheck);
	controlPane.add(gapSuppressCheck);

        // Put everything together.
	// The order in which these are added CANNOT be changed without
	// also changing the indexing in actionPerformed().
        add(alignAndZoom);
        add(Box.createRigidArea(new Dimension(0, 10)));
        add(labelAndTextField);
        add(alignSlider);
        add(controlPane);
        setBorder(BorderFactory.createEmptyBorder(10,10,10,10));

	update(maxSliderVal);
    }

    // Catch changes to sliders
    public void stateChanged(ChangeEvent e) {
    	
        JSlider source = (JSlider)e.getSource();
	if ( source == alignSlider && (boolean) ok2update.get(alignSlider)){
	    update((int)alignSlider.getValue());
	}
	else if (source == zoomSlider){
	    alignPanel.defaultFont = alignPanel.defaultFont.deriveFont(1 + (float)zoomSlider.getValue());
	    alignPanel.iResized();
	    //repaint();
	}
    }

    // Changes the alignment frame
    private void update(int index) {
	Alignment align = aligns.get(index);
	updatePicture(index);
	ok2update.put(alignSlider, false);
	alignSlider.setValue(index);
	ok2update.put(alignSlider, true);
	indexField.setText(Integer.toString(index));

	// update accuracy field
	accField.setText(getAccFieldString(align));
	
	// update delta accuracy field
	dAccLabel.setText(getDeltaAccFieldString(align, index));

	// update gap factor field
	gfField.setText(getGapFactorFieldString(align, index));
    }

    /// Get value for the accuracy field (varies based on color scheme).
    private String getAccFieldString(Alignment align){

	// displayed value depends on the chosen color scheme
	double accFieldValue = 0;
	switch (align.colorScheme) {
	case 0: // Accuracy
	    accFieldValue = align.getAccNorm();
	    break;
	case 1: // SPS
	    accFieldValue = align.getSPSNorm();
	    break;
	case 2: // PPV
	    accFieldValue = align.getPPVNorm();
	    break;
	case 3: // Certainty
	    accFieldValue = align.getCertNorm();
	    break;
	case 4: // Consistency
	    accFieldValue = align.getConsNorm();
	    break;
	default:
	    System.exit(1);
	}

	// test for undefined value
	// (-1 is returned if undefined)
	if (accFieldValue < 0)
	    return "Undefined";
	else {
	    accFieldValue = Math.rint(accFieldValue*10000)/10000;
	    return Double.toString(accFieldValue);
	}

    }

    /// Get value for the delta accuracy field (varies based on color scheme).
    /*
     * @param index index of the alignment
     * \return 0 if first alignment, N/A if alternate alignment and delta accuracy otherwise
     */
    private String getDeltaAccFieldString(Alignment align, int index){

	// displayed value depends on the chosen color scheme
	double delta = 0;
	switch (align.colorScheme) {
	case 0: // Accuracy
	    delta = align.getDeltaAccNorm();
	    break;
	case 1: // SPS
	    delta = align.getDeltaSPSNorm();
	    break;
	case 2: // PPV
	    delta = align.getDeltaPPVNorm();
	    break;
	case 3: // Certainty
	    delta = align.getDeltaCertNorm();
	    break;
	case 4: // Consistency
	    delta = align.getDeltaConsNorm();
	    break;
	default:
	    System.exit(1);
	}

	// handle case of first alignment
	if (index == 0)
	    delta = 0;

	// now format nicely
	NumberFormat formatter;			
	if (delta < 0)
	    formatter = new DecimalFormat("#0.0000");
	else
	    formatter = new DecimalFormat("  #0.0000");

	String entry = "\u0394 Column " + colorSchemesList[align.colorScheme] + ": " + formatter.format(delta);

	// handle case of alternate alignment
	if (index == -1)
	    entry = "\u0394 Column " + colorSchemesList[align.colorScheme] + ":    N/A     ";

	return entry;
    }

    /// Get value for the gap factor field.
    /*
     * @param index index of the alignment
     * \return Infinity if first alignment, N/A if alternate alignment and implicit gap factor otherwise
     */
    private String getGapFactorFieldString(Alignment align, int index){

	double gf = align.getGapFactor();
	gf = Math.rint(gf * 10)/10;

	String entry;
	// handle case of alternate alignment
	if (index == -1)
	    return "N/A";
	// handle case of big gap factor
	else if (gf > 100 || gf < 0)
	    return ">100";
	else
	    return Double.toString(gf);

    }

    // Switch between alternate alignment and FSA alignment
    private void toggleAltAlign(boolean toAlt){
	if (toAlt){
	    inAlt = true;

	    Alignment align = aligns.get(-1);
	    indexField.setText("ALT");
	    updatePicture(-1);

	    // update accuracy field
	    accField.setText(getAccFieldString(align));

	    // update delta accuracy field
	    dAccLabel.setText(getDeltaAccFieldString(align, -1));

	    // update gap factor field
	    gfField.setText(getGapFactorFieldString(align, -1));

	    altAlignButton.setText("FSA Alignment");
	    if (parent != null)
		parent.setTitle("Multiple Alignment Display: " + altPath + "  (Alternate alignment)");
	} else{
	    inAlt = false;

	    update((int)alignSlider.getValue());
	    altAlignButton.setText("Alt Alignment");
	    if (parent != null)
		parent.setTitle("Multiple Alignment Display: " + fsaPath + "  (FSA alignment)");
	}

	prevButton.setEnabled(!toAlt);
	nextButton.setEnabled(!toAlt);
	startStopButton.setEnabled(!toAlt);
	accField.setEnabled(!toAlt);
	gfField.setEnabled(!toAlt);
	indexField.setEnabled(!toAlt);
	alignSlider.setEnabled(!toAlt);
	
    }

    // Catch changes to text boxes
    public void propertyChange(PropertyChangeEvent e) {
       	if (!"value".equals(e.getPropertyName())) 
	    return;

	// alignment index field
	if ( e.getSource() == indexField && ok2update.get(indexField) ) {
	    Number value = (Number)e.getNewValue();
	    if (value != null) {
		update(value.intValue());
		ok2update.put(indexField,false);
	    }
	}

	// accuracy field
	else if ( e.getSource() == accField && ok2update.get(accField) ) {
	    //	    Float acc = (Float) e.getNewValue();
	    int index = -1;
	    int scheme = colorSchemes.getSelectedIndex();
	    switch (scheme) {
	    case 0:
		index = getIndexFromAcc ((Float)e.getNewValue());
		break;
	    case 1:
		index = getIndexFromSPS ((Float)e.getNewValue());
		break;
	    case 2:
		index = getIndexFromPPV ((Float)e.getNewValue());
		break;
	    case 3:
		index = getIndexFromCert ((Float)e.getNewValue());
		break;
	    case 4:
		index = getIndexFromCons ((Float)e.getNewValue());
		break;
	    default:
		System.exit(1);
	    }

	    // if runs off the left, go to first alignment
	    if (index < 0) {
		update(0);
		ok2update.put(accField, false);
	    }
	    // else go to the requested alignment
	    else if (index >= 0 && index <= maxSliderVal) {
		update(index);
		ok2update.put(accField ,false);
	    } 
	}

	// gap factor field
	else if ( e.getSource() == gfField && (boolean) ok2update.get(gfField) ) {

	    int index = getIndexFromGF ((Float)e.getNewValue());

	    // if runs off the left, go to first alignment
	    if (index < 0) {
		update(0);
		ok2update.put(gfField, false);
	    }
	    // else go to the requested alignment
	    else if (index >= 0 && index <= maxSliderVal) {
		update(index);
		ok2update.put(gfField ,false);
	    } 
	}

    }

    // Catch button clicks and timer
    public void actionPerformed(ActionEvent e) {

	int curr = !inAlt ? (int)alignSlider.getValue() : -1;

	// next alignment
	if ( e.getSource() == nextButton) {
	    if ( curr == maxSliderVal )
		update(0);
	    else
		update(curr+1);
	}

	// timer
	else if (e.getSource() == timer){
	    if (curr == maxSliderVal ){
		timer.stop();
		frozen = true;
		startStopButton.setText("Start Animation");
	    }
	    else
		update(curr+1);		
	}
	// previous alignment
	else if ( e.getSource() == prevButton ) {
	    if ( curr == 0 )
		update(maxSliderVal);
	    else
		update(curr-1);			
	}
	// start/stop animation
	else if ( e.getSource() == startStopButton ) {
	    if ( frozen ) {
		if (curr == maxSliderVal)
		    update(0);
		timer.start();
		startStopButton.setText("Stop Animation");
		altAlignButton.setEnabled(false);
	    } else {
		timer.stop();
		startStopButton.setText("Start Animation");
		if (aligns.altAlign != null)
		    altAlignButton.setEnabled(true);
	    }
	    frozen = !frozen;
	}
	// alternate alignment
	else if ( e.getSource() == altAlignButton ) {
	    toggleAltAlign(altAlignButton.getText() == "Alt Alignment");			
	}
	// color schemes
	else if ( e.getSource() == colorSchemes ){

	    JComboBox cb = (JComboBox)e.getSource();

	    // get the integer index of the color scheme
	    int scheme = cb.getSelectedIndex();

	    // update the color scheme
	    aligns.setColorScheme(scheme);

	    // update the legend label
	    legLabel.setText(colorSchemesLabels[scheme]);

	    // update the accuracy field label
	    accLabel.setText(colorSchemesList[scheme] + ": ");

	    // to do: this needs to check whether in the alternate alignment;
	    // currently uses curr, which is always the FSA alignment
	    // (not updated properly for alternate alignment)

	    // update accuracy field value
	    // NB aligns.get(curr) returns the current Alignment object
	    accField.setText(getAccFieldString(aligns.get(curr)));

	    // update the delta accuracy label
	    dAccLabel.setText(getDeltaAccFieldString(aligns.get(curr), curr));

	    // update the gap factor label
	    gfField.setText(getGapFactorFieldString(aligns.get(curr), curr));

	    // redraw the panel
	    alignPanel.iResized();
	}

    }

    // Catch check box change
    public void itemStateChanged(ItemEvent e) {
	
	Object source = e.getItemSelectable();

	// colored
	if (source == coloredCheck){
	    if (e.getStateChange() == ItemEvent.DESELECTED) 
		alignPanel.setColored(false);
	    else
		alignPanel.setColored(true);
	    alignPanel.iResized();
	}

	// gap suppression
	else if (source == gapSuppressCheck){
	    if (e.getStateChange() == ItemEvent.DESELECTED) 
		aligns.setGapSuppress(false);
	    else
		aligns.setGapSuppress(true);
	    alignPanel.iResized();
	    //repaint();
	}
	    
    }


	public void alignmentToMov(){
		
		// Get save path for movie file.
		FileDialog chooseFile = new FileDialog(new Frame(), "Export as Quicktime Movie", FileDialog.SAVE);
		chooseFile.setVisible(true);
		String movFile = chooseFile.getFile();
		if (movFile == null)
			return;
		if (!movFile.endsWith(".mov"))
			movFile = movFile.concat(".mov");
		String dirPath = chooseFile.getDirectory();
		String movPath = dirPath + movFile;
		
		Vector<String> files = new Vector<String>();
		String fName;
			
		
		// Ceate JPEGs
		for (int i = 0; i < aligns.size(); i++){
			update(i);
			fName = String.format(dirPath + "/%08d.jpg",i);
			files.add(fName);
			alignPanel.saveAsJPEG(fName);
		}
		
		// Create movie
		JpegImagesToMovie imageToMovie = new JpegImagesToMovie();
		imageToMovie.doIt2(alignPanel.getWidth(), alignPanel.getHeight(), 5, files, movPath);

		// Delete temporary JPEGs
		Iterator iter = files.iterator();
		while (iter.hasNext()){
			fName = (String) iter.next();
			(new File(fName)).delete();
		}
	}
	
    protected void updatePicture(int frameNum) {
	alignPanel.setIndex(frameNum);
    }

    protected int calcTickSpacing(int val) {
	return ((val/10) - (val%10));
    }

    /// Get index for the alignment with specified accuracy.
    private int getIndexFromAcc(Float a) {
	int index = -1;
	double minDiff = (double) Double.MAX_VALUE;
	double diff;
 		
	for (int i = 0; i < aligns.alignAccs.length; i++){
	    diff = Math.abs(aligns.alignAccs[i].doubleValue() - a);
	    if (diff < minDiff){
		minDiff = diff;
		index = i;
	    }
	}
 		
	return index;
    }

    /// Get index for the alignment with specified SPS.
    private int getIndexFromSPS(Float a) {
	int index = -1;
	double minDiff = Double.MAX_VALUE;
	double diff;
 		
	for (int i = 0; i < aligns.alignSPSs.length; i++){
	    diff = Math.abs(aligns.alignSPSs[i].doubleValue() - a);
	    if (diff < minDiff){
		minDiff = diff;
		index = i;
	    }
	}
 		
	return index;
    }

    /// Get index for the alignment with specified PPV.
    private int getIndexFromPPV(Float a) {
	int index = -1;
	double minDiff = Double.MAX_VALUE;
	double diff;
 		
	for (int i = 0; i < aligns.alignPPVs.length; i++){
	    diff = Math.abs(aligns.alignPPVs[i].doubleValue() - a);
	    if (diff < minDiff){
		minDiff = diff;
		index = i;
	    }
	}
 		
	return index;
    }

    /// Get index for the alignment with specified Certainty.
    private int getIndexFromCert(Float a) {
	int index = -1;
	double minDiff = Double.MAX_VALUE;
	double diff;
 		
	for (int i = 0; i < aligns.alignCerts.length; i++){
	    diff = Math.abs(aligns.alignCerts[i].doubleValue() - a);
	    if (diff < minDiff){
		minDiff = diff;
		index = i;
	    }
	}
 		
	return index;
    }

    /// Get index for the alignment with specified Consistency.
    private int getIndexFromCons(Float a) {
	int index = -1;
	double minDiff = Double.MAX_VALUE;
	double diff;
 		
	for (int i = 0; i < aligns.alignConss.length; i++){
	    diff = Math.abs(aligns.alignConss[i].doubleValue() - a);
	    if (diff < minDiff){
		minDiff = diff;
		index = i;
	    }
	}
 		
	return index;
    }

    /// Get index for the alignment with specified gap factor.
    private int getIndexFromGF(Float a) {
	int index = -1;
	double minDiff = Double.MAX_VALUE;
	double diff;
 		
	for (int i = 0; i < aligns.alignGFs.length; i++){
	    diff = Math.abs(aligns.alignGFs[i].doubleValue() - a);
	    if (diff < minDiff){
		minDiff = diff;
		index = i;
	    }
	}
 		
	return index;
    }


	//React when the user presses Enter.
	private void handleEnterKeyStroke( final JFormattedTextField field ) {
        field.getInputMap().put(KeyStroke.getKeyStroke( KeyEvent.VK_ENTER, 0), "check");
        field.getActionMap().put("check", new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
				//The text is invalid.
                if (!field.isEditValid()) { 
                    Toolkit.getDefaultToolkit().beep();
                    field.selectAll();
				//The text is valid, so use it.
                } else try {
					ok2update.put(field,true);
                    field.commitEdit();     
                } catch (java.text.ParseException exc) { }
            }
        });
	}
}

