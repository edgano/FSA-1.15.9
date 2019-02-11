
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Robert Bradley.
 */

package mad;

import java.util.*;
import java.beans.*;
import java.io.*;
import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.Component;

public class SaveAsTiffAction extends AbstractAction implements PropertyChangeListener {

    AlignmentPanel parentAlignPanel;
    MadPanel parent;
    JFileChooser fc;
	
    SaveAsTiffAction(MadPanel parent ) {
	super("Save as TIFF file");
	this.parent = parent;
	fc = new JFileChooser();
	PropertyChangeHandler.getPropertyChangeSupport()
	    .addPropertyChangeListener(PropertyChangeIDs.CHANGE_ALIGNMENT.toString(), this);
	parentAlignPanel =  parent.alignPanel;
    }

    public void actionPerformed(ActionEvent e) {

	System.out.print("Saving as TIFF image...");

	try {
	    int returnVal = fc.showSaveDialog(parent);
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		File file = fc.getSelectedFile();
		parentAlignPanel.saveAsTiff (file.getPath());
		System.out.println("done");
	    }
	} catch (Exception ioe) {
	    JOptionPane.showMessageDialog(parent,
					  "Failed to save TIFF file:\n\n" + ioe.getMessage(),
					  "Error saving file.",
					  JOptionPane.ERROR_MESSAGE);
	    ioe.printStackTrace();
	}
	
    }

    public void propertyChange(PropertyChangeEvent e) {}
}
