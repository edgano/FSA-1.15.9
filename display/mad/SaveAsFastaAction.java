
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Michael Smoot and Adam Roberts.
 */

package mad;

import java.util.*;
import java.beans.*;
import java.io.*;
import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.Component;


public class SaveAsFastaAction extends AbstractAction implements PropertyChangeListener {

	Alignments aligns;
	MadPanel parent;
	JFileChooser fc;
	
	SaveAsFastaAction(MadPanel parent ) {
		super("Save as Multi-FASTA");
		this.parent = parent;
		fc = new JFileChooser();
		PropertyChangeHandler.getPropertyChangeSupport()
		                     .addPropertyChangeListener(PropertyChangeIDs.CHANGE_ALIGNMENT.toString(), this);
		aligns =  parent.aligns;
	}

	public void actionPerformed(ActionEvent e) {
		Alignment a = aligns.get(aligns.index);
		System.out.print ("Saving as Multi-FASTA file...");
		if ( a != null ) {
			try {
				int returnVal = fc.showSaveDialog(parent);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					FileWriter fw = new FileWriter(file);
					String mfa = a.toMultiFasta();
					fw.write(mfa,0,mfa.length());
					fw.close();
					System.out.println ("done");
				}
			} catch (Exception ioe) {
				JOptionPane.showMessageDialog(parent,
				                              "Failed to save Multi-FASTA file:\n\n" + ioe.getMessage(),
				                              "Error saving file.",
				                              JOptionPane.ERROR_MESSAGE);
				ioe.printStackTrace();
			}
		} else
			System.err.println("No alignment found.");
			
	}

	public void propertyChange(PropertyChangeEvent e) {}
}
