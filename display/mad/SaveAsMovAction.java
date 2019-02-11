
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 */


package mad;

import java.util.*;
import java.beans.*;
import java.io.*;
import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.*;


public class SaveAsMovAction extends AbstractAction {

	MadPanel parent;
	
	SaveAsMovAction(MadPanel parent ) {
		super("Export as Quicktime Movie");
		this.parent = parent;
	}

	public void actionPerformed(ActionEvent e) {
		this.parent.alignmentToMov();	
	}

}
