
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 */

package mad;

import java.io.*;
import java.util.*;
import javax.swing.*;
import java.applet.*;
import java.awt.Dimension;

public class MadApplet extends JApplet{
	public void init() {
	  	System.out.println("Applet initializing");
		String fsaFile = getParameter("fsaFile");
		Alignments aligns = new Alignments(fsaFile, fsaFile + ".gui", fsaFile + ".probs");
	  	MadPanel panel = new MadPanel(aligns, fsaFile, "");
		panel.setPreferredSize( new Dimension(900,700) );
		getContentPane().add(panel);
	}

	public void start() {
		System.out.println("Applet starting");
	}
	public void stop() {
		System.out.println("Applet stopping");
	}

	public void destroy() {
		System.out.println("Applet destroyed");
	}

	
}
