
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Michael Smoot and Adam Roberts.
 */

package mad;

import java.io.*;
import java.util.*;
import javax.swing.*;
import java.awt.Dimension;

public class MAD{
	
	public static void main(String[] args){
		
		System.out.print("Reading data and building structures...");
		try {
			Alignments aligns = null;
			String fsaFile = "";
			String altFile = "";
			if (args.length == 1){
				fsaFile = args[0];
				aligns = new Alignments(fsaFile, fsaFile + ".gui", fsaFile + ".probs");
				
			} else if (args.length == 2){
				fsaFile = args[0];
				altFile = args[1];
				aligns = new Alignments(fsaFile, fsaFile + ".gui", fsaFile + ".probs", altFile);
			} else{
				System.out.println("USAGE: java -jar MAD.jar sequencefile [alternatealignfile]");
				System.exit(1);
			}
			
			System.out.println("done");

        	final JFrame frame = new JFrame("Multiple Alignment Display: " + args[0] + "  (FSA alignment)");
        	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        	MadPanel panel = new MadPanel(frame, aligns, fsaFile, altFile);
        	panel.setPreferredSize( new Dimension(900,700) );

			JMenuBar menuBar = new JMenuBar();
			JMenu fileMenu = new JMenu("File");
			menuBar.add(fileMenu);
			JMenuItem saveAsFasta = new JMenuItem( new SaveAsFastaAction(panel) );
			JMenuItem saveAsTiff = new JMenuItem( new SaveAsTiffAction(panel) );
			JMenuItem saveAsMov = new JMenuItem(new SaveAsMovAction(panel) );
			fileMenu.add(saveAsFasta);
			fileMenu.add(saveAsTiff);
			fileMenu.add(saveAsMov);
			frame.setJMenuBar(menuBar);


        	//Display the window.
        	frame.setContentPane(panel);
        	frame.pack();
        	frame.setVisible(true);
        	
        	
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}

	}
}
