
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Michael Smoot.
 */

package mad;

import java.util.*;
import java.beans.*;


public class PropertyChangeHandler {

	private static PropertyChangeSupport pcs;
	private static Object o;

	static {
		o = new Object();
		pcs = new PropertyChangeSupport( o );
	}

	public static PropertyChangeSupport getPropertyChangeSupport() {
		return pcs;
	}

	public static void firePropertyChange(String id, Object oldValue, Object newValue) {
		try {
			PropertyChangeIDs.valueOf(id);
		} catch (IllegalArgumentException ex) {
			System.err.println("Illegal PropertyChangeEvent ID: " + id + "   ignoring!");
			return;
		}
		PropertyChangeEvent e = new PropertyChangeEvent(pcs, id, oldValue, newValue);
		pcs.firePropertyChange(e);
	}
}
