package ebmoll.lipid_droplet_counter.droplet_finder;

import ebmoll.lipid_droplet_counter.filters.Bandpass3D;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
/*
 * Author: Samuel Moll (moll@biochem.mpg.de)
 * released under GNU GPL version 3 or above
 */
public class DF_Bandpass implements PlugIn {
	static String		pluginname = "Bandpass";
	static String		hprad_str = "Maximal_feature_size";
	static String		lprad_str = "Minimal_feature_size";
	static String		zdx_str = "Z/X_aspect_ratio";
	public double		hprad;
	public double		lprad;
	public double		zdx;
	
	public void run(String args) {
		// get current image
		ImagePlus stack = IJ.getImage();
		if(stack==null)
			return;
		// show dialog and get input parameters
		GenericDialog gd = new GenericDialog(pluginname);
		gd.addMessage("Apply bandpass on image \""+stack.getTitle()+"\"");
		gd.addMessage("This plugin takes a while...");
		// get preferences for settings
		hprad = Prefs.get(pluginname+"."+hprad_str, 7.0);
		lprad = Prefs.get(pluginname+"."+lprad_str, 2.5);
		zdx = Prefs.get(pluginname+"."+zdx_str, 1.0);
		gd.addNumericField(hprad_str,hprad,2);
		gd.addNumericField(lprad_str,lprad,2);
		gd.addMessage("The next field should be 1.0 if the z-resolution (\"depth\" of one image " + "\n" +
				"of the stack) is the same as the x/y-resolution. It should be greater than one " + "\n" +
				"if your lateral resolution is better than the z-resolution. Probably you want " + "\n" +
				"this to be greater than one." );
		gd.addNumericField(zdx_str,zdx,3);
		gd.showDialog();
		if(gd.wasCanceled())
			return;
		hprad = gd.getNextNumber();
		lprad = gd.getNextNumber();
		zdx = gd.getNextNumber();
		Prefs.set(pluginname+"."+hprad_str,hprad);
		Prefs.set(pluginname+"."+lprad_str,lprad);
		Prefs.set(pluginname+"."+zdx_str,zdx);
		
		// filter image with 3D-bandpass
		Bandpass3D bp3d = new Bandpass3D();
		bp3d.in_hprad = hprad;
		bp3d.in_lprad = lprad;
		bp3d.in_xz_ratio = zdx;
		bp3d.in_image = stack.getStack();
		if(!bp3d.checkInputParams().equals("")) {
			IJ.showMessage(bp3d.checkInputParams());
			return;
		}
		bp3d.filterit();
		ImagePlus impOut= new ImagePlus("BP_" + stack.getTitle(), bp3d.out_result);
		impOut.show();
	}
}
