package ebmoll.lipid_droplet_counter.droplet_finder;
/*
 * Author: Samuel Moll (moll@biochem.mpg.de)
 * released under GNU GPL version 3 or above
 */
import ebmoll.lipid_droplet_counter.filters.Watershed3D;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class DF_Watershed3D implements PlugInFilter {
	ImagePlus			imp;
	static String		pluginname = "Watershed3D";
	static String		radx_str = "Radius_x";
	static String		rady_str = "Radius_y";
	static String		radz_str = "Radius_z";
	static String		invert_str = "Invert";
	float				radx;
	float				rady;
	float				radz;
	boolean				invert;
	
	public int setup(String arg, ImagePlus imp) {
		
		this.imp = imp;
		
		// show Dialog: radx, rady, radz
		GenericDialog gd = new GenericDialog(pluginname);
		radx = (float)Prefs.get(pluginname+"."+radx_str, 1.5);
		rady = (float)Prefs.get(pluginname+"."+rady_str, 1.5);
		radz = (float)Prefs.get(pluginname+"."+radz_str, 1.5);
		invert = Prefs.get(pluginname+"."+invert_str, false);
        gd.addNumericField(radx_str,radx,2);
        gd.addNumericField(rady_str,rady,2);
        gd.addNumericField(radz_str,radz,2);
        gd.addCheckbox(invert_str,invert);
        gd.showDialog();
        if (gd.wasCanceled())
            return DONE;
        radx = (float)gd.getNextNumber();
        rady = (float)gd.getNextNumber();
        radz = (float)gd.getNextNumber();
        invert = gd.getNextBoolean();
		Prefs.set(pluginname+"."+radx_str, radx);
		Prefs.set(pluginname+"."+rady_str, rady);
		Prefs.set(pluginname+"."+radz_str, radz);
		Prefs.set(pluginname+"."+invert_str, invert);
		
		return DOES_ALL;
	}
	
	public void run(ImageProcessor ip) {
		Watershed3D ws3d = new Watershed3D();
		ws3d.in_image = imp.getStack();
		ws3d.in_invert = invert;
		ws3d.in_radx = radx;
		ws3d.in_rady = rady;
		ws3d.in_radz = radz;
		ws3d.filterit();
		ImagePlus impOut= new ImagePlus("WS_" + imp.getTitle(), ws3d.out_transform);
		impOut.show();
	}
}
