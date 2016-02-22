package droplet_finder;
import ij.*;
import ij.plugin.*;
import ij.plugin.filter.PlugInFilter;
/*
 * Author: Samuel Moll (moll@biochem.mpg.de)
 * released under GNU GPL version 3 or above
 */
/*
 * This is a convenience plugin that first calls the Bandpass and the watershed
 * plugins and then the segment analyzer
 */


public class DF_Filterstacker implements PlugIn {
	
	public void run(String args) {
		// get current image
		ImagePlus stack = IJ.getImage();
		if(stack==null)
			return;
		
		DF_Bandpass dfbp = new DF_Bandpass();
		dfbp.run("");
		
		DF_Watershed3D dfws = new DF_Watershed3D();
		if(dfws.setup("",IJ.getImage()) != PlugInFilter.DONE) {
			dfws.run(IJ.getImage().getProcessor());
			DF_Segment_Analyzer dfsa = new DF_Segment_Analyzer();
			dfsa.run("");
		}
	}
}
