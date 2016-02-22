package droplet_finder;
import ij.*;
import ij.process.*;
import ij.plugin.*;
import java.util.*;

import java.awt.geom.*;
import java.awt.*;
import ij.measure.*;

import ij.gui.*;
import ij_ImagePlusOverlay.*;
/*
 * Author: Samuel Moll (moll@biochem.mpg.de)
 * released under GNU GPL version 3 or above
 */
/*
 * This plugin analyzes 3D-Images that have been segmented. It takes as input
 * a non-color source image and a 8/16/32(24)-bit mask.
 * The mask has all pixels labeled with a unique integer  from [0 - numregions]
 * according to which region they belong to.
 * 
 * there is another filter to convert a black/white segmented image
 * to the required format.
 * 
 * The plugin performs various statistics on each region and can then do
 * processing on the original image according to the properties of the
 * regions, e.g. reject regions that contain mainly noise.
 * 
 * It can also do a gauss-fit on a region which is useful for finding
 * particles.
 * 
 * It can output the results in various formats.
 * 
 * TODO: actually implement all of this ;)
 * 
 */

/*
 * TODO: Bugs:
 *  - Leaks Memory somewhere. Probably this is ImageJ's fault ;). When previewing
 *    a lot, memory accumulates that doesn't get freed even when closing all
 *    images.
 *  - The preview scheme is very weak (because it uses the AWT thread...) and
 *    bug-prone.
 *  
 *  TODO:
 *  This needs a big rewrite....
 *  
 *  A good new structure would be:
 *  1. read image and mask (this is already done...)
 *  2. correct for oversegmentation, possibly with preview and
 *     tweakable parameters, different meaningful algorithms and
 *     so on...
 *     This is currently done in Network.threshold() and uses a very
 *     simple heuristic
 *  3. decide which segments are droplets and which are not.
 *  
 */

public class DF_Segment_Analyzer implements PlugIn, DialogListener {
	Network				network = null;
	byte[][]			particle = null;
	ImagePlusOverlay	preview_image;
	int					dimx;
	int					dimy;
	int					dimz;
	
	private float		at_last;
	private float		mt_last;
	private float		connectbias_last;
	private int			sn_last;
	
	private boolean		previewrunning = false;
	
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e){
    	// is the network already generated and the preview window available?
    	if(network==null || particle==null || preview_image==null)
    		return true;
    	// "prevent" another instance of this function to be runned
    	// this probably doesn't guarantee that only one instance of this
    	// function runs at any one time (which doesn't matter).
    	if(previewrunning)
    		return true;
    	previewrunning = true;
    	// get parameters
    	float at = (float)gd.getNextNumber();
    	float mt = (float)gd.getNextNumber();
    	float connectbias = (float)gd.getNextNumber();
    	int sn = (int)gd.getNextNumber();
    	if(gd.invalidNumber())
    		return true;
    	if(sn<1) sn = 1;
    	if(sn>dimz+1) sn = dimz+1;
    	if(at_last!=at || mt_last!=mt || connectbias_last!=connectbias) {
    		network.threshold(at, mt, connectbias);
    		at_last = at;
    		mt_last = mt;
    		connectbias_last = connectbias;
    	}
    	if(sn_last!=sn) {
    		// redraw current slice in preview window
    		float[] pix = (float[])preview_image.getProcessor().getPixels();
    		for(int i=0;i<dimx*dimy;i++) {
    			pix[i] = network.image[sn-1][i];
    		}
    		preview_image.updateAndDraw();
    	}
    	/*preview_image.getCanvas().setDisplayList(network.getParticleMarkers(sn-1),Color.red,new BasicStroke(1));
    	preview_image.repaintWindow();*/
    	network.drawParticleMarkers(preview_image.getOverlay(),sn-1);
    	previewrunning = false;
    	return true;
    }
    
	public void run(String args) {
		// show image and mask selecting dialog
		GenericDialog gd = new GenericDialog("Segment Analyzer");
 		int[] wList = WindowManager.getIDList();
        if (wList == null){
            IJ.noImage();
            return;
        }
		String[] imagelist = new String[wList.length];
        for (int i = 0; i < wList.length; i++){
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null)
                imagelist[i] = imp.getTitle();
            else
                imagelist[i] = "";
        }
		gd.addChoice("Image", imagelist, imagelist[0]);
		gd.addChoice("Mask", imagelist, imagelist[(wList.length>1)?1:0]);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        ImagePlus imp_image = WindowManager.getImage(wList[gd.getNextChoiceIndex()]);
        ImagePlus imp_mask = WindowManager.getImage(wList[gd.getNextChoiceIndex()]);
        
        // TODO: maybe lock image and mask
        
        // TODO: maybe be a bit memory-friendlier
		// convert image and mask to float[][] / int[][] stacks
		IJ.showStatus("converting image to float...");
		ImageStack stack = imp_image.getStack();
		dimx = stack.getWidth();
		dimy = stack.getHeight();
		dimz = stack.getSize();
		// this is the variable we load the image into
		float[][] image = new float[dimz][];
		if(imp_image.getProcessor() instanceof FloatProcessor) {
			for(int i = 0; i < dimz; i++){
				image[i] = (float[])stack.getProcessor(i+1).getPixels();
			}
		} else {
			for(int i = 0; i < dimz; i++){
				image[i] = (float[])stack.getProcessor(i+1).convertToFloat().getPixels();
			}
		}
		IJ.showStatus("converting mask to int...");
		stack = imp_mask.getStack();
		if(dimx!=stack.getWidth() || dimy!=stack.getHeight() || dimz!=stack.getSize()) {
			IJ.showMessage("Image and Mask have to be of the same size!");
			return;
		}
		if(!(imp_mask.getProcessor() instanceof ShortProcessor || imp_mask.getProcessor() instanceof ByteProcessor || imp_mask.getProcessor() instanceof ColorProcessor)) {
			IJ.showMessage("Mask has to be either 8-bit, 16-bit or 32-bit Color");
			return;
		}
		// this is the variable we load the mask into
		int[][] mask = new int[dimz][];
		for(int i = 0; i < dimz; i++){
			mask[i] = new int[dimx*dimy];
			Object slice = stack.getProcessor(i+1).getPixels();
			if(imp_mask.getProcessor() instanceof ShortProcessor) {
				for(int j=0;j<dimx*dimy;j++){
					mask[i][j] = ((short[])slice)[j];
					if(mask[i][j]<0)
						mask[i][j]=0x10000-mask[i][j];
				}
			}
			else if(imp_mask.getProcessor() instanceof ByteProcessor) {
				for(int j=0;j<dimx*dimy;j++){
					mask[i][j] = ((byte[])slice)[j];
					if(mask[i][j]<0)
						mask[i][j]=0x100-mask[i][j];
				}
			}
			else if(imp_mask.getProcessor() instanceof ColorProcessor) {
				for(int j=0;j<dimx*dimy;j++){
					mask[i][j] = ((int[])slice)[j] & 0x00FFFFFF;
				}
			}
		}
		// create particle bitmask
		particle = new byte[dimz][];
		for(int i = 0; i < dimz; i++){
			particle[i] = new byte[dimx*dimy];
			for(int j=0;j<dimx*dimy;j++){
				particle[i][j] = 0;
			}
		}
		
		// creating network
		IJ.showStatus("creating network...");
		network = new Network(image,mask,particle,dimx,dimy,dimz);
		
		// parameters dialog
		at_last = (float)Prefs.get("segment_analyzer.at_last", 0.2);
		mt_last = (float)Prefs.get("segment_analyzer.mt_last", 0.3);
		connectbias_last = (float)Prefs.get("segment_analyzer.connectbias_last", 0.9);
		sn_last = imp_image.getCurrentSlice();
		network.threshold(at_last, mt_last, connectbias_last);
		gd = new GenericDialog("Segment Analyzer");
		gd.addNumericField("Area_threshold", at_last, 3);
		gd.addNumericField("Maximum_threshold", mt_last,2);
		gd.addNumericField("Connect_threshold", connectbias_last,2);
		gd.addSlider("Slice", 1, dimz, imp_image.getCurrentSlice());
		gd.addDialogListener(this);
		// create preview window
		preview_image = new ImagePlusOverlay("Preview: Segment Analyzer", new FloatProcessor(dimx,dimy));
		preview_image.setDisplayRange(imp_image.getDisplayRangeMin(),imp_image.getDisplayRangeMax());
		ImageOverlay overlay = new ImageOverlay(preview_image);
		overlay.show();
		preview_image.setOverlay(overlay);
		preview_image.show();
		float[] pix = (float[])preview_image.getProcessor().getPixels();
		for(int i=0;i<dimx*dimy;i++) {
			pix[i] = network.image[sn_last-1][i];
		}
		preview_image.updateAndDraw();
		// show dialog
		gd.showDialog();
		preview_image.close();
		if(!gd.wasCanceled()) {
			//save presets
			Prefs.set("segment_analyzer.at_last", at_last);
			Prefs.set("segment_analyzer.mt_last", mt_last);
			Prefs.set("segment_analyzer.connectbias_last", connectbias_last);

			// put result in a new image
			/*IJ.showStatus("creating bitmask...");
			ImageStack stackOut = new ImageStack(dimx,dimy);
			for (int k = 0; k < dimz; k++){
				ImageProcessor oip = new ByteProcessor(dimx,dimy);
				byte[] px = (byte[])oip.getPixels();
				for (int i = 0; i < dimx*dimy; i++){
					px[i] = particle[k][i];
				}
				stackOut.addSlice(null,oip);
			}*/
			
			// measure particles...
			ResultsTable mes = network.measure();
			mes.show("Particle Results");
		}
	}
}


class Network {
	/*
	 * A network of regions as defined by the image and mask.
	 * mask has to contain all regions labeled with a unique
	 * integer. Every integer from zero to some maximal value
	 * has to be present.
	 * 
	 * TODO: add ability to release image and mask to save memory.
	 */
	private ArrayList<Region>		regions;
	private ArrayList<Particle>		particles;
	/*
	 * TODO: rewrite this, so the PCA gets generated
	 * in Calculate() or so.
	 * Caveat: it probably has a big influence on
	 * the outcome of the PCA how the different variables
	 * of the region characterization are normalized! 
	 */

	float[][]				image;
	int[][]					mask;
	byte[][]				particle;
	private int				dimx;
	private int				dimy;
	private int				dimz;
	public AutoHistogram	sumhist = null;
	public AutoHistogram	imagehist = null;
	
	public Network(float[][] image, int[][] mask, byte[][] particle, int dimx, int dimy, int dimz) {
		/*
		 * generates a Region-Network from source image and mask.
		 */
		regions = new ArrayList<Region>();
		this.image = image;
		this.mask = mask;
		this.particle = particle;
		this.dimx = dimx;
		this.dimy = dimy;
		this.dimz = dimz;
		Calculate();
	}

	public void Calculate() {
		/*
		 * calculates the Network from scratch
		 */
		regions.clear();
		/*
		 * loop through the image and initialize the Region array,
		 * find neighbors.
		 */
		imagehist = new AutoHistogram();
		for( int z = 0; z < dimz; z++) {
			IJ.showProgress(z,dimz);
			for( int y = 0; y < dimy; y++) {
				for( int x = 0; x < dimx; x++) {
					int pmask = mask[z][y*dimx+x];
					// create Regions if neccesary
					while(pmask >= regions.size())
						regions.add(new Region(regions.size(),image,mask,particle,dimx,dimy,dimz));
					regions.get(pmask).add(x,y,z);
					// test for connectivity
					ArrayList<int[]> n = NBH.getNeighborhood3Dx6(x,y,z,dimx,dimy,dimz);
					for(int i=0;i<n.size();i++) {
						int[] pos = n.get(i);
						int val = mask[pos[2]][pos[1]*dimx+pos[0]];
						while(val >= regions.size())
							regions.add(new Region(regions.size(),image,mask,particle,dimx,dimy,dimz));
						if(val!=pmask)
							regions.get(pmask).addneighbor(regions.get(val));
					}
					// calculate image histogram
					imagehist.add(image[z][y*dimx+x]);
				}
			}
		}
		// calculate histogram of region sums
		sumhist = new AutoHistogram();
		for(int i=0;i<regions.size();i++) {
			Region reg = regions.get(i);
			sumhist.add(reg.hist.getSum());
		}
		// HACK: recalculate region histograms...
		/*for(int i=0;i<regions.size();i++)
			regions.get(i).hist = new AutoHistogram(imagehist.getMin(),imagehist.getMax());
		imagehist = new AutoHistogram(imagehist.getMin(),imagehist.getMax());
		for( int z = 0; z < dimz; z++) {
			IJ.showProgress(z,dimz);
			for( int y = 0; y < dimy; y++) {
				for( int x = 0; x < dimx; x++) {
					int pmask = mask[z][y*dimx+x];
					regions.get(pmask).hist.add(image[z][y*dimx+x]);
					imagehist.add(image[z][y*dimx+x]);
				}
			}
		}*/
		
		/*
		 * TODO:
		 *  Why this doesn't work:
		 *   - Maybe the Histograms have to be normalized
		 *     by area or so
		 *   - Maybe the Histograms contain not enough data
		 *   - Maybe it works, only it looks strange
		 *   - Maybe some sort of color-correction has to
		 *     be applied, so all droplets use the same colors
		 *   - Maybe the data has to be passed in a way that
		 *     allows rescaling and translation by the PCA
		 * 
		 * 
		 * 
		 */
		
		/*ImagePlus img = imagehist.display();
		//regions.get(866).hist.displayInHistogram(img);
		regions.get(287).hist.displayInHistogram(img);
		regions.get(6198).hist.displayInHistogram(img);
		IJ.log("Min:"+imagehist.getMin());
		IJ.log("Max:"+imagehist.getMax());
		
		// HACK: create data from Histograms...
		double data[][] = new double[regions.size()][256];
		float stepsize = (imagehist.getMax()-imagehist.getMin())/256;
		for(int i=0;i<regions.size();i++) {
			for(int m=0;m<256;m++) {
				float t = imagehist.getMin()+m*stepsize;
				data[i][m] = regions.get(i).hist.getBin(t,t+stepsize);
			}
		}
		PCA pca = new PCA(data);
		double[][] EigenVec = pca.getEigenVec();
		double[] EigenVal = pca.getEigenVal();
		//for(int m=0;m<10;m++) {
			AutoHistogram HISTOBLA = new AutoHistogram(0.0f,256.0f);
			for(int i=0;i<256;i++) {
				int qnumb = (int)(EigenVec[0][i]*300);
				for(int q=0;q<qnumb;q++) {
					HISTOBLA.add(i);
				}
			}
			HISTOBLA.displayInHistogram(img);
			//IJ.log("Eigval:"+EigenVal[m]);
		//}
		 */
		/*LMAffineGaussian func = new LMAffineGaussian();
			double[] params = {276.0, 243.0, 19.0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 1000.0};
			BoundingBox bbox = new BoundingBox();
			bbox.expand(0, 0, 0);
			bbox.expand(dimx-1, dimy-1, dimz-1);
			LMfit(image,dimx,dimy,dimz,bbox,func,params);
			for( int z = bbox.minz; z <= bbox.maxz; z++) {
				for( int y = bbox.miny; y <= bbox.maxy; y++) {
					for( int x = bbox.minx; x <= bbox.maxx; x++) {
						double[] pos = {x,y,z};
						image[z][y*dimx+x] -= func.val(pos, params);
					}
				}
			}*/
	}
	
	private void CalculateParticleList() {
		// flood-fill all regions and generate particle lists 
		particles = new ArrayList<Particle>();
		
		boolean[] processed = new boolean[regions.size()];
		LinkedList<Region> fifo = new LinkedList<Region>();
		for(int i=0;i<regions.size();i++) {
			if(processed[i])
				continue;
			//start new particle...
			Particle newparticle = new Particle();
			particles.add(newparticle);
			fifo.clear();
			fifo.addLast(regions.get(i));
			processed[i]=true;
			// process fifo until empty...
			while(!fifo.isEmpty()) {
				Region reg = fifo.removeFirst();
				newparticle.addRegion(reg);
				// add all unprocessed connected regions to fifo
				for(int k=0;k<reg.connected.size();k++) {
					Region con = reg.connected.get(k);
					if(!processed[con.index]) {
						fifo.addLast(con);
						processed[con.index]=true;
					}
				}
			}
		}
	}

	/**
	 * decides which voxels belong to a particle.
	 * minArea is the minimal luminosity a region must have
	 * minMax is the minimal value of its maximum 
	 */
	public void threshold(float minArea, float minMax, float fwhmbias) {
		// loop through the regions and calculate stuff...
		
		// first thresholds each individual region with a
		// "biased" FWHM (not half maximum, but 0.0-1.0
		// times the maximum) then connects all regions
		// that touch each other.
		
		// then it selects those regions that fulfill the two
		// min...-thresholds and masks them with FWHM
		
		// biased threshold of all regions
		// & clear connections between regions
		for(int i=0;i<regions.size();i++) {
			Region reg = regions.get(i);
			reg.threshold(reg.hist.getMax()*fwhmbias, false);
			reg.disconnectAll();
		}
		
		// loop over all voxels and connect regions...
		for( int z = 0; z < dimz; z++) {
			for( int y = 0; y < dimy; y++) {
				for( int x = 0; x < dimx; x++) {
					if(particle[z][y*dimx+x]==0)
						continue;
					int pmask = mask[z][y*dimx+x];
					// test for connectivity
					ArrayList<int[]> n = NBH.getNeighborhood3Dx6(x,y,z,dimx,dimy,dimz);
					for(int i=0;i<n.size();i++) {
						int[] pos = n.get(i);
						if(particle[pos[2]][pos[1]*dimx+pos[0]]==0)
							continue;
						int val = mask[pos[2]][pos[1]*dimx+pos[0]];
						if(val!=pmask)
							regions.get(pmask).connect(regions.get(val));
					}
				}
			}
		}
		
		// calculate particle list by floodfilling
		CalculateParticleList();
		
		//float thres = sumhist.getMean()+4.0f*sumhist.getDeviation();
		float thres = sumhist.getMin()+minArea*(sumhist.getMax()-sumhist.getMin());
		float minmaxthres = imagehist.getMin()+minMax*(imagehist.getMax()-imagehist.getMin());
		for(int i=0;i<regions.size();i++) {
			Region reg = regions.get(i);
			// reject regions based on thresholds
			// do a FWHM on the regions that pass the test
			if(reg.hist.getSum()>=thres && reg.hist.getMax()>=minmaxthres)
				reg.threshold(reg.hist.getMax()/2.0f, false);
			// set remaining regions to zero...
			else
				//reg.threshold(reg.hist.getMax()+1.0f, false);
				reg.setparticlenull();
		}
	}

	public void threshold_PCA() {
		/*
		 * TODO: write a thresholding function based on
		 * the PCA-data
		 */

	}
	
	/*
	 * return a GeneralPath that marks all thresholded particles with a red circle
	 * that approximates their size in the slice z.
	 */
	public GeneralPath getParticleMarkers(int z) {
		// TODO: rewrite so that the x/z-ratio is correctly considered
		GeneralPath path = new GeneralPath();
		for(int i=0;i<particles.size();i++) {
			Particle reg = particles.get(i);
			if(reg.bbox.minz<=z && reg.bbox.maxz>=z) {
				float[] pos = reg.position();
				int volume = reg.volume();
				if( volume > 0 ) {
					float rad;
					if(dimz == 1)
						rad = (float)java.lang.Math.pow(volume/3.1415,1.0/2.0);
					else
						rad = (float)java.lang.Math.pow(volume/4.1,1.0/3.0);
					rad = (rad*rad-(z-pos[2])*(z-pos[2]))+1;
					if(rad>0) {
						rad=(float)java.lang.Math.sqrt(rad);
						path.append(new Ellipse2D.Float(pos[0]-rad, pos[1]-rad, rad*2f, rad*2f), false);
					}
				}
			}
		}
		return path;
	}
	
	public void drawParticleMarkers(ImageOverlay oi, int slice) {
		// create LUT for particle <--> region correlation
		int[] particleindex = new int[regions.size()];
		for(int i=0;i<particles.size();i++) {
			Particle part = particles.get(i);
			for(int k=0;k<part.regions.size();k++) {
				Region reg = part.regions.get(k);
				particleindex[reg.index] = i;
			}
		}
		for( int y = 0; y < dimy; y++) {
			for( int x = 0; x < dimx; x++) {
				boolean border = false;
				int partnum = particleindex[mask[slice][y*dimx+x]];
				// look at neighbors
				ArrayList<int[]> n = NBH.getNeighborhood2Dx4(x,y,dimx,dimy);
				for(int i=0;i<n.size();i++) {
					int[] pos = n.get(i);
					int npartnum = particleindex[mask[slice][pos[1]*dimx+pos[0]]];
					if(npartnum!=partnum)
						border = true;
				}
				if(border) {
					oi.getPixels()[y*dimx+x]=0x00000000;//0x4400FF00;
				} else {
					if(particle[slice][y*dimx+x]!=0)
						oi.getPixels()[y*dimx+x]=0xFFFF6000;
					else
						oi.getPixels()[y*dimx+x]=0x00000000;
				}
			}
		}
		// redraw overlay
		oi.reset();
	}
	
	public ResultsTable measure() {
		/*
		 * measure all kinds of things about the particles and output them in a
		 * ij.measure.ResultsTable
		 * things to be measured:
		 *  - volume
		 *  - position
		 *  - surface area
		 *  - shape
		 */
		ResultsTable mes = new ResultsTable();
		if(particles.size()==0)
			return mes;
		int i_volume = mes.getFreeColumn("volume");
		int i_posx = mes.getFreeColumn("position_x");
		int i_posy = mes.getFreeColumn("position_y");
		int i_posz = mes.getFreeColumn("position_z");
		int i_surface = mes.getFreeColumn("surface area");
		//int i_IQ = mes.getFreeColumn("isoperimetric quotient");
		IJ.showStatus("measuring particles...");
		for(int i=0;i<particles.size();i++) {
			IJ.showProgress(i,particles.size());
			Particle part = particles.get(i);
			float vol = part.volume();
			if(vol>0) {
				mes.incrementCounter();
				mes.addValue(i_volume,vol);
				float[] pos = part.position();
				mes.addValue(i_posx,pos[0]);
				mes.addValue(i_posy,pos[1]);
				mes.addValue(i_posz,pos[2]+1.0f); // stacks begin with 1
				float area = part.surfacearea();
				mes.addValue(i_surface,area);
				//mes.addValue(i_IQ, 36*3.1415926536*vol*vol/area/area/area);
			}
		}
		return mes;
	}
}

class NBH {
	static ArrayList<int[]> getNeighborhood3Dx6(int x, int y, int z, int dimx, int dimy, int dimz) {
		ArrayList<int[]> result = new ArrayList<int[]>();
		int i; int j; int k;
		i = x-1; j=y; k=z;
		if(i>=0 && i<dimx && j>=0 && j<dimy && k>=0 && k<dimz)
			result.add(new int[] {i,j,k});
		i = x+1; j=y; k=z;
		if(i>=0 && i<dimx && j>=0 && j<dimy && k>=0 && k<dimz)
			result.add(new int[] {i,j,k});
		i = x; j=y-1; k=z;
		if(i>=0 && i<dimx && j>=0 && j<dimy && k>=0 && k<dimz)
			result.add(new int[] {i,j,k});
		i = x; j=y+1; k=z;
		if(i>=0 && i<dimx && j>=0 && j<dimy && k>=0 && k<dimz)
			result.add(new int[] {i,j,k});
		i = x; j=y; k=z-1;
		if(i>=0 && i<dimx && j>=0 && j<dimy && k>=0 && k<dimz)
			result.add(new int[] {i,j,k});
		i = x; j=y; k=z+1;
		if(i>=0 && i<dimx && j>=0 && j<dimy && k>=0 && k<dimz)
			result.add(new int[] {i,j,k});
		return result;
	}
	
	static ArrayList<int[]> getNeighborhood2Dx8(int x, int y, int dimx, int dimy) {
		ArrayList<int[]> result = new ArrayList<int[]>();
		int i; int j;
		i = x-1; j=y;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x-1; j=y-1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x-1; j=y+1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x+1; j=y;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x+1; j=y-1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x+1; j=y+1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x; j=y-1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x; j=y+1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		return result;
	}
	
	static ArrayList<int[]> getNeighborhood2Dx4(int x, int y, int dimx, int dimy) {
		ArrayList<int[]> result = new ArrayList<int[]>();
		int i; int j;
		i = x-1; j=y;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x+1; j=y;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x; j=y-1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		i = x; j=y+1;
		if(i>=0 && i<dimx && j>=0 && j<dimy)
			result.add(new int[] {i,j});
		return result;
	}
}


/*
 * One particle consists of several regions.
 * This is to correct for Oversegmentation.
 */
class Particle {
		ArrayList<Region>		regions;
		BoundingBox				bbox;
		
		public Particle() {
			regions = new ArrayList<Region>();
			bbox = new BoundingBox();
		}
		
		public void addRegion(Region region) {
			regions.add(region);
			bbox.expand(region.bbox);
		}
		
		public int volume() {
			// returns number of voxels the particle contains
			int volume=0;
			for(int i=0;i<regions.size();i++) {
				volume +=regions.get(i).volume();
			}
			return volume;
		}
		
		public float[] position() {
			// returns the mean position of the particle
			float[] pos = {0.0f,0.0f,0.0f};
			int sum = 0;
			for(int i=0;i<regions.size();i++) {
				Region reg = regions.get(i);
				for( int z = reg.bbox.minz; z <= reg.bbox.maxz; z++)
					for( int y = reg.bbox.miny; y <= reg.bbox.maxy; y++)
						for( int x = reg.bbox.minx; x <= reg.bbox.maxx; x++)
							if(reg.mask[z][y*reg.dimx+x] == reg.index && reg.particle[z][y*reg.dimx+x] == -1) {
								pos[0]+=x;
								pos[1]+=y;
								pos[2]+=z;
								sum++;
							}
			}
			if(sum==0) {
				pos[0] = -1.0f;
				pos[1] = -1.0f;
				pos[2] = -1.0f;
				return pos;
			}
			pos[0] /= sum;
			pos[1] /= sum;
			pos[2] /= sum;
			return pos;
		}
		
		public float surfacearea() {
			/*
			 * estimate surface area
			 * classify each voxel according to its surface sides and assign
			 * each class a surface weight taken from the literature:
			 * see "Voxel-based surface area estimation: from theory to practice"
			 * G. Windreich et.al.
			 * 
			 * TODO:
			 * currently particle surface that is on the surface of
			 * the voxel volume (e.g. of particles that are only
			 * partially in the picture) is not counted. Maybe make
			 * an option to include this part of the surface.
			 * 
			 * TODO: Limitation:
			 * probably the surface area is quite biased for
			 * 2D-images (i.e. 1-image-stacks)
			 */
			float area = 0.0f;
			for( int z = bbox.minz; z <= bbox.maxz; z++)
				for( int y = bbox.miny; y <= bbox.maxy; y++)
					for( int x = bbox.minx; x <= bbox.maxx; x++)
						if(isInParticle(x,y,z)) {
							// classify voxel 
							int countx = 0;
							if(!isInParticle(x-1,y,z))
								countx++;
							if(!isInParticle(x+1,y,z))
								countx++;
							int county = 0;
							if(!isInParticle(x,y-1,z))
								county++;
							if(!isInParticle(x,y+1,z))
								county++;
							int countz = 0;
							if(!isInParticle(x,y,z-1))
								countz++;
							if(!isInParticle(x,y,z+1))
								countz++;
							// add area according to classification:
							//class 1: W1=0.894f
							if(countx+county+countz==1) {
								area+=0.894f;
							}
							//class 2: W2=1.3409f
							else if(countx+county+countz==2) {
								if(countx!=2&&county!=2&&countz!=2)
									area+=1.3409f;
							}
							//class 3: W3=1.5879f
							else if(countx+county+countz==3) {
								if(countx!=2&&county!=2&&countz!=2)
									area+=1.5879f;
							}
							//class 4: W4=2f
							else if(countx+county+countz==3) {
								if(countx==2||county==2||countz==2)
									area+=2.0f;
							}
							//class 5: W5=8.f/3.f
							else if(countx+county+countz==4) {
								if(countx!=0&&county!=0&&countz!=0)
									area+=2.6667f;
							}
							//class 6: W6=10.f/3.f
							else if(countx+county+countz==5) {
								area+=3.3333f;
							}
							//class 7: W7=1.79f
							else if(countx+county+countz==2) {
								if(countx==2||county==2||countz==2)
									area+=1.79f;
							}
							//class 8: W8=2.68f
							else if(countx+county+countz==4) {
								if(countx==0||county==0||countz==0)
									area+=2.68f;
							}
							//class 9: W9=4.08f
							else if(countx+county+countz==6) {
								area+=4.08f;
							}
						}
			return area;
		}
		
		/*
		 * helper function for surface estimation
		 * returns true if (x,y,z) is part of the particle
		 */
		boolean isInParticle(int x, int y, int z) {
			if(regions.size()==0)
				return false;
			if(x<0 || x>=regions.get(0).dimx)
				return false;
			if(y<0 || y>=regions.get(0).dimy)
				return false;
			if(z<0 || z>=regions.get(0).dimz)
				return false;
			for(int i=0;i<regions.size();i++) {
				Region reg = regions.get(i);
				int[][] mask = reg.mask;
				byte[][] particle = reg.particle;
				if(mask[z][y*reg.dimx+x] == reg.index && particle[z][y*reg.dimx+x] == -1)
					return true;
			}
			return false;
		}
}

class Region {
	/*
	 * A Region in the Image as represented in the mask.
	 * Contains various facts about the Region
	 */
	float[][]		image; // references to image and mask
	int[][]			mask; // so Region can do e.g. a Gauss-fit
	byte[][]		particle; // so region can set particle yes/no
	int				dimx;
	int				dimy;
	int				dimz;
	int				index;	// the corresponding color in the mask
	
	BoundingBox		bbox;	// Bounding Box for efficiency purposes
	AutoHistogram	hist;
	ArrayList<Region>	neighbors;
	ArrayList<Region>	connected;

	public Region(int index, float[][] image, int[][] mask, byte[][] particle, int dimx, int dimy, int dimz) {
		this.image = image;
		this.mask = mask;
		this.particle = particle;
		this.dimx = dimx;
		this.dimy = dimy;
		this.dimz = dimz;
		this.index = index;
		reset();
	}
	public void reset() {
		bbox = new BoundingBox();
		hist = new AutoHistogram();
		neighbors = new ArrayList<Region>();
		connected = new ArrayList<Region>();
	}
	public void add(int x, int y, int z) {
		bbox.expand(x, y, z);
		hist.add(image[z][y*dimx+x]);
	}
	public void disconnectAll() {
		connected = new ArrayList<Region>();
	}
	public void connect(Region region) {
		if(!connected.contains(region))
			connected.add(region);
		if(!region.connected.contains(this))
			region.connected.add(this);
	}
	public void addneighbor(Region region) {
		if(!neighbors.contains(region))
			neighbors.add(region);
		if(!region.neighbors.contains(this))
			region.neighbors.add(this);
	}
	public void eat(Region region) {
		/*
		 * sets all mask-pixels of region to this.index, updates
		 * connectivities. Sets region to empty.
		 */
		for( int z = region.bbox.minz; z <= region.bbox.maxz; z++) {
			for( int y = region.bbox.miny; y <= region.bbox.maxy; y++) {
				for( int x = region.bbox.minx; x <= region.bbox.maxx; x++) {
					if(mask[z][y*dimx+x] == region.index) {
						add(x,y,z);
						mask[z][y*dimx+x] = index;
					}
				}
			}
		}
		for(int i=0;i<region.neighbors.size();i++) {
			addneighbor(region.neighbors.get(i));
		}
		region.reset();
	}

	public void setNull() {
		for( int z = bbox.minz; z <= bbox.maxz; z++)
			for( int y = bbox.miny; y <= bbox.maxy; y++)
				for( int x = bbox.minx; x <= bbox.maxx; x++) {
					if(mask[z][y*dimx+x] == index)
						image[z][y*dimx+x] = 0.0f;
				}
	}
	public void setparticlenull() {
		for( int z = bbox.minz; z <= bbox.maxz; z++)
			for( int y = bbox.miny; y <= bbox.maxy; y++)
				for( int x = bbox.minx; x <= bbox.maxx; x++)
					if(mask[z][y*dimx+x] == index) {
						particle[z][y*dimx+x] = 0;
					}
	}
	public void threshold(float thres, boolean invert) {
		// if invert is set, all pixels below the threshold are set to particle
		float sign = invert?-1.0f:1.0f;
		for( int z = bbox.minz; z <= bbox.maxz; z++)
			for( int y = bbox.miny; y <= bbox.maxy; y++)
				for( int x = bbox.minx; x <= bbox.maxx; x++) {
					if(mask[z][y*dimx+x] == index) {
						if(sign*image[z][y*dimx+x]>sign*thres)
							particle[z][y*dimx+x] = -1;
						else
							particle[z][y*dimx+x] = 0;
					}
				}
	}
	public int volume() {
		// returns number of voxels the particle contains
		int sum = 0;
		for( int z = bbox.minz; z <= bbox.maxz; z++)
			for( int y = bbox.miny; y <= bbox.maxy; y++)
				for( int x = bbox.minx; x <= bbox.maxx; x++)
					if(mask[z][y*dimx+x] == index && particle[z][y*dimx+x] == -1)
						sum++;
		return sum;
	}
	public float[] position() {
		// returns the mean position of the particle
		float[] pos = {0.0f,0.0f,0.0f};
		int sum = 0;
		for( int z = bbox.minz; z <= bbox.maxz; z++)
			for( int y = bbox.miny; y <= bbox.maxy; y++)
				for( int x = bbox.minx; x <= bbox.maxx; x++)
					if(mask[z][y*dimx+x] == index && particle[z][y*dimx+x] == -1) {
						pos[0]+=x;
						pos[1]+=y;
						pos[2]+=z;
						sum++;
					}
		if(sum==0) {
			pos[0] = -1.0f;
			pos[1] = -1.0f;
			pos[2] = -1.0f;
			return pos;
		}
		pos[0] /= sum;
		pos[1] /= sum;
		pos[2] /= sum;
		return pos;
	}
	public float surfacearea() {
		/*
		 * estimate surface area
		 * classify each voxel according to its surface sides and assign
		 * each class a surface weight taken from the literature:
		 * see "Voxel-based surface area estimation: from theory to practice"
		 * G. Windreich et.al.
		 * 
		 * TODO:
		 * currently particle surface that is on the surface of
		 * the voxel volume (e.g. of particles that are only
		 * partially in the picture) is not counted. Maybe make
		 * an option to include this part of the surface.
		 * 
		 * TODO: Limitation:
		 * probably the surface area is quite biased for
		 * 2D-images (i.e. 1-image-stacks)
		 */


		float area = 0.0f;
		for( int z = bbox.minz; z <= bbox.maxz; z++)
			for( int y = bbox.miny; y <= bbox.maxy; y++)
				for( int x = bbox.minx; x <= bbox.maxx; x++)
					if(mask[z][y*dimx+x] == index && particle[z][y*dimx+x] == -1) {
						// classify voxel 
						int countx = 0;
						if(x-1>=0) {
							if(mask[z][y*dimx+x-1]!=index || particle[z][y*dimx+x-1]!=-1) {
								countx++; }
						} else countx++;
						if(x+1<dimx) {
							if(mask[z][y*dimx+x+1]!=index || particle[z][y*dimx+x+1]!=-1) {
								countx++; }
						} else countx++;
						int county = 0;
						if(y-1>=0) {
							if(mask[z][(y-1)*dimx+x]!=index || particle[z][(y-1)*dimx+x]!=-1) {
								county++; }
						} else county++;
						if(y+1<dimy) {
							if(mask[z][(y+1)*dimx+x]!=index || particle[z][(y+1)*dimx+x]!=-1) {
								county++; }
						} else county++;
						int countz = 0;
						if(z-1>=0) {
							if(mask[z-1][y*dimx+x]!=index || particle[z-1][y*dimx+x]!=-1) {
								countz++; }
						} else countz++;
						if(z+1<dimz) {
							if(mask[z+1][y*dimx+x]!=index || particle[z+1][y*dimx+x]!=-1) {
								countz++; }
						} else countz++;
						// add area according to classification:
						//class 1: W1=0.894f
						if(countx+county+countz==1) {
							area+=0.894f;
						}
						//class 2: W2=1.3409f
						else if(countx+county+countz==2) {
							if(countx!=2&&county!=2&&countz!=2)
								area+=1.3409f;
						}
						//class 3: W3=1.5879f
						else if(countx+county+countz==3) {
							if(countx!=2&&county!=2&&countz!=2)
								area+=1.5879f;
						}
						//class 4: W4=2f
						else if(countx+county+countz==3) {
							if(countx==2||county==2||countz==2)
								area+=2.0f;
						}
						//class 5: W5=8.f/3.f
						else if(countx+county+countz==4) {
							if(countx!=0&&county!=0&&countz!=0)
								area+=2.6667f;
						}
						//class 6: W6=10.f/3.f
						else if(countx+county+countz==5) {
							area+=3.3333f;
						}
						//class 7: W7=1.79f
						else if(countx+county+countz==2) {
							if(countx==2||county==2||countz==2)
								area+=1.79f;
						}
						//class 8: W8=2.68f
						else if(countx+county+countz==4) {
							if(countx==0||county==0||countz==0)
								area+=2.68f;
						}
						//class 9: W9=4.08f
						else if(countx+county+countz==6) {
							area+=4.08f;
						}
					}
		return area;
	}
}

class BoundingBox {
	public int		minx;
	public int		miny;
	public int		minz;
	public int		maxx; // max is part of the BB
	public int		maxy;
	public int		maxz;
	private boolean initialized;
	public BoundingBox() {
		initialized = false;
	}
	public void expand(int x, int y, int z) {
		// expand BB to include (x,y,z)
		if(!initialized) {
			minx = x;
			miny = y;
			minz = z;
			maxx = x;
			maxy = y;
			maxz = z;
			initialized = true;
			return;
		}
		if(x < minx)
			minx = x;
		else if(x > maxx)
			maxx = x;
		if(y < miny)
			miny = y;
		else if(y > maxy)
			maxy = y;
		if(z < minz)
			minz = z;
		else if(z > maxz)
			maxz = z;
	}
	public void expand(BoundingBox bbox) {
		if(!initialized) {
			minx = bbox.minx;
			miny = bbox.miny;
			minz = bbox.minz;
			maxx = bbox.maxx;
			maxy = bbox.maxy;
			maxz = bbox.maxz;
			initialized = true;
			return;
		}
		if(bbox.minx < minx)
			minx = bbox.minx;
		if(bbox.maxx > maxx)
			maxx = bbox.maxx;
		if(bbox.miny < miny)
			miny = bbox.miny;
		if(bbox.maxy > maxy)
			maxy = bbox.maxy;
		if(bbox.minz < minz)
			minz = bbox.minz;
		if(bbox.maxz > maxz)
			maxz = bbox.maxz;
	}
	public int volume() {
		return (maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1);
	}
}

class AutoHistogram {
	/*
	 * Automatic Histogram Class. No need to set min and max,
	 * just set the number of bins (default=256, min=2) and
	 * add values. The histogram will expand as neccessary.
	 */
	private int[] bins;
	private int numberOfBins;
	private int init;
	private float lower; // the lower border of bin 0
	private float increment; // the bin size
	private boolean fixed; // are lower and increment set by the user?
	
	private void init(int numberOfBins) {
		numberOfBins = (numberOfBins<2)?2:numberOfBins;
		bins = new int[numberOfBins];
		for(int i=0;i<numberOfBins;i++) {
			bins[i] = 0;
		}
		this.numberOfBins = numberOfBins;
		lower = 0.0f;
		increment = 0.0f;
		init = 0;
	}
	private void setfixed(float lower, float upper) {
		this.fixed = true;
		this.lower = lower;
		this.increment = (upper-lower)/numberOfBins;
		this.init = 2;
	}
	public AutoHistogram() {
		init(256);
	}
	public AutoHistogram(int numberOfBins) {
		init(numberOfBins);
	}
	public AutoHistogram(int numberOfBins,float lower,float upper) {
		init(numberOfBins);
		setfixed(lower,upper);
	}
	public AutoHistogram(float lower,float upper) {
		init(256);
		setfixed(lower,upper);
	}
	public void add(float value) {
		/*
		 * add value to the histogram
		 */
		if(init == 0) {
			bins[0] = 1;
			lower = value;
			increment = 0.0f;
			init = 1;
		} else if(init==1) {
			if(value>lower){
				bins[numberOfBins-1]=1;
				increment = (value-lower)/(numberOfBins-1);
				init = 2;
			} else if(value<lower) {
				increment = (lower-value)/(numberOfBins-1);
				bins[numberOfBins-1]=bins[0];
				bins[0]=1;
				lower = value;
				init = 2;
			} else {
				bins[0]+=1;
			}
		} else {
			int n = (int)java.lang.Math.floor((value-lower)/increment);
			if(n<0) {
				if(fixed) {
					n = 0;
				} else {
					int count = (-n-1)/numberOfBins+2;
					rebinlower(count);
					n = (n+(count-1)*numberOfBins)/count;
				}
			} else if(n>=numberOfBins) {
				if(fixed) {
					n = numberOfBins-1;
				} else {
					int count = n/numberOfBins+1;
					rebinupper(count);
					n = n/count;
				}
			}
			bins[n] += 1;
		}
	}
	public int getBin(float value) {
		/*
		 * returns the number of entries in the bin that belongs to value
		 */
		if(init==1)
			return bins[0];
		int n = (int)java.lang.Math.floor((value-lower)/increment);
		if(n<0 || n >=numberOfBins)
			return 0;
		return bins[n];
	}
	public float getBin(float Min,float Max) {
		/*
		 * returns the number of entries between Min and Max
		 * half included bins are linearly interpolated
		 */
		if(Max<Min)
			return 0.0f;
		if(init==0)
			return 0.0f;
		if(init==1) {
			if(Min<=lower && Max>lower)
				return bins[0];
			else
				return 0.0f;
		}
		if(Min<lower) {
			Min=lower;
		}
		if(Max>lower+increment*numberOfBins) {
			Max = lower+increment*numberOfBins;
		}
		
		int Min_n = (int)java.lang.Math.floor((Min-lower)/increment);
		int Max_n = (int)java.lang.Math.floor((Max-lower)/increment);
		if(Min_n<0) Min_n=0;
		if(Min_n>=numberOfBins) return 0.0f;
		if(Max_n>=numberOfBins) Max_n=numberOfBins-1;
		if(Max_n<0) return 0.0f;
		float sum = 0.0f;
		for(int n=Min_n;n<=Max_n;n++) {
			sum += bins[n];
		}
		sum-=((Min-lower)/increment - Min_n) * bins[Min_n];
		sum-=(Max_n - (Max-lower)/increment + 1) * bins[Max_n];
		return sum;
	}
	public float getMin() {
		/*
		 * get (approx.) minimum value that was ever put into
		 * the AutoHistogram.
		 */
		if(init == 0)
			return Float.POSITIVE_INFINITY;
		int i = 0;
		while(bins[i]==0)
			i+=1;
		return lower+i*increment;
	}
	public float getMax() {
		/*
		 * get the (approx.) maximum value that was ever put into
		 * the AutoHistogram
		 */
		if(init == 0)
			return Float.NEGATIVE_INFINITY;
		int i = numberOfBins-1;
		while(bins[i]==0)
			i-=1;
		return lower+(i+1)*increment;
	}
	public float getMean() {
		/*
		 * get the (approx.) mean of the distribution.
		 */
		float sum = 0.0f;
		float val = lower+0.5f*increment;
		int number = 0;
		for(int i=0;i<numberOfBins;i++) {
			sum+=bins[i]*val;
			number+=bins[i];
			val += increment;
		}
		return sum/number;
	}
	public float getDeviation() {
		/*
		 * get the (approx.) standard deviation of the distribution.
		 */
		float mean = getMean();
		float val = lower+0.5f*increment;
		float quadsum = 0.0f;
		int number = 0;
		for(int i=0;i<numberOfBins;i++) {
			quadsum+=bins[i]*(val-mean)*(val-mean);
			number+=bins[i];
			val += increment;
		}
		return (float)java.lang.Math.sqrt(quadsum/number);
	}
	public float getMedian() {
		/*
		 * returns the median of the distribution
		 */
		int number = getNumber()/2;
		int count = 0;
		float median = 0.0f;
		for(int i=0;i<numberOfBins;i++) {
			count+=bins[i];
			if(count>=number) {
				median = lower+i*increment;
				break;
			}
		}
		return median;
	}
	public int getNumber() {
		/*
		 * gets the number of values that have been put into
		 * the AutoHistogram
		 */
		int number = 0;
		for(int i=0;i<numberOfBins;i++) {
			number+=bins[i];
		}
		return number;
	}
	public int getMaxCount() {
		/*
		 * get the highest number of entries any bin has.
		 */
		int max = 0;
		for (int i = 0; i < numberOfBins; i++){
			if(bins[i]>max)
				max = bins[i];
		}
		return max;
	}
	public float getSum() {
		/*
		 * returns the (approx.) sum of all values in the distribution
		 */
		float sum = 0.0f;
		for (int i = 0; i < numberOfBins; i++){
			sum += bins[i]*(lower+((float)i+0.5f)*increment);
		}
		return sum;
	}
	private void rebinupper(int count) {
		int sum = 0;
		int p = 0;
		int c = 0;
		for(int i=0;i<numberOfBins;i++) {
			sum += bins[i];
			c++;
			if(c == count) {
				bins[p] = sum;
				sum = 0;
				p++;
				c = 0;
			}
		}
		while(p<numberOfBins) {
			bins[p] = sum;
			sum = 0;
			p++;
		}
		increment*=count;
	}
	private void rebinlower(int count) {
		int sum = 0;
		int p = numberOfBins-1;
		int c = 0;
		for(int i=numberOfBins-1;i>=0;i--) {
			sum += bins[i];
			c++;
			if(c == count) {
				bins[p] = sum;
				sum = 0;
				p--;
				c = 0;
			}
		}
		while(p>=0) {
			bins[p] = sum;
			sum = 0;
			p--;
		}
		lower-=increment*(count-1)*numberOfBins;
		increment*=count;
	}
	public ImagePlus display() {
		// put result in a new image
		ImageProcessor ip = new ByteProcessor(numberOfBins,numberOfBins);
		byte[] px = (byte[])ip.getPixels();
		for (int i = 0; i < numberOfBins*numberOfBins; i++){
			px[i] = 0;
		}
		float max = (float)java.lang.Math.log(getMaxCount());
		int yold = 0;
		for (int i = 0; i < numberOfBins; i++){
			float yf = (float)java.lang.Math.log(bins[i]);
			yf/=(max/numberOfBins);
			int y = (int)yf;
			y = (y<0)?0:y;
			y = (y>=numberOfBins)?numberOfBins-1:y;
			if(y<yold) {
				for(int k=y;k<=yold;k++)
					px[k*numberOfBins+i] = -1;
			} else {
				for(int k=yold;k<=y;k++)
					px[k*numberOfBins+i] = -1;
			}
			yold = y;
		}
		ImagePlus impOut = new ImagePlus(WindowManager.makeUniqueName("SegmentHistogram"),ip);
		impOut.setProperty("AutoHistogram", this);
		impOut.show();
		return impOut;
	}
	public void displayInHistogram(ImagePlus imp) {
		AutoHistogram hist = (AutoHistogram)imp.getProperty("AutoHistogram");
		// put result in existing image
		ImageProcessor ip = imp.getProcessor();
		byte[] px = (byte[])ip.getPixels();
		float max = (float)java.lang.Math.log(hist.getMaxCount());
		int yold = 0;
		for (int i = 0; i < hist.numberOfBins; i++){
			float yf = (float)java.lang.Math.log(getBin(hist.lower+((float)i+0.5f)*hist.increment));
			yf/=(max/numberOfBins);
			int y = (int)yf;
			y = (y<0)?0:y;
			y = (y>=hist.numberOfBins)?hist.numberOfBins:y;
			if(y<yold) {
				for(int k=y;k<=yold;k++)
					px[k*numberOfBins+i] = 127;
			} else {
				for(int k=yold;k<=y;k++)
					px[k*numberOfBins+i] = 127;
			}
			yold = y;
		}
		imp.draw();
	}
}