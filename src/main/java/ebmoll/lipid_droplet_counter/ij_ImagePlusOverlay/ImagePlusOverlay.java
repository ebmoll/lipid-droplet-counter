package ebmoll.lipid_droplet_counter.ij_ImagePlusOverlay;

import ij.*;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import java.awt.*;
//import java.awt.image.*;
import ij.gui.*;
import ij.macro.Interpreter;


/*
 * This is a hacked-up class that extends ImagePlus with a very useful overlay,
 * derived from Peter Sebastian Masny's ImageJ-Branch. I modified some stuff
 * that seemed to be broken or incorrect on his website. For the original
 * visit: http://www.masny.dk/imagej/overlay/index.php
 * 
 * This class is useful if you want to provide a preview window for some
 * filter that needs tweaking parameters. It can be used everywhere where 
 * you normally would use the ImagePlus.getCanvas().setDisplayList()
 * interface, which is insufficient for all but the most simple things.
 * 
 * The ImagePlusOverlay supports color with alpha, meaning you can draw
 * semi-transparent markers in different colors, or whatever you want.
 * 
 * I release this code under the GNU GPL. This code contains modified parts of
 * Peter Sebastian Masny's original code, which doesn't explicitly state a
 * license, but which I assume is released under GNU GPL, because ImageJ is,
 * and he compiled and released a version with his changes.
 * 
 * How to use this class:
 * ----------------------
 * 
 * import ij_ImagePlusOverlay.*;
 * 
 * [...]
 * 
 * dimx = 100;
 * dimy = 150;
 * image = new ImagePlusOverlay("Cool Image With Overlay", new FloatProcessor(dimx,dimy));
 * ImageOverlay overlay = new ImageOverlay(image);
 * overlay.show(); // This needs to be called because the overlay is invisible by default
 * image.setOverlay(overlay);
 * image.show();
 * // How to draw a semitransparent green pixel at [42,21] on the overlay:
 * oi.getPixels()[21*dimx+42] = 0x4400FF00;
 * oi.reset(); // shows the changes
 * 
 * 
 * 
 * Copyright 2009 Samuel Moll (moll@biochem.mpg.de)
 * Max-Planck-Institute of Biochemistry
 * Research Group for Organelle Architecture and Dynamics / Tobias Walther
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */


public class ImagePlusOverlay extends ImagePlus {

	private ImageJ ij = IJ.getInstance();
	private boolean activated;
	private ImageOverlay overlay;


	public ImagePlusOverlay(String arg0) {
		super(arg0);
		// TODO Auto-generated constructor stub
	}

	public ImagePlusOverlay(String arg0, Image arg1) {
		super(arg0, arg1);
		// TODO Auto-generated constructor stub
	}

	public ImagePlusOverlay(String arg0, ImageProcessor arg1) {
		super(arg0, arg1);
		// TODO Auto-generated constructor stub
	}

	public ImagePlusOverlay(String arg0, ImageStack arg1) {
		super(arg0, arg1);
		// TODO Auto-generated constructor stub
	}

	public ImageOverlay getImageOverlay() {
		return overlay;
	}
	public void setOverlay(ImageOverlay overlay) {
		// TODO: sanity checks for size, etc
		this.overlay = overlay;
	}

	public void show(String statusMessage) {
		if (win!=null) return;
		if ((IJ.isMacro() && ij==null) || Interpreter.isBatchMode()) {
			ImagePlus img = WindowManager.getCurrentImage();
			if (img!=null) img.saveRoi();
			WindowManager.setTempCurrentImage(this);
			Interpreter.addBatchModeImage(this);
			return;
		}
		// TODO: This is just fucking stupid... the class ImagePlus is probably broken
		//if (Prefs.useInvertingLut && getBitDepth()==8 && ip!=null && !ip.isInvertedLut()&& !ip.isColorLut())
		//	invertLookupTable();
		img = getImage();
		if ((img!=null) && (width>=0) && (height>=0)) {
			activated = false;
			int stackSize = getStackSize();
			//if (compositeImage) stackSize /= nChannels;
			ImageCanvasOverlay ico = new ImageCanvasOverlay(this);
			if (stackSize>1)
				win = new StackWindow(this, ico);
			else
				win = new ImageWindow(this, ico);
			if (roi!=null) roi.setImage(this);
			draw();
			IJ.showStatus(statusMessage);
			if (IJ.isMacro()) { // wait for window to be activated
				long start = System.currentTimeMillis();
				while (!activated) {
					IJ.wait(5);
					if ((System.currentTimeMillis()-start)>2000) {
						WindowManager.setTempCurrentImage(this);
						break; // 2 second timeout
					}
				}
			}
			notifyListeners(OPENED);
		}
	}
	
	// TODO: This is also very stupid...
	public void setActivated() {
		activated = true;
	}

}


/** This is a Canvas used to display images in a Window.
 * It is modified to also draw an Overlay */
class ImageCanvasOverlay extends ImageCanvas {
	
	// this is sooo stupid...
	static final long		serialVersionUID = 23562357788L;
	
	public ImageCanvasOverlay(ImagePlusOverlay imp) {
		super(imp);
	}
	
	public void paint(Graphics g) {
		ImageOverlay overlay = ((ImagePlusOverlay)imp).getImageOverlay();

		super.paint(g);
		try {
			if (overlay != null) {
				if (overlay.display)
					g.drawImage(overlay.getImage(), 0,0,
							(int)(srcRect.width*magnification),
							(int)(srcRect.height*magnification),
							srcRect.x, srcRect.y,srcRect.x+srcRect.width,
							srcRect.y+srcRect.height, null);
			}
		}
		catch(OutOfMemoryError e) {IJ.outOfMemory("Paint");}

	}
}


