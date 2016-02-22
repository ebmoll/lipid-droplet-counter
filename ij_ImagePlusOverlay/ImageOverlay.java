package ij_ImagePlusOverlay;

import ij.ImagePlus;
import java.awt.*;
import java.awt.image.*;

/*
 * This is a hacked-up class that extends ImagePlus with a very useful overlay,
 * derived from Peter Sebastian Masny's ImageJ-Branch. I modified some stuff
 * that seemed to be broken or incorrect on his website. For the original
 * visit: http://www.masny.dk/imagej/overlay/index.php
 * 
 * I release this code under the GNU GPL. This code contains modified parts of
 * Peter Sebastian Masny's original code, which doesn't explicitly state a
 * license, but which I assume is released under GNU GPL, because ImageJ is,
 * and he compiled and released a version with his changes.
 * 
 * Copyright 2008 Samuel Moll (moll@biochem.mpg.de)
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


/** Other data structures might be more space efficient, but
 * for simplicity with draw I used a standard 32-bit structure */
public class ImageOverlay {
    
    boolean display = false;
    protected int width, height;
    protected ImagePlus imp;
    protected MemoryImageSource source;
    protected ColorModel cm;
    protected Image img;
    protected int[] pixels; // 32-bit
    
    
    public ImageOverlay(ImagePlus imp) {
        this.imp = imp;
        width = imp.getWidth();
        height = imp.getHeight();
        cm = ColorModel.getRGBdefault();
        pixels = new int[width*height];
        source = new MemoryImageSource(width, height, cm, pixels, 0, width);
        source.setAnimated(true);
        source.setFullBufferUpdates(true);
        img = Toolkit.getDefaultToolkit().createImage(source);
        reset();
        show();
    }
    
    public void show() {
        display = true;
    }
    
    public void hide() {
        display = false;
    }

 /*   public void setPixel(int x, int y, int i) {
        if ((x>=0) && (y>=0) && (x=0?dx:-dx;
        int absdy = dy>=0?dy:-dy;
        int n = absdy>absdx?absdy:absdx;
        double xinc = (double)dx/n;
        double yinc = (double)dy/n;
        n++;
        double x = x1<0?x1-0.5:x1+0.5;
        double y = y1<0?y1-0.5:y1+0.5;
        do {
            setPixel((int)x, (int)y,i);
            x += xinc;
            y += yinc;
        } while (--n>0);
    }*/
    
    public void reset(){
        img = Toolkit.getDefaultToolkit().createImage(source);
        if(imp.getWindow()!=null)
        	if(imp.getWindow().getCanvas()!=null)
        		imp.getWindow().getCanvas().repaint();
        //Necessary to update from within a mouse event, but seems wrong
    }
    
    public int[] getPixels() {
    	return pixels;
    }
    
    public int getPixel(int x, int y){
        return pixels[x+ width*y];
    }
    public Image getImage() {
        return img;
    }
}

