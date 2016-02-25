package ebmoll.lipid_droplet_counter.filters;

import java.util.*;

import ij.ImageStack;
import ij.process.*;

public class Watershed3D {
	public ImageStack		in_image;
	public ImageStack		out_transform;
	public float			in_radx = 2.0f;
	public float			in_rady = 2.0f;
	public float			in_radz = 2.0f;
	public boolean			in_invert = true;
	
	
	public void filterit() {
		// read image stack into float[][] array
		int dimx = in_image.getWidth();
		int dimy = in_image.getHeight();
		int dimz = in_image.getSize();
		float[][] data = new float[dimz][];
		if(in_image.getProcessor(1) instanceof FloatProcessor) {
			for(int i = 0; i < dimz; i++){
				data[i] = (float[])in_image.getProcessor(i+1).getPixels();
			}
		} else {
			for(int i = 0; i < dimz; i++){
				data[i] = (float[])in_image.getProcessor(i+1).convertToFloat().getPixels();
			}
		}
		
		int[][] transform = (new watershed()).watershed_transform(data,dimx,dimy,dimz,in_invert,in_radx,in_rady,in_radz);
		
		out_transform = new ImageStack(dimx,dimy);
		for (int k = 0; k < dimz; k++){
			ImageProcessor oip = new ColorProcessor(dimx,dimy,transform[k]);
			out_transform.addSlice(null,oip);
		}
	}

	public String getDescription() {
		return "Makes a watershed transform of the image. Takes as input any" +
				"greyscale image or image stack, and outputs a RGB color" +
				"image/stack where every color is a index for a region.";
	}

}

/*
 * This plugin implements the Watershed algorithm of Osma-Ruiz et al. as
 * described in "An improved watershed algorithm based on efficient computation
 * of shortest paths", applied to 3D-image stacks.
 * 
 * TODO:
 * BUGS:
 *   the plugin currently makes a float copy of the image stack for
 *   the transformation. This works, but uses memory and is not really
 *   necessary...
 * 
 * 
 * CHANGELOG:
 *      - Added "Invert" option
 *      - 2008-11-05: Output is now in 32-bit RGB-image format instead
 *        of 16-bit grayscale. Overflow therefore is no longer a problem.
 * 
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
 */


class watershed {
	
	public int[][] watershed_transform(float[][] stack, int dimx, int dimy, int dimz, boolean invert, float radx, float rady, float radz) {
		/*
		 * Does a watershed transform of stack.
		 * Returns an int stack of equal size where every point has a positive
		 * number according to the catchment basin. If invert is true,
		 * the image is processed inverted.
		 * 
		 * BUGS: may do something weird if the number of catchment basins
		 * exceeds 0x7FFFFFFF. This should never happen except with VERY
		 * large and noisy stacks. 
		 */
		/*
		 * TODO: maybe implement this also for integer/short arrays
		 * 
		 * TODO: optimize: maybe modify the algorithm not to flood-fill
		 * the whole plateau, but to fill it from the edge:
		 *  - ditch qInner
		 *  - create a qEdgeUp and qEdgeDown and process the plateau from
		 *    whichever Edge is encountered first. This should perform much
		 *    better with large plateaus (like completely uniform images)
		 *    
		 * TODO: optimize the algorithm not to need FIFO, so vectors or
		 * simple arrays can be used for queues instead of linked lists
		 */
		/*
		 * Explanation of the algorithm:
		 * loop through all voxels and let every voxel point to its steepest
		 *   descending neighbor. Label the voxel if it is a minimum.
		 * if the voxel belongs to a plateau (i.e. it has a neighbor with the
		 *   same value) then flood-fill the whole plateau and label every
		 *   point either as inner or edge point. Loop through all edge
		 *   points (this is where FIFO gets important) and let the
		 *   adjacent inner points point to the edge.
		 *   If there are no edge points, label the whole plateau as
		 *   a minimum. 
		 * now all voxels point to someone, and all minima are labeled, so
		 *   just follow the pointers to the minimum and assign it.
		 */
		
		// TODO: maybe use bigger values, because if NumCatchment reaches
		// one of these numbers, something goes wrong...
		final int UNVISITED = -0x01000000;
		final int PENDING = -0x02000000;
		
		// should the image be inverted?
		float mult = 1.0f;
		if(invert)
			mult = -1.0f;
		
		// create the four empty queues
		PosQueue qPending = new PosQueue();
		PosQueue qEdge = new PosQueue();
		PosQueue qInner = new PosQueue();
		PosQueue qDescending = new PosQueue();
		
		// TODO: this should really be an argument to the function, not read from
		// class variables
		ArrayList<int[]> neighborhood = generateNeighborhood(radx,rady,radz);

		// create and initialize transform array
		int[][] transform = new int[dimz][dimx*dimy];
		for(int d=0;d<dimz;d++){
			for(int i=0;i<dimx*dimy;i++){
				transform[d][i] = UNVISITED;
			}
		}
		int NumCatchment = 0;
		/* first pass: find minima and steepest descending paths;
		 * resolve plateaus
		 */
		for( int z = 0; z < dimz; z++) {
			//IJ.showProgress(z*3,dimz*4);
			for( int y = 0; y < dimy; y++) {
				for( int x = 0; x < dimx; x++) {
					// if the point has not been analyzed yet, process it
					if(transform[z][y*dimx+x]==UNVISITED) {
						float pVal = mult * stack[z][y*dimx+x];
						// loop through neighborhood
						float MinVal = pVal;
						int MinPointer = packpos(0,0,0);
						int[] Min = {x,y,z};
						ArrayList<Integer> N = getNeighborhood(x,y,z,dimx,dimy,dimz,neighborhood);
						for( int i=0 ; i < N.size() ; i++) {
							int pnPointer = N.get(i);
							int[] pnPos = unpackpos(pnPointer);
							pnPos[0]+=x;pnPos[1]+=y;pnPos[2]+=z;
							float pnVal = mult * stack[pnPos[2]][pnPos[1]*dimx+pnPos[0]];
							// put all plateau points in qPending
							if(pnVal == pVal) {
								if(qPending.isEmpty()) {
									int[] pPos = {x,y,z};
									qPending.put(pPos);
									transform[z][y*dimx+x] = PENDING;
								}
								qPending.put(pnPos);
								transform[pnPos[2]][pnPos[1]*dimx+pnPos[0]] = PENDING;
							}
							// find the minimum
							else if(pnVal < MinVal) {
								MinVal = pnVal;
								MinPointer = pnPointer;
								Min[0]=pnPos[0];Min[1]=pnPos[1];Min[2]=pnPos[2];
							}
						}
						// if the point doesn't belong to a plateau, point
						// it to Min, if Min does not exist, mark point as Min
						if( qPending.isEmpty() ) {
							if(Min[0]==x && Min[1]==y && Min[2]==z) {
								transform[z][y*dimx+x] = NumCatchment;
								NumCatchment++;
							}
							else {
								transform[z][y*dimx+x] = MinPointer;
							}
						}
						// point belongs to a plateau -> process the plateau
						else {
							// flood-fill the plateau
							while(!qPending.isEmpty()) {
								int[] pnPos = qPending.get();
								// [x,y,z] has already been processed
								if(pnPos[0]!=x || pnPos[1]!=y || pnPos[2]!=z) {
									// loop through neighborhood
									float pnVal = mult * stack[pnPos[2]][pnPos[1]*dimx+pnPos[0]];
									MinVal = pnVal;
									MinPointer = packpos(0,0,0);
									Min[0]=pnPos[0];Min[1]=pnPos[1];Min[2]=pnPos[2];
									N = getNeighborhood(pnPos[0],pnPos[1],pnPos[2],dimx,dimy,dimz,neighborhood);
									for( int i=0 ; i < N.size() ; i++) {
										int pnnPointer = N.get(i);
										int[] pnnPos = unpackpos(pnnPointer);
										pnnPos[0]+=pnPos[0];pnnPos[1]+=pnPos[1];pnnPos[2]+=pnPos[2];
										float pnnVal = mult * stack[pnnPos[2]][pnnPos[1]*dimx+pnnPos[0]];
										// put all plateau points in qPending
										if(pnnVal == pnVal) {
											if(transform[pnnPos[2]][pnnPos[1]*dimx+pnnPos[0]] == UNVISITED) {
												qPending.put(pnnPos);
												transform[pnnPos[2]][pnnPos[1]*dimx+pnnPos[0]] = PENDING;
											}
										}
										// find the minimum
										else if(pnnVal < MinVal) {
											MinVal = pnnVal;
											MinPointer = pnnPointer;
											Min[0]=pnnPos[0];Min[1]=pnnPos[1];Min[2]=pnnPos[2];
										}
									}
									
								}
								// classify point as inner or edge
								if(Min[0]==pnPos[0] && Min[1]==pnPos[1] && Min[2]==pnPos[2]) {
									qInner.put(pnPos);
								} else {
									qEdge.put(pnPos);
									transform[pnPos[2]][pnPos[1]*dimx+pnPos[0]] = MinPointer;
								}
							}
							// now qPending is empty, qInner and qEdge contain the plateau
							// no edge points -> mark plateau as minimum
							if(qEdge.isEmpty()) {
								while(!qInner.isEmpty()) {
									int[] pnPos = qInner.get();
									transform[pnPos[2]][pnPos[1]*dimx+pnPos[0]] = NumCatchment;
								}
								NumCatchment++;
							}
							// if there are edge points, recursively point
							// adjacent inner points to the edge
							else {
								if(!qInner.isEmpty()) {
									while(!qEdge.isEmpty()) {
										int[] pnPos = qEdge.get();
										// loop through neighborhood
										N = getNeighborhood(pnPos[0],pnPos[1],pnPos[2],dimx,dimy,dimz,neighborhood);
										for( int i=0 ; i < N.size() ; i++) {
											int pnnPointer = N.get(i);
											int[] pnnPosRel = unpackpos(pnnPointer);
											int[] pnnPos = pnnPosRel.clone();
											pnnPos[0]+=pnPos[0];pnnPos[1]+=pnPos[1];pnnPos[2]+=pnPos[2];
											// point all PENDING neighbors to pnPos
											if(transform[pnnPos[2]][pnnPos[1]*dimx+pnnPos[0]] == PENDING) {
												int pnPointer = packpos(-pnnPosRel[0],-pnnPosRel[1],-pnnPosRel[2]);
												transform[pnnPos[2]][pnnPos[1]*dimx+pnnPos[0]] = pnPointer;
												qEdge.put(pnnPos);
											}
										}
									}
									qInner.clear();
								} else {
									qEdge.clear();
								}
							}
						}
					}
				}
			}
		}
		/*
		 * second pass: resolve steepest descending paths and assign
		 * a catchment basin to every point
		 */
		for( int z = 0; z < dimz; z++) {
			//IJ.showProgress(z+dimz*3,dimz*4);
			for( int y = 0; y < dimy; y++) {
				for( int x = 0; x < dimx; x++) {
					int[] pnPos = {x,y,z};
					int pnVal = transform[z][y*dimx+x];
					// follow the Pointer to the minimum and put points along the
					// way into qDescending
					while(pnVal < 0) {
						qDescending.put(pnPos.clone());
						int[] pnnPos = unpackpos(pnVal);
						pnPos[0]+=pnnPos[0];pnPos[1]+=pnnPos[1];pnPos[2]+=pnnPos[2];
						pnVal = transform[pnPos[2]][pnPos[1]*dimx+pnPos[0]];
					}
					// mark all points in qDescending as minimum
					while(!qDescending.isEmpty()) {
						pnPos = qDescending.get();
						transform[pnPos[2]][pnPos[1]*dimx+pnPos[0]] = pnVal;
					}
				}
			}
		}
		// return the transform
		return transform;
	}
	
	private ArrayList<Integer> getNeighborhood(int x, int y, int z, int dimx, int dimy, int dimz, ArrayList<int[]> nbh) {
		/*
		 * returns a copy of nbh, where all points that lie outside the
		 * voxelstack are removed
		 */
		ArrayList<Integer> result = new ArrayList<Integer>();
		
		for(int c=0;c<nbh.size();c++) {
			int[] p = nbh.get(c);
			int i = p[0];
			int j = p[1];
			int k = p[2];
			if(x+i>=0 && x+i<dimx && y+j>=0 && y+j<dimy && z+k>=0 && z+k<dimz)
				result.add(packpos(i,j,k));
		}
		return result;
	}
	
	private ArrayList<int[]> generateNeighborhood(float rx, float ry, float rz) {
		/*
		 * Returns all positions that are contained within an ellipse of
		 * radii rx,ry,rz around the origin. The positions are ordered according
		 * to their distance from the origin
		 */
		//rx = (rx<1.0f) ? 1.0f : rx;
		//ry = (ry<1.0f) ? 1.0f : ry;
		//rz = (rz<1.0f) ? 1.0f : rz;
		ArrayList<int[]> result = new ArrayList<int[]>();
		int xmax = (int)Math.ceil(rx);
		int ymax = (int)Math.ceil(ry);
		int zmax = (int)Math.ceil(rz);
		for(int k=-zmax;k<=zmax;k++){
			for(int j=-ymax;j<=ymax;j++){
				for(int i=-xmax;i<=xmax;i++){
					if((i*i)/(rx*rx)+(j*j)/(ry*ry)+(k*k)/(rz*rz) <= 1.0f &&
							( i!=0 || j!=0 || k!=0 ) ) {
						int[] res = {i,j,k};
						result.add(res);
					}
				}
			}
		}
		Collections.sort(result,new DistanceComparator());
		return result;
	}
	
	class DistanceComparator implements Comparator<int[]>
	{
		public int compare(int[] x, int[] y)
		{
			float dx = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
			float dy = y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
			if(dx<dy)
				return -1;
			if(dx>dy)
				return 1;
			return 0;
		}
	}

	private int packpos(int x, int y, int z) {
		/*
		 * Encodes a relative position in a single negative value
		 */
		return -((x+128)+((y+128)<<8)+((z+128)<<16));
	}
	
	private int[] unpackpos(int value) {
		int[] pos = new int[3];
		int val = -value;
		pos[0] = (val&0xFF)-128;
		val >>>= 8;
		pos[1] = (val&0xFF)-128;
		val >>>= 8;
		pos[2] = (val&0xFF)-128;
		return pos;
	}

	private class PosQueue{
		// a FIFO queue that holds positions as int[3]
		private LinkedList<int[]> q;
		public PosQueue() {
			q = new LinkedList<int[]>();
		}
		public void put(int[] pos) {
			q.addLast(pos);
		}
		public int[] get() {
			return q.removeFirst();
		}
		public void clear() {
			q.clear();
		}
		public boolean isEmpty() {
			return q.isEmpty();
		}
		public int getLength() {
			return q.size();
		}
		public boolean contains(int[] Entry) {
			for(int i=0;i<q.size();i++) {
				int[] x = q.get(i);
				if(x[0]==Entry[0] && x[1]==Entry[1] && x[2]==Entry[2])
					return true;
			}
			return false;
		}
	}
}
