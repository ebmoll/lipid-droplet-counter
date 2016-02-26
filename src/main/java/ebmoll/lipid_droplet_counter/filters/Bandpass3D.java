package ebmoll.lipid_droplet_counter.filters;

import ij.IJ;
//import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

public class Bandpass3D {
	public ImageStack			in_image;
	public double				in_hprad;
	public double				in_lprad;
	public double				in_xz_ratio;
	public ImageStack			out_result;
	
	public String getDescription() {
		return "Executes a Bandpass on the image." +
				"The Bandpass function is made of two Gaussians." +
				"xz_ratio is normally bigger than one (i.e. z" +
				"-Resolution is worse than x/y-resoution)";
	}
	
	public String checkInputParams() {
		String error = "";
		if(!(in_image.getProcessor(1) instanceof FloatProcessor))
			error += "only 32bit float images are supported.\n";
		if(in_hprad <= 0.0f)
			error += "hprad cannot be <= 0!\n";
		if(in_lprad <= 0.0f)
			error += "lprad cannot be <= 0!\n";
		if(in_lprad >= in_hprad)
			error += "hprad must be greater than lprad!\n";
		if(in_xz_ratio <= 0.0f)
			error += "xz_ratio cannot be <= 0!\n";
		return error;
	}
	
	public void filterit() {
		// create lowpass filter
		Convolve3D conv3d = new Convolve3D();
		conv3d.in_image = in_image;
		conv3d.in_kernel = createfilter();
		conv3d.in_normalize = false;
		//(new ImagePlus("blub",conv3d.in_kernel)).show();
		if(!conv3d.checkInputParams().equals("")) {
			IJ.showMessage(conv3d.checkInputParams());
			return;
		}
		conv3d.filterit();
		out_result = conv3d.out_result;
	}
	
	ImageStack createfilter() {
		int in_image_x = in_image.getWidth();
		int in_image_y = in_image.getHeight();
		int in_image_z = in_image.getSize();

		
		ImageStack out_gaussian = new ImageStack(in_image_x,in_image_y);
		
		int x_cent = in_image_x/2;
		int y_cent = in_image_y/2;
		int z_cent = in_image_z/2;
		
		double h_width = in_hprad*0.5;
		double h_height = in_hprad*0.5;
		double h_depth = in_hprad/in_xz_ratio*0.5;
		double l_width = in_lprad*0.5;
		double l_height = in_lprad*0.5;
		double l_depth = in_lprad/in_xz_ratio*0.5;
		
		// find high and low normalization factor
		double hsum = 0.0;
		double lsum = 0.0;
		for(int z=0;z<in_image_z;z++) {
			for(int x=0;x<in_image_x;x++) {
				for(int y=0;y<in_image_y;y++) {
					double hdist_sq = ((x-x_cent)/h_width)*((x-x_cent)/h_width) + ((y-y_cent)/h_height)*((y-y_cent)/h_height) + ((z-z_cent)/h_depth)*((z-z_cent)/h_depth);
					hsum += Math.pow(2.0,-hdist_sq);
					double ldist_sq = ((x-x_cent)/l_width)*((x-x_cent)/l_width) + ((y-y_cent)/l_height)*((y-y_cent)/l_height) + ((z-z_cent)/l_depth)*((z-z_cent)/l_depth);
					lsum += Math.pow(2.0,-ldist_sq);
				}
			}
		}
		// create image...
		for(int z=0;z<in_image_z;z++) {
			float[] arr = new float[in_image_x*in_image_y];
			for(int x=0;x<in_image_x;x++) {
				for(int y=0;y<in_image_y;y++) {
					double hdist_sq = ((x-x_cent)/h_width)*((x-x_cent)/h_width) + ((y-y_cent)/h_height)*((y-y_cent)/h_height) + ((z-z_cent)/h_depth)*((z-z_cent)/h_depth);
					double ldist_sq = ((x-x_cent)/l_width)*((x-x_cent)/l_width) + ((y-y_cent)/l_height)*((y-y_cent)/l_height) + ((z-z_cent)/l_depth)*((z-z_cent)/l_depth);
					// filter = low - high
					arr[y*in_image_x+x] = (float)((Math.pow(2.0,-ldist_sq)/lsum)-(Math.pow(2.0,-hdist_sq)/hsum));
					//arr[y*in_image_x+x] = (float)(Math.pow(2.0,-ldist_sq)/lsum)+(float)((1-Math.pow(2.0,-hdist_sq))/hsum);
				}
			}
			out_gaussian.addSlice(null,new FloatProcessor(in_image_x,in_image_y,arr,null));
		}
		return out_gaussian;
	}
}
