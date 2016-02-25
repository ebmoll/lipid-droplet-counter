package ebmoll.lipid_droplet_counter.filters;

import ij.ImageStack;
import ij.process.FloatProcessor;


public class Gaussian3D {
	public float			in_width;
	public float			in_height;
	public float			in_depth;
	public int				in_image_x;
	public int				in_image_y;
	public int				in_image_z;
	public ImageStack		out_gaussian;
	
	public String getDescription() {
		return "Calculates a centered Gaussian of the specified" +
				"FWHM width, height and depth. The resulting image size" +
				"can also be chosen. Creates a 32bit greyscale image. " +
				"If image x,y or z is 0, then image size will be chosen" +
				"automatically based on gaussian size. There will always" +
				"be one center pixel with intensity 1.0.";
	}
	
	public String checkInputParams() {
		String error = "";
		
		if(in_width <= 0.0f)
			error += "width cannot be <= 0!\n";
		if(in_height <= 0.0f)
			error += "height cannot be <= 0!\n";
		if(in_depth <= 0.0f)
			in_depth = 0.000001f;
		if(in_image_x <= 0)
			in_image_x = new Float(in_width * 4.5).intValue();
		if(in_image_y <= 0)
			in_image_y = new Float(in_height * 4.5).intValue();
		if(in_image_z <= 0)
			in_image_z = new Float(in_depth * 4.5).intValue();
		// all dimensions must be at least 1
		in_image_x = (in_image_x <= 0) ? 1 : in_image_x;
		in_image_y = (in_image_y <= 0) ? 1 : in_image_y;
		in_image_z = (in_image_z <= 0) ? 1 : in_image_z;
		
		return error;
	}
	
	public void filterit() {
		out_gaussian = new ImageStack(in_image_x,in_image_y);
		
		int x_cent = in_image_x/2;
		int y_cent = in_image_y/2;
		int z_cent = in_image_z/2;
		
		for(int z=0;z<in_image_z;z++) {
			// use a float[] array instead!!!!
			float[] arr = new float[in_image_x*in_image_y];
			for(int x=0;x<in_image_x;x++) {
				for(int y=0;y<in_image_y;y++) {
					double dist_sq = ((x-x_cent)/in_width)*((x-x_cent)/in_width) + ((y-y_cent)/in_height)*((y-y_cent)/in_height) + ((z-z_cent)/in_depth)*((z-z_cent)/in_depth);
					arr[y*in_image_x+x] = (float)Math.pow(2.0,-dist_sq);
				}
			}
			out_gaussian.addSlice(null,new FloatProcessor(in_image_x,in_image_y,arr,null));
		}
	}
}