package ch.epfl.biop.coloc.utils;

import java.awt.Color;

import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.Thresholder;
import ij.process.Blitter;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.ShortProcessor;

public class ColocOutput {

	/*
	 * This class should pick up the results from a coloc experiment and output pretty images
	 * 
	 * Images to output v0.1
	 * 
	 * - RGB ImageA Original with current brightness and contrast settings and ROI overlay
	 * - RGB ImageB Original with current brightness and contrast settings and ROI overlay
	 * - Mask of Thresholded ImageA with ROI
	 * - Mask of Thresholded ImageB with ROI
	 * - Fluorogram of Channels A/B, showing thresholds and regression
	 * - LI-ICQ plots
	 * 
	 * Offer possibility to montage whatever we need
	 */
	
	

	public static ImagePlus Fluorogram(Plot fp, ImagePlus impA, ImagePlus impB, int min, int max, int thrA, int thrB) {

		float[] valA = fp.getXValues();
		float[] valB = fp.getYValues();
		
		int nbins = 256;
		int lut_size = 15;
		// Keep number of bins of 256 and we can scale the image as needed later
		
		// Processor with the bins
		ImageProcessor fluo_ip = new FloatProcessor(nbins, nbins);
		fluo_ip.set(0.0);
		// Processor with the X Scale
		ImageProcessor grad_A = new ShortProcessor(nbins, lut_size);
		
		// Processor with the Y Scale
		ImageProcessor grad_B = new ShortProcessor(lut_size, nbins);
		
		// Final Color processor with the fluorogram
		ImageProcessor fluo = new ColorProcessor(nbins+lut_size+1, nbins+lut_size+1);
		
		//Normalize the values and bin them
		for (int i=0; i<valA.length; i++) {
			
			int binA = (int) Math.round( (((double) valA[i] - (double) min) / (double) max * (double) (nbins-1)));
			int binB = (int) Math.round( (((double) valB[i] - (double) min) / (double) max * (double) (nbins-1)));
			//IJ.log("A:"+binA+" B:"+binB);
			
			if (binA < nbins && binB < nbins && binA >= 0 && binB >= 0) {
				fluo_ip.setf(binA, nbins-binB-1, fluo_ip.getPixelValue(binA, nbins-binB-1)+1);
			}
		}
		//Log of the values
		fluo_ip.log();
		
		// Build the grayLevels
		for(int i=0; i<nbins; i++) {
			for(int j=0; j<lut_size; j++) {
				grad_A.set(i,j,i);
				grad_B.set(j,nbins-i-1,i);
			}
		}
		
		grad_A.setLut(impA.getLuts()[0]);
		grad_B.setLut(impB.getLuts()[0]);
		
		// Somehow load the Fire LUT
		fluo_ip.setLut(fireLUT());
		fluo_ip.setMinAndMax(0, 6);
		
		fluo_ip.convertToRGB();
		
		//fluo_ip.setColor(impA.getLuts()[0].getRGB(0));
		fluo_ip.setColor(new Color(255,255,255));
		fluo_ip.drawLine(thrA, nbins, thrA, 0);

		//fluo_ip.setColor(impB.getLuts()[0].getColorModel().);
		fluo_ip.drawLine(0, nbins-thrB-1, nbins, nbins-thrB-1);
		
		// Build fluorogram
		fluo.copyBits(grad_A.convertToRGB(), lut_size+1, nbins+1, Blitter.ADD);
		fluo.copyBits(grad_B.convertToRGB(), 0, 0, Blitter.ADD);
		fluo.copyBits(fluo_ip, lut_size+1, 0, Blitter.ADD);
		
		// Add threshold lines
		
		return new ImagePlus(fp.getTitle(), fluo);
	}
	
	public static ImagePlus[] getImageMasks(ImagePlus impA, ImagePlus impB, Roi roi, int min, int max, int thrA, int thrB) {
		
		// Make masks from the images and make an image that also contains the overlap. Careful to add the ROI
		Thresholder t = new Thresholder();
		
		
		return null;
	}
		
		public static LUT fireLUT() {
	        int[] r = {0,0,1,25,49,73,98,122,146,162,173,184,195,207,217,229,240,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
	        int[] g = {0,0,0,0,0,0,0,0,0,0,0,0,0,14,35,57,79,101,117,133,147,161,175,190,205,219,234,248,255,255,255,255};
	        int[] b = {0,61,96,130,165,192,220,227,210,181,151,122,93,64,35,5,0,0,0,0,0,0,0,0,0,0,0,35,98,160,223,255};
			
	        byte[] reds   = new byte[256];
	        byte[] greens = new byte[256];
	        byte[] blues  = new byte[256];
	        
	        for (int i=0; i<r.length; i++) {
	            reds[i] = (byte)r[i];
	            greens[i] = (byte)g[i];
	            blues[i] = (byte)b[i];
	        }
	        interpolate(reds, greens, blues, r.length);
	        return new LUT(reds, greens, blues);
		}
		
		private static void interpolate(byte[] reds, byte[] greens, byte[] blues, int nColors) {
	        byte[] r = new byte[nColors]; 
	        byte[] g = new byte[nColors]; 
	        byte[] b = new byte[nColors];
	        System.arraycopy(reds, 0, r, 0, nColors);
	        System.arraycopy(greens, 0, g, 0, nColors);
	        System.arraycopy(blues, 0, b, 0, nColors);
	        double scale = nColors/256.0;
	        int i1, i2;
	        double fraction;
	        for (int i=0; i<256; i++) {
	            i1 = (int)(i*scale);
	            i2 = i1+1;
	            if (i2==nColors) i2 = nColors-1;
	            fraction = i*scale - i1;
	            //IJ.write(i+" "+i1+" "+i2+" "+fraction);
	            reds[i] = (byte)((1.0-fraction)*(r[i1]&255) + fraction*(r[i2]&255));
	            greens[i] = (byte)((1.0-fraction)*(g[i1]&255) + fraction*(g[i2]&255));
	            blues[i] = (byte)((1.0-fraction)*(b[i1]&255) + fraction*(b[i2]&255));
	        }
	    }
		
		
}