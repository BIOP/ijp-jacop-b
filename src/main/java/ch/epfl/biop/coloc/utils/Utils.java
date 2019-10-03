package ch.epfl.biop.coloc.utils;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.plugin.Duplicator;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
/*
 * This class contains several utilities needed to manage images for colocalization
 * 
 * These are either features that have issues in Fiji and just need to be self-implemented
 * Or just convenience functions
 */
public class Utils {

	/*
	 * Reimplemented clearOutside as ip.fillOutside(roi) is causing a shift of (1,1) of the ROI
	 * when running it as a plugin versus running it as a macro or GUI command...
	 * This in turn causes a huge error with threshold methods later, because there is still background in the ROI
	 */
	public static synchronized void clearOutside(ImageProcessor ip, ImageProcessor mask) {
		for(int i=0; i<ip.getWidth(); i++) {
			for(int j=0; j<ip.getHeight(); j++) {
				ip.setf(i, j, (mask.getf(i, j) > 0 ? ip.getf(i, j) : 0));

			}
		}
	}

	/*
	 * Get a clean binary mask from a given threshold value
	 */
	public static ImageProcessor binarize(ImageProcessor ip, int lowerThr) {
		ByteProcessor bp = new ByteProcessor(ip.getWidth(), ip.getHeight());
		for(int x=0; x<ip.getWidth(); x++) {
			for(int y=0; y<ip.getHeight(); y++) {
				if (ip.getf(x, y) >= lowerThr) {
					bp.set(x,y, 255);
				}
			}
		}
		return bp;
	}

	/*
	 * Binarize for image stacks 
	 */
	public static ImagePlus binarize(ImagePlus imp, int lowerThr) { 

		if (imp.getStackSize() == 1) {
			ImagePlus imp2 = new ImagePlus("Mask "+imp.getTitle(), Utils.binarize(imp.getProcessor(), lowerThr));
			imp2.setCalibration(imp.getCalibration());
			return imp2;
		}

		// Handle stack
		ImageProcessor ip;
		ImageStack stack1 = imp.getStack();
		ImageStack stack2 = imp.createEmptyStack();
		int nSlices = imp.getStackSize();
		String label;

		for(int i=1; i<=nSlices; i++) {
			label = 
					stack1.getSliceLabel(i);
			ip = stack1.getProcessor(i);
			stack2.addSlice(label, Utils.binarize(ip, lowerThr));
		}

		ImagePlus imp2 = new ImagePlus("Mask "+imp.getTitle(), stack2);
		imp2.setCalibration(imp.getCalibration()); //update calibration

		return imp2;
	}

	/*
	 * Because images get their display range reset sometimes when manipulating them
	 * it is useful to have a way to get and reuse the display range
	 */
	public static double[][] getDisplayRange(ImagePlus imp) {
		double[][] drange = new double[imp.getNChannels()][];

		for(int i=0;i<drange.length;i++) {
			imp.setC(i+1);
			drange[i] = new double[2];
			drange[i][0] = imp.getDisplayRangeMin();
			drange[i][1] = imp.getDisplayRangeMax();
		}
		return drange;
	}

	public static  void setDisplayRange(ImagePlus impr, double[][] range) {
		// TODO Auto-generated method stub
		for(int i=0; i<range.length;i++) {
			impr.setC(i+1);
			impr.setDisplayRange(range[i][0], range[i][1]);
		}
	}
	
	/*
	 * Rescaling a stack or a single slice
	 */
    public static ImagePlus scale(ImagePlus fluo, int width) {
    	if(fluo.getStackSize() == 1) {
			ImageProcessor ip = fluo.getProcessor();
		    ip.setInterpolationMethod(ImageProcessor.BILINEAR);

		    return new ImagePlus(fluo.getTitle(), ip.resize(width));
		}
		
		
		ImageStack fscaledStack = fluo.createEmptyStack();
		
		for(int i=0; i<fluo.getStackSize(); i++) {
			ImageProcessor ip = fluo.getProcessor();
			fscaledStack.addSlice(ip.resize(width));
		}
		return new ImagePlus(fluo.getTitle(), fscaledStack);
	}
    
	
	/* 
	 * Convenience method to crop slices as a bunch of ImagePluses
	 * Note that unlike cropTime, we don't need to handle ROIs here
	 */
	public static ImagePlus[] cropSlices(ImagePlus singleTimeImp, Roi roi, Boolean doSeparateZ) {
		
		int nC = singleTimeImp.getNChannels();
		int nZ = singleTimeImp.getNSlices();

		Duplicator dup = new Duplicator();
		
		if(doSeparateZ) {
			ImagePlus[] zSlices = new ImagePlus[nZ];
			for(int z=1; z<=nZ; z++) {
				zSlices[z-1] = dup.run(singleTimeImp, 1, nC, z, z, 1,1);
				zSlices[z-1].setTitle(singleTimeImp.getTitle()+" Z"+z);
				if(roi !=null) zSlices[z-1].setRoi(roi);
			}

			return zSlices;
		}

		else {
			ImagePlus[] allSlices = new ImagePlus[1];
			allSlices[0] = singleTimeImp;
			allSlices[0].setTitle(singleTimeImp.getTitle());
			if(roi !=null) allSlices[0].setRoi(roi);

			return allSlices;
		}
	}
	
	/*
	 * Convenience method to crop timepoints
	 * Ensures ROIs are properly cared for
	 */
	public static ImagePlus cropTime(ImagePlus imp, Roi roi, int timepoint, Boolean is_crop) {

		
		int nC = imp.getNChannels();
		int nZ = imp.getNSlices();
		
		// We are working on imp, which is unchanged
		Duplicator dup = new Duplicator();

		
		// if we want to crop the ROI, do it here
		if(is_crop) {
			imp.setRoi(roi);
		} else {
			imp.killRoi();
		}

		ImagePlus tmp = dup.run(imp, 1, nC, 1, nZ , timepoint, timepoint);			
		tmp.setTitle(imp.getTitle()+" T"+timepoint);
		// If there was a roi, we should update it
		if(roi!= null && is_crop) {
			Roi roi2 = (Roi)roi.clone();
			roi2.setLocation(0, 0);
			tmp.setRoi(roi2);
			roi = roi2;
		}
		return tmp;
	}
	/*
	 * Because getMask() returns a cropped version, we need to implement our own...
	 */
	public static ImageProcessor getMask(ImagePlus imp, Roi roi) {
        ImageProcessor ip = new ByteProcessor(imp.getWidth(), imp.getHeight());
        ip.setRoi(roi);
        ip.setValue(255);
        ip.fill(ip.getMask());
        return ip;
	}
}

