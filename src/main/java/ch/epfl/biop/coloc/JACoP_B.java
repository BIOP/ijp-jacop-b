package ch.epfl.biop.coloc;

/*JACoP: "Just Another Colocalization Plugin..." v1, 13/02/06
    Fabrice P Cordelieres, fabrice.cordelieres at curie.u-psud.fr
    Susanne Bolte, Susanne.bolte@isv.cnrs-gif.fr
 
    Copyright (C) 2006 Susanne Bolte & Fabrice P. Cordelieres
  
    License:
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *
 *
*/

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ch.epfl.biop.coloc.utils.ImageColocalizer;
import ch.epfl.biop.montage.StackMontage;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.MontageMaker;
import ij.plugin.PlugIn;
import ij.plugin.StackCombiner;
import ij.plugin.Thresholder;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

public class JACoP_B implements PlugIn {
    
	int channelA, channelB;
	
	ImagePlus impA, impB;
	
	String thrA, thrB;
	
	int mThrA, mThrB;
	
	int fluo_min=0;
	
	int fluo_max=255;
	
	Boolean doCostesThr=false, doPearsons=false, doOverlap=false, doManders=false, doFluorogram=false, doICA=false, doCostesRand=false;
	
	Calibration calib;
	
	int xyBlock, zBlock, nbRand, fitMeth;
	double binWidth;
	Boolean xyRand, zRand, showRand;
	
	File imageFolder = null;
	
	Boolean doSeparateZ;
	
	Boolean hasRoiSets = false;
	
	ImagePlus imp;

	private boolean is_stack_hist_z;
	@Override
	public void run(String arg) {
		
		if (!showDialog()) return;
	
		// Make some sense of it all
		
		//Switch from folder to single image mode
		if (imageFolder == null) {
			
			// Setup the image and all settings
			imp = IJ.getImage();
			runColoc(imp);
			
		}
	}

		
	private void runColoc(ImagePlus im) {
		// check dimensions and split in Z and T accordingly
		
		int nC = im.getNChannels();
		int nZ = im.getNSlices();
		int nT = im.getNFrames();
		
		// Check which thresholds to use
		if(thrA.matches("Use Manual Threshold Below")) thrA = String.valueOf(mThrA);
		if(thrB.matches("Use Manual Threshold Below")) thrB = String.valueOf(mThrB);
		
		if(thrA.matches("Costes Auto-Threshold") || thrA.matches("Costes Auto-Threshold"))  doCostesThr = true;

		
		Duplicator dup = new Duplicator();
		ImagePlus singleImp; 
		
		ImageColocalizer ic;
		
		// Store results
		List<ImagePlus> results = new ArrayList<ImagePlus>();
		
		// Get eventual ROIs
		RoiManager rm = RoiManager.getInstance2();
		if(rm == null) {
			rm = new RoiManager(false);
			IJ.log("The ROI Manager was empty");			
			
		}
		IJ.log("There are "+rm.getCount()+" ROIs available");

		if (rm.getCount() > 0) {
			hasRoiSets = true;
		}
		
		Roi roi = im.getRoi();

		int rcount = (rm.getCount()>0) ? rm.getCount() : 1;
		
		for(int r=0; r<rcount; r++) {
			if(hasRoiSets) {
				im.setRoi(rm.getRoi(r));
				roi = rm.getRoi(r);
				IJ.log("Selected ROI number "+r);

			}

			im.deleteRoi();
			
			if (!doSeparateZ) {
				// Work on just time
				for(int t=1; t<=nT; t++) {
					
					
					singleImp = dup.run(im, 1, nC, 1, nZ , t, t);
					singleImp.setTitle(im.getTitle());
					
					if(roi != null) {
						singleImp.setRoi(roi);

					}
					
					ic = new ImageColocalizer(singleImp, channelA, channelB);
					IJ.log("ThrA="+ thrA);
					IJ.log("ThrB="+ thrB);
					if(doCostesThr) {
						ic.CostesAutoThr();
					} else {
						ic.setThresholds(thrA, thrB);
					}
					
					// Do the analysis
					ImagePlus tmp = runAnalysis(ic);
					tmp.setTitle(tmp.getTitle()+" T"+t);
					results.add(tmp);
					ic.addResult("Time", t);
					ic.showResults();	

	
				}
				
			} else {
				
				for(int t=1; t<=nT; t++) {
					
					// Get the stack thresholds in case we need them
					singleImp = dup.run(im, 1, nC, 1, nZ , t, t);
					if(roi != null) {
						singleImp.setRoi(roi);
					}
										
					ic = new ImageColocalizer(singleImp, channelA, channelB);
					
					if(doCostesThr) {
						ic.CostesAutoThr();
					} else {
						ic.setThresholds(thrA, thrB);
					}
					int thrA_stack = ic.getThresholdA();
					int thrB_stack = ic.getThresholdB();
					ic.removeLastRow();
					
					ImageStack allZ = null;
					ImagePlus tmpImg = null;
					for(int z=1; z<=nZ; z++) {

						singleImp = dup.run(im, 1, nC, z, z , t, t);
						singleImp.setTitle(im.getTitle());
						
						if(roi != null) {
							singleImp.setRoi(roi);
						}
						
						ic = new ImageColocalizer(singleImp, channelA, channelB);
						
						if(is_stack_hist_z) {
							ic.setThresholds(thrA_stack, thrB_stack);
						} else {
							if(doCostesThr) {
								ic.CostesAutoThr();
							} else {
								ic.setThresholds(thrA, thrB);
							}
						}
						// Do the analysis
						
						tmpImg = runAnalysis(ic);
						
						if(z==1) {
							allZ = tmpImg.createEmptyStack();
						}
						
						allZ.addSlice(tmpImg.getProcessor());
						allZ.setSliceLabel(imp.getTitle()+" Z"+z+" T"+t, z);
						ic.addResult("Slice", z);
						ic.addResult("Time", t);
						
					}
					results.add(new ImagePlus(tmpImg.getTitle()+"Individual Z T"+t, allZ));
					ic.showResults();	
					

				}
				// Because here we have exactly one slide and one time, we can make them into a stack at the end

				
			}
			if(roi != null) {
				imp.setRoi(roi);
			}
		}
		
		// Show the results...
	 	for(ImagePlus i : results) {
	 		i.show();
	 	}
	 	

		
	 	
	}
	
	ImagePlus runAnalysis(ImageColocalizer ic) {
		if(doPearsons) ic.Pearson();
		if(doManders) ic.MM();
		if(doOverlap) ic.Overlap();
		if(doICA) ic.ICA();
		if(doFluorogram) ic.CytoFluo();
		
		// Always get the images and make a montage
		ArrayList<ImagePlus> imgs = new ArrayList<ImagePlus>();
		
		imgs.add(ic.getRGBImageA(false));
		imgs.add(ic.getRGBMaskA());
		
		imgs.add(ic.getRGBImageB(false));
		imgs.add(ic.getRGBMaskB());
		
		imgs.add(ic.getRGBColocImage());
		imgs.add(ic.getRGBANDMask());

        ImagePlus montage;
        
        
        if(imgs.get(0).getNSlices() == 1) {
        	ImageStack montagestk = imgs.get(0).createEmptyStack();
            // Make a normal montage
        	 	for(ImagePlus i : imgs) {
        	 		montagestk.addSlice(i.getProcessor().convertToRGB());
        	 	}
        	 	MontageMaker mm = new MontageMaker();
            	montage = mm.makeMontage2(new ImagePlus("for montage",montagestk), 2, 3, 1.0, 1, imgs.size(), 1, 0, false);

        } else {
        	montage = StackMontage.montageImages(imgs,3, 2);
        }
    	

		//Eventually add the fluorogram
		if(doFluorogram) {
			ImagePlus fluo = ic.getFluorogramImage(fluo_min, fluo_max);
			ImagePlus scaledFluo = scale(fluo, montage.getWidth());
			
			ImageStack flst = scaledFluo.getStack();
			// Make same number of dimensions as the montage (Z slices eventually)
			for(int i=1; i<montage.getStackSize(); i++) {
				flst.addSlice(scaledFluo.getProcessor().duplicate());
			}
			scaledFluo.setStack(flst);
			// Finally assemble them
			StackCombiner sc = new StackCombiner();
			ImageStack montagefluo = sc.combineVertically(montage.getStack(), scaledFluo.getStack());
			montage = new ImagePlus(imp.getTitle()+" Colocalization Report", montagefluo);
			//new ImagePlus("Test", montagefluo).show();
		}
		return montage;
	}

	

    private ImagePlus scale(ImagePlus fluo, int width) {
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
    
    
	private Boolean showDialog() {
		// If no images, dialog with folder otherwise normal dialog
		int nImages = WindowManager.getImageCount();
		
		// Make the GUI, old school
		GenericDialogPlus d = new GenericDialogPlus("Colocalization Parameters");
		
		if (nImages == 0) d.addDirectoryField("Image_Folder", "");
		d.addNumericField("Channel_A", 1, 0);
		d.addNumericField("Channel_B", 2, 0);
		
		List<String> thrs= new ArrayList<String>(Arrays.asList(Thresholder.methods));
		thrs.add("Costes Auto-Threshold");
		thrs.add("Use Manual Threshold Below");
		
		String[] thrsArr = new String[thrs.size()]; 
		thrsArr = thrs.toArray(thrsArr);
		
		d.addChoice("Threshold_for_Channel_A", thrsArr, thrsArr[0]);
		d.addChoice("Threshold_for_Channel_B", thrsArr, thrsArr[0]);
		
		d.addNumericField("Manual_Threshold_A", 0, 0);
		d.addNumericField("Manual_Threshold_B", 0, 0);

        
		d.addCheckbox("Consider_Z_Slices_Separately", false);
		d.addCheckbox("Set_Auto_Thresholds_On_Stack_Histogram", true);
		
		d.addMessage("Colocalization Result Options");
		
		d.addCheckbox("Get_Pearsons Correlation", true);
		d.addCheckbox("Get_Manders Coefficients", true);
		d.addCheckbox("Get_Overlap Coefficients", true);
		d.addCheckbox("Get_Li_ICA ", true);
		d.addCheckbox("Get_Fluorogram ", true);
		
		d.addCheckbox("Perform_Costes_Randomization (Not implemented)", false);

		
		//d.addChoiceMessage("Report Choice");
		
		d.showDialog();
		if(d.wasCanceled()) {
			return false;
		}
		
		// Pick up the data
		if (nImages == 0) imageFolder = new File(d.getNextString());
		channelA = (int) d.getNextNumber();
		channelB = (int) d.getNextNumber();
		thrA = d.getNextChoice();
		thrB = d.getNextChoice();
		
		mThrA = (int) d.getNextNumber();
		mThrB = (int) d.getNextNumber();
		
		doSeparateZ = d.getNextBoolean();
		is_stack_hist_z = d.getNextBoolean();
		doPearsons = d.getNextBoolean();
		doManders = d.getNextBoolean();
		doOverlap = d.getNextBoolean();
		doICA = d.getNextBoolean();
		doFluorogram = d.getNextBoolean();
		
		doCostesRand = d.getNextBoolean();
		
		IJ.log("ThrA="+ thrA);
		IJ.log("ThrB="+ thrB);
		return true;
    }
    
    
    public static void main(String[] args) {     
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = JACoP_B.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);
		
		ImageJ ij = new ImageJ();
		ij.exitWhenQuitting(true);
		
		// Make some nice images
		//ImagePlus imp = IJ.openImage("http://wsr.imagej.net/images/confocal-series.zip");
		//imp.show();
	//	ImagePlus imp = IJ.openImage("http://wsr.imagej.net/images/FluorescentCells.zip");
		ImagePlus imp = IJ.openImage("F:\\People\\Nadine Schmidt\\20170712_KN35_IF_A1_2_LUT_BC.tif");
		imp.show();
		IJ.openImage("F:\\People\\Nadine Schmidt\\ROI Sets\\20170712_KN35_IF_A1_2_LUT_BC.zip");
		
		//int[] xpoints = {166,150,258,253,212};
		//int[] ypoints = {256,208,137,244,291};
		//imp.setRoi(new PolygonRoi(xpoints,ypoints,5,Roi.POLYGON));

		IJ.run("JACoP_B", "");
    }
}



