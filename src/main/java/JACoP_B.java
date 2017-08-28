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
	Boolean doCropRois;
	
	Boolean hasRoiSets = false;
	
	
	ImagePlus imp=null;

	private boolean is_stack_hist_z;

	private Roi roi;

	private String imageName;

	private boolean is_montage_vertical;
	
	
	@Override
	public void run(String arg) {
		
		if (!showDialog()) return;
	
		// Make some sense of it all
		
		//Switch from folder to single image mode			
		// Setup the image and all settings
		this.imp = IJ.getImage();
		runColoc();
			
	}
	

		
	private void runColoc() {
		// check dimensions and split in Z and T accordingly
		
		int nC = imp.getNChannels();
		int nZ = imp.getNSlices();
		int nT = imp.getNFrames();
		
		// Check which thresholds to use
		if(thrA.matches("Use Manual Threshold Below")) thrA = String.valueOf(mThrA);
		if(thrB.matches("Use Manual Threshold Below")) thrB = String.valueOf(mThrB);
		
		if(thrA.matches("Costes Auto-Threshold") || thrA.matches("Costes Auto-Threshold"))  doCostesThr = true;

		
		// He who makes the magic happen
		ImageColocalizer ic;
		
		// Store results
		List<ImagePlus> results = new ArrayList<ImagePlus>();
		
		
		// Get eventual ROIs from the Manager
		RoiManager rm = RoiManager.getInstance();
		if(rm == null) {
			rm = new RoiManager(false);
		}
		if (rm.getCount() > 0) {
			hasRoiSets = true;
		}
		
		// 	The ROI we need to haul around in case it's there
		if(imp.getRoi() != null)
			roi = (Roi) imp.getRoi().clone();
		
		int rcount = (rm.getCount()>0) ? rm.getCount() : 1;
		
		// Master loop is for all ROIs
		for(int r=0; r<rcount; r++) {
			ImagePlus singleTimeImp;
			// Choose the ROI here
			if(hasRoiSets) {
				roi = (Roi) rm.getRoi(r).clone();
			}
			String roiName = null;
			
			imp.setRoi(roi);
			
			// Get the name of the ROI	
			if (roi != null) {
				roiName = roi.getName();
				if (roiName == null) {
	    			roiName = "ROI";
	    		}
			}
						
			imageName = imp.getTitle()+" ("+roiName+")";
			
			
			// Loop through time
			for(int t=1; t<=nT; t++) {
				singleTimeImp = cropTime(t, doCropRois);
				
				// Do Z separately, maybe
				ImagePlus[] zImages = cropSlices(singleTimeImp, doSeparateZ);
				int thrAH = 0;
				int thrBH = 0;
				// Make sure to get the right thresholds if we are working on the stack histogram...
				if(is_stack_hist_z) {
					// Need to run the colocalizer to get the thresholds
					ic = new ImageColocalizer(singleTimeImp, channelA, channelB);
					if(doCostesThr) {
						ic.CostesAutoThr();
					} else {
						ic.setThresholds(thrA, thrB);
					}
					
					// Save the thresholds for use below
					thrAH = ic.getThresholdA();
					thrBH = ic.getThresholdB();
					
					// Because initializing the colocalizer writes a result already,
					// We need to remove this ghost value
					ic.removeLastRow();
					
				}
				
				ImageStack res = null;
				
				// Iterate through all the images and do the actual work
				ImagePlus tmp = null;
				
				for(int i=0; i< zImages.length;i++) {
					if(roi != null) {
						roi.setPosition(0);
						zImages[i].setRoi(roi);
					}
					
					ic = new ImageColocalizer( zImages[i], channelA, channelB);
					
					if(doCostesThr) {
						ic.CostesAutoThr();
					} else {
						ic.setThresholds(thrA, thrB);
					}
					
					ic.addResult("Using Stack Histogram", "False");
					
					// Set the thresholds again in case they are to be made on the stack
					if( is_stack_hist_z) {
						ic.setThresholds(thrAH, thrBH);

						// Let the user know we used the stack histogram
						ic.addResult("Using Stack Histogram", "True");
					}
					// Do the analysis
					tmp = runAnalysis(ic);
					if( i==0 ) {
						res = tmp.createEmptyStack();
					}
					
					if(zImages.length > 1) { 
						res.addSlice(tmp.getProcessor());
					} else {
						res = tmp.getStack();
					}
					//res.setSliceLabel(tmp.getTitle(), res.getSize()-1);
		
					ic.addResult("Timepoint", t);
					if(zImages.length > 1) ic.addResult("Slice", i+1);
					
					ic.showResults();
					
				}
				tmp = new ImagePlus(tmp.getTitle(), res);
				
				results.add(tmp);
			}
		}
		
	 	for(ImagePlus i : results) {
	 		i.show();
	 	}
	}
	
	private ImagePlus[] cropSlices(ImagePlus singleTimeImp, Boolean doSeparateZ) {
		
		int nC = imp.getNChannels();
		int nZ = imp.getNSlices();

		Duplicator dup = new Duplicator();

		if(doSeparateZ) {
			ImagePlus[] zSlices = new ImagePlus[nZ];
			for(int z=1; z<=nZ; z++) {
				zSlices[z-1] = dup.run(singleTimeImp, 1, nC, z, z, 1,1);
				zSlices[z-1].setTitle(singleTimeImp.getTitle()+" Z"+z);

			}

			return zSlices;
		}

		else {
			ImagePlus[] allSlices = new ImagePlus[1];
			allSlices[0] = singleTimeImp;
			allSlices[0].setTitle(singleTimeImp.getTitle());

			return allSlices;
		}
	}


	private ImagePlus cropTime(int timepoint, Boolean is_crop) {
		
		
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
		if(is_crop && roi!= null) {
			roi.setLocation(0, 0);
			tmp.setRoi(roi);
		}		
		return tmp;
	}


	ImagePlus runAnalysis(ImageColocalizer ic) {
		if(doPearsons) ic.Pearson();
		if(doManders) ic.MM();
		if(doOverlap) ic.Overlap();
		if(doICA) ic.ICA();
		if(doFluorogram) ic.CytoFluo();
		int rows = 3;
		int columns = 2;
		// Always get the images and make a montage
		ArrayList<ImagePlus> imgs = new ArrayList<ImagePlus>();
        if(is_montage_vertical) {
			imgs.add(ic.getRGBImageA(false));
			imgs.add(ic.getRGBMaskA());
			
			imgs.add(ic.getRGBImageB(false));
			imgs.add(ic.getRGBMaskB());
			
			imgs.add(ic.getRGBColocImage());
			imgs.add(ic.getRGBANDMask());
        } else {
			imgs.add(ic.getRGBImageA(false));
			imgs.add(ic.getRGBImageB(false));
			imgs.add(ic.getRGBColocImage());

			imgs.add(ic.getRGBMaskA());
			imgs.add(ic.getRGBMaskB());
			imgs.add(ic.getRGBANDMask());
			columns = 3;
			rows = 2;
			
        }
        ImagePlus montage;
        
        // Make montage either vertical or horizontal

        
        
        if(imgs.get(0).getNSlices() == 1) {
        	ImageStack montagestk = imgs.get(0).createEmptyStack();
            // Make a normal montage
        	 	for(ImagePlus i : imgs) {
        	 		montagestk.addSlice(i.getProcessor().convertToRGB());
        	 	}
        	 	MontageMaker mm = new MontageMaker();
            	montage = mm.makeMontage2(new ImagePlus("for montage",montagestk), columns, rows, 1.0, 1, imgs.size(), 1, 0, false);

        } else {
        	montage = StackMontage.montageImages(imgs,rows, columns);
        }
    	
        

		//Eventually add the fluorogram
		if(doFluorogram) {
			ImagePlus fluo = ic.getFluorogramImage(fluo_min, fluo_max);
			
			ImagePlus scaledFluo = null;
			
	        if(is_montage_vertical) {
	        	scaledFluo = scale(fluo, montage.getWidth());
	        } else {
	        	scaledFluo = scale(fluo, montage.getHeight());

	        }
			ImageStack flst = scaledFluo.getStack();
			// Make same number of dimensions as the montage (Z slices eventually)
			for(int i=1; i<montage.getStackSize(); i++) {
				flst.addSlice(scaledFluo.getProcessor().duplicate());
			}
			scaledFluo.setStack(flst);
			// Finally assemble them
			StackCombiner sc = new StackCombiner();
			ImageStack montagefluo = null;
	        if(is_montage_vertical) {
	        	montagefluo = sc.combineVertically(montage.getStack(), scaledFluo.getStack());
	        } else {
	        	montagefluo = sc.combineHorizontally(montage.getStack(), scaledFluo.getStack());

	        }
	        
			montage = new ImagePlus(imageName+" Report", montagefluo);
		}
		
		montage.setTitle(imageName+" Report");
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

		d.addCheckbox("Crop_ROIs", true);

		d.addCheckbox("Consider_Z_Slices_Separately", false);
		d.addCheckbox("Set_Auto_Thresholds_On_Stack_Histogram", true);
		
		d.addMessage("Colocalization Result Options");
		
		d.addCheckbox("Get_Pearsons Correlation", true);
		d.addCheckbox("Get_Manders Coefficients", true);
		d.addCheckbox("Get_Overlap Coefficients", true);
		d.addCheckbox("Get_Li_ICA", true);
		d.addCheckbox("Get_Fluorogram", true);
		d.addCheckbox("Report_As_Vertical_Montage", false);
		
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
		
		doCropRois = d.getNextBoolean();
		doSeparateZ = d.getNextBoolean();
		
		is_stack_hist_z = d.getNextBoolean();
		doPearsons = d.getNextBoolean();
		doManders = d.getNextBoolean();
		doOverlap = d.getNextBoolean();
		doICA = d.getNextBoolean();
		doFluorogram = d.getNextBoolean();
		is_montage_vertical = d.getNextBoolean();
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
		
		IJ.run("JACoP B", "");
    }
}



