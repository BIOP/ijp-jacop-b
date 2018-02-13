/*JACoP: "Just Another Colocalization Plugin..." v1, 13/02/06
    Fabrice P Cordelieres, fabrice.cordelieres at curie.u-psud.fr
    Susanne Bolte, Susanne.bolte@isv.cnrs-gif.fr
 
    Copyright (C) 2006 Susanne Bolte & Fabrice P. Cordelieres
    
    Readaptation of code by Olivier Burri v1, 20.09.2017
    
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
import ch.epfl.biop.coloc.utils.Utils;
import ch.epfl.biop.montage.StackMontage;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.MontageMaker;
import ij.plugin.PlugIn;
import ij.plugin.StackCombiner;
import ij.plugin.Thresholder;
import ij.plugin.frame.RoiManager;

public class JACoP_B implements PlugIn {
    
	int channelA, channelB;
	
	ImagePlus impA, impB;
	
	String thrA, thrB;
	
	int mThrA, mThrB;
	
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


	private boolean use_advanced;

	private boolean is_auto_fluo = true;
	private int fluo_bins = 256;
	private int fluo_min=0;
	private int fluo_max=255;
	
	@Override
	public void run(String arg) {
		
		if (!mainDialog()) return;
	
		// Make some sense of it all
		
		//Switch from folder to single image mode			
		// Setup the image and all settings
		this.imp = IJ.getImage();
		runColoc();
			
	}
	

	/*
	 * As it currently stands, the ImageColocalizer will do a coloc analysis of
	 * Time, Z, all bundled together, so we work to split them BEFORE
	 * giving them to the ImageColocalizer
	 * STEPS: 
	 *  A - Split timepoints, always
	 *  B - Split Z slices, as needed.
	 *  
	 *  In both cases, care should be taken to manage ROIs properly
	 *  STEPS: 
	 *  A - If there is a ROI and we want to crop the images, duplicate the ROI and offset it
	 *  B - If there is no ROI, just hand the image over
	 *  C - If there is a ROI but no crop, duplicate the image without the ROI and add it again
	 *  
	 *  These handle whether the ROI is on the RoiManager or on the image
	 *  NOTE that the ROIManager will always win. So if there is a ROI on the image that is not 
	 *  in the ROI Manager too, it will be ignored (and lost)
	 */
	private void runColoc() {
		
		// Check which thresholds to use
		if(thrA.matches("Use Manual Threshold Below")) thrA = String.valueOf(mThrA);
		if(thrB.matches("Use Manual Threshold Below")) thrB = String.valueOf(mThrB);
		
		if(thrA.matches("Costes Auto-Threshold") || thrA.matches("Costes Auto-Threshold"))  doCostesThr = true;

		// He who makes the magic happen
		ImageColocalizer ic;
		
		// Store results, as there might be many runs of ImageColocalizer
		List<ImagePlus> results = new ArrayList<ImagePlus>();
		
		// 	The ROI we need to haul around in case it's there
		if(imp.getRoi() != null)
			roi = (Roi) imp.getRoi().clone();
		
		// Get eventual ROIs from the Manager
		RoiManager rm = RoiManager.getInstance();
		if(rm == null) rm = new RoiManager(false);
		// If this flag is set, then we will look for ROIs in the Roi Manager
		if (rm.getCount() > 0) hasRoiSets = true;
		
		// Because we will be looping ROIs no matter what
		// we need to have a ROI counter equal to at least 1
		// Otherwise it will exit the loop without having done anythin
		int rcount = (rm.getCount()>0) ? rm.getCount() : 1;
		
		// LOOP 1: All ROIs, or single ROI if hasRoiSets is false
		for(int r=0; r<rcount; r++) {
			// Choose the ROI here
			if(hasRoiSets) roi = (Roi) rm.getRoi(r).clone(); 
			
			// Get some data on the ROI, is there one? Does it have a name?
			String roiName = null;
			if (roi != null) {
				roiName = roi.getName();
				if (roiName == null) {
	    			roiName = "ROI";
	    		}
				
				// Use to name final results wit the ROI name
				imageName = imp.getTitle()+" ("+roiName+")";
			} else {
				imageName = imp.getTitle();
			}
						
			// Set the ROI, could be null but that's not a problem
			imp.setRoi(roi);
			
			// Single timepoint temporary ImagePlus
			ImagePlus singleTimeImp;
			
			// Loop through timepoints
			int nT = imp.getNFrames();
		
			// LOOP 2: Time
			for(int t=1; t<=nT; t++) {
				
				singleTimeImp = Utils.cropTime(imp, roi, t, doCropRois);
				
				if(roi != null) roi = (Roi) singleTimeImp.getRoi().clone();
				
				// Do Z separately, maybe
				ImagePlus[] zImages = Utils.cropSlices(singleTimeImp, roi, doSeparateZ);
				
				int thrAH = 0;
				int thrBH = 0;
				
				
				// Little Hacky: Make sure to get the right thresholds if we are working on the stack histogram...
				// We need to initialize ImageColocalizer, get the thresholds and save them
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
					IJ.log("Stack Histogram Threshold B Set: "+thrB);
					IJ.log("Calculated Stack Histogram Threshold B: "+thrBH);
					
					// Because initializing the colocalizer writes a result already,
					// We need to remove this ghost value
					ic.removeLastRow();
				} // END CONDITION IS STACK Z for getting threhsolds
				
				
				// The resulting image will be RGB and potentially a stack of T and Z
				ImageStack res = null;
				
				// This final loop will do the actual work
				
				// Temp placeholder for the resulting image that will need to be appended
				ImagePlus tmp = null;
								
				// LOOP 3: Slices
				for(int i=0; i< zImages.length;i++) {
					// re-add the ROI to the dataset
					if(roi != null) {
						// OHMYGODDONTASKMEABOUTHIS
						roi.setPosition(0);  
						// ... no what this does is ensure the ROI has no position data
						// otherwise it does not get copied properly...
						
						// Normally there would only be the line below
						zImages[i].setRoi(roi);
					}
					
					
					// Run some magic
					ic = new ImageColocalizer( zImages[i], channelA, channelB);
					
					// Set the threshold to use, depending on the user choice
					if(doCostesThr) {
						ic.CostesAutoThr();
					} else {
						ic.setThresholds(thrA, thrB);
					}
					
					// Set the thresholds again in case they are to be made on the stack histogram.
					if( is_stack_hist_z) {
						ic.setThresholds(thrAH, thrBH);

						// Let the user know we used the stack histogram
						ic.addResult("Using Stack Histogram", "True");
					} else { 
						ic.addResult("Using Stack Histogram", "False");
					}
					
					// Do the analysis, finally...
					tmp = runAnalysis(ic);
					
					// Create the empty stack from the result of tmp...
					// Looks a bit ugly but we are unsure of the size of the image
					// Before the first run
					if( i==0 ) {
						res = tmp.createEmptyStack();
					}
					
					// In the case we do Z, we need to catch a single imageProcessor and add it as a slice
					// Otherwise we take the stack directly0
					if(zImages.length > 1) { 
						res.addSlice(tmp.getProcessor());
					} else {
						res = tmp.getStack();
					}
					
					// Some extra data results to add, for the user to locate his data slice in the results
					ic.addResult("Timepoint", t);
					if(zImages.length > 1) ic.addResult("Slice", i+1);
					
					ic.showResults();
					
				} // END LOOP 3: Slices
				
				// Reuse tmp to make the final ImagePlus for all Z
				tmp = new ImagePlus(tmp.getTitle(), res); 
				
				// This list contains all the results for all ROIs and timepoints
				results.add(tmp);
			} // END LOOP 2: Time
			
			// Perhaps here we could make hyperstack in time as well, but it does not seem too useful...
			
		} // END LOOP 1: ROIs
		
		
		// Show all images
	 	for(ImagePlus i : results) {
	 		i.show();
	 	}
	}
	
	

	/*
	 * This method actually runs the analysis and outputs a pretty report
	 * At least... Oli thinks it's pretty
	 */
	private ImagePlus runAnalysis(ImageColocalizer ic) {
		if(doPearsons) ic.Pearson();
		if(doManders) ic.MM();
		if(doOverlap) ic.Overlap();
		if(doICA) ic.ICA();
		if(doFluorogram) ic.CytoFluo();
		// Add Areas
		ic.Areas();
		
		// Define rows and columns here in case we want a vertical report
		int rows = 3;
		int columns = 2;
		
		
		// Get the images and make a montage
		ArrayList<ImagePlus> imgs = new ArrayList<ImagePlus>();
        
		// Vertical montage
		if(is_montage_vertical) {
			imgs.add(ic.getRGBImageA(false));
			imgs.add(ic.getRGBMaskA());
			
			imgs.add(ic.getRGBImageB(false));
			imgs.add(ic.getRGBMaskB());
			
			imgs.add(ic.getRGBColocImage());
			imgs.add(ic.getRGBANDMask());
			
        // Horizontal Montage
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
        
        // Need to do a different thing if we are using stack or a single slice.        
        if(imgs.get(0).getNSlices() == 1) {
        	ImageStack montagestk = imgs.get(0).createEmptyStack();
            // Make a normal montage
        	 	for(ImagePlus i : imgs) {
        	 		montagestk.addSlice(i.getProcessor());
        	 	}
        	 	// Use Montage Maker
        	 	MontageMaker mm = new MontageMaker();
            	montage = mm.makeMontage2(new ImagePlus("for montage",montagestk), columns, rows, 1.0, 1, imgs.size(), 1, 0, false);

        } else {
        	// We can use Oli's Stack Montage for convenience
        	montage = StackMontage.montageImages(imgs,rows, columns);
        }
    	
        //Eventually add the fluorogram
		if(doFluorogram) {
			ImagePlus fluo = (is_auto_fluo) ? ic.getFluorogramImage() : ic.getFluorogramImage(fluo_bins, fluo_min,fluo_max);
			
			ImagePlus scaledFluo = null;
			
			// Scale the fluorogram to the width or height of the image
	        if(is_montage_vertical) {
	        	scaledFluo = Utils.scale(fluo, montage.getWidth());
	        } else {
	        	scaledFluo = Utils.scale(fluo, montage.getHeight());

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
		
		// Make sure it's the title we want.
		montage.setTitle(imageName+" Report");
		
		// And VIOLA
		return montage;
	}

	private Boolean mainDialog() {
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
		d.addCheckbox("Set Advanced Parameters", false);

		
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
		
		use_advanced = d.getNextBoolean();
		
		if(use_advanced) {
			return advancedDialog();
		}
		
		return true;
    }
	
	
    
    
    private Boolean advancedDialog() {
		// Advanced features, like Fluorogram bins
		GenericDialogPlus d = new GenericDialogPlus("Advanced Parameters");
    	d.addCheckbox("Auto-Adjust Fluorogram Per Image", true);
    	d.addMessage("Otherwise, use parameters below");
    	d.addNumericField("Fluorogram_Bins", 256, 0);
		d.addNumericField("Fluorogram_Min", 0, 0);
		d.addNumericField("Fluorogram_Max", 255, 0);
		
		d.showDialog();
		if(d.wasCanceled()) {
			return false;
		}
		
		is_auto_fluo = d.getNextBoolean();
		fluo_bins = (int) d.getNextNumber();
		fluo_min  = (int) d.getNextNumber();
		fluo_max  = (int) d.getNextNumber();

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
	//	ImagePlus imp = IJ.openImage("F:\\People\\Nadine Schmidt\\20170712_KN35_IF_A1_2_LUT_BC.tif");
		
		ImagePlus imp = IJ.openImage("E:\\JACOP\\20170816_KN47_IF_B1_2_LUT_BC.tif");
		RoiManager rm =  new RoiManager();
		rm.runCommand("Open", "E:\\JACOP\\singleCell.roi");
		imp.show();

		IJ.run("JACoP B", "");
    }
}