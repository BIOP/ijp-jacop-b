/*-
 * #%L
 * An update to the JACoP Plugin that helps in the management of ROIs, Z sections and helps generate cleaner reports
 * %%
 * Copyright (C) 2009 - 2024 Susanne Bolte, Fabrice P. Cordelières, Olivier Burri
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package ch.epfl.biop.coloc.utils;

import java.awt.Color;

import ij.gui.Overlay;
import ij.process.LUT;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.type.numeric.integer.IntType;
import sc.fiji.coloc.algorithms.AutoThresholdRegression;
import sc.fiji.coloc.algorithms.Histogram2D;
import sc.fiji.coloc.algorithms.MissingPreconditionException;
import sc.fiji.coloc.algorithms.PearsonsCorrelation;
import sc.fiji.coloc.algorithms.SpearmanRankCorrelation;

/*  JACoP: "Just Another Colocalization Plugin..." v1, 13/02/06
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
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.Blitter;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;
import sc.fiji.coloc.gadgets.DataContainer;
import sc.fiji.coloc.results.ResultHandler;
import sc.fiji.coloc.results.Warning;

/**
 *
 * @author Fabrice Cordelieres
 * @author Olivier Burri
 * @author Nicolas Chiaruttini
 */
public class ImageColocalizer {
    int width, height, nbSlices, depth, length, widthCostes, heightCostes, nbsliceCostes;
    String titleA, titleB;
    int[] A, B, M; // M for Mask
    int Amin, Amax, Bmin, Bmax;
    double Amean, Bmean;
    Calibration cal, micronCal;
    ResultsTable rt;
    
    //Values for stats
    boolean doThat;
    double sumA, sumB, sumAB, sumsqrA, Aarraymean, Barraymean;
    
    // OB: Add output for fluorogram ICQs and the likes to be separated from the actual processing (avoid .show())
    private Plot fluorogram_plot, cff_plot;
	private Plot ica_a_plot;
	private Plot ica_b_plot;
    
	// OB: ROI Handling
	private Roi roi=null;
    String roiName="";
	// OB: Get the images
	private ImagePlus impA, impB;
	
	// OB: Keep the Thresholds
	private int thrA, thrB;
	private ImageProcessor mask;
	    
    /** Creates a new instance of ImageColocalizer */
    public ImageColocalizer(ImagePlus ipA, ImagePlus ipB) {
        setup(ipA, ipB, new Calibration());
    }

    public ImageColocalizer(ImagePlus ipA, ImagePlus ipB, Calibration cal) {
        setup(ipA, ipB, cal);
    }
    
    /*
     * Convenience method to call the colocalizer using
     * a single multichannel image and the channels to be used for coloc (1-based)
     */
    public ImageColocalizer(ImagePlus imp, int channelA, int channelB) {
    	
    	
		Roi roi = imp.getRoi();
		imp.deleteRoi();
		ImagePlus[] channels = ChannelSplitter.split(imp);
		if(roi != null) {
			channels[channelA-1].setRoi(roi);
			channels[channelB-1].setRoi(roi);
		}
		
		
		setup(channels[channelA-1], channels[channelB-1], imp.getCalibration());
		
    }


    //double totalNumberPixInRoi = 0;
    
    /**
     * This method sets up the ROIs and the data for each image that is given.
     */
	public void setup(ImagePlus oipA, ImagePlus oipB, Calibration cal) {
    	// Add functionality to handle ROIs
    	// If either image has a ROI, clear the pixels outside the ROI
    	if (oipA.getRoi() != null || oipB.getRoi() != null) {
    		roi = (oipA.getRoi() != null ? oipA.getRoi() : oipB.getRoi());
    		
    		roiName = roi.getName();
    	
    		oipA.deleteRoi();
    		oipB.deleteRoi();
    		
    		impA = oipA.duplicate();
    		impB = oipB.duplicate();

    		oipA.setRoi(roi);
    		oipB.setRoi(roi);

    		ImageProcessor tip = oipA.getProcessor();
    		tip.setRoi(roi);
    		mask = Utils.getMask(impA, roi);

    		// Clear outside for both
    		for(int i = 1; i<=impA.getNSlices(); i++) {
    			Utils.clearOutside(impA.getStack().getProcessor(i), mask);
    			Utils.clearOutside(impB.getStack().getProcessor(i), mask);
    		}
    		// Set the title of the images
    		impA.setTitle(oipA.getTitle());
    		impB.setTitle(oipB.getTitle());
    	
    	} else { // If there are no ROIs, it's rather simpler

    		impA = oipA;
    		impB = oipB;

            mask = Utils.getMask(impA, roi);
    	}
    	
    	/*
    	 * Convenience initializations
    	 */
    	this.width=impA.getWidth();
        this.height=impA.getHeight();
        this.nbSlices=impA.getNSlices();
        this.depth=impA.getBitDepth();
        
        if (this.width!=impB.getWidth() || this.height!=impB.getHeight() || this.nbSlices!=impB.getNSlices() || this.depth!=impB.getBitDepth()){
            IJ.error("ImageColocalizer expects both images to have the same size and depth");
            return;
        }
        this.length=this.width*this.height*this.nbSlices;
        this.A=new int[this.length];
        this.B=new int[this.length];
        this.M=new int[this.length];
        this.titleA=impA.getTitle();
        this.titleB=impB.getTitle();
        this.cal=cal;
        this.micronCal=(Calibration) cal.clone();
        
        this.micronCal.pixelDepth/=1000;
        this.micronCal.pixelHeight/=1000;
        this.micronCal.pixelWidth/=1000;
        this.micronCal.setUnit("um");
        
        buildArray(impA, impB, mask);
        
        IJ.log("**************************************************\nImage A: "+this.titleA+"\nImage B: "+this.titleB);
        
        
        // OB: Build results Table
        rt = ResultsTable.getResultsTable();
        if(rt == null) { rt = new ResultsTable(); }
        
        // Increment the counter already
        rt.incrementCounter();
        rt.addValue("Image A", this.titleA);
        rt.addValue("Image B", this.titleB);
        if (roi != null) {
            rt.addValue("ROI", roiName);
        }
     }
   
    public void Pearson() {
        this.doThat=true;
        // OB: Add Pearsons Coefficient to Results Table
        double pearsons_coeff = linreg(A,B,0,0)[2];
        IJ.log("\nPearson's Coefficient:\nr="+round(pearsons_coeff,3));
        rt.addValue("Pearson's Coefficient", pearsons_coeff);
    }

    public void SpearmanRank() {
        this.doThat=true;

        SpearmanRankCorrelation<IntType> corr = new SpearmanRankCorrelation<>();

        int countPix = 0;

        // No threshold, takes ROI and mask into account

        for (int i=0; i<M.length; i++){
            if (M[i]>0) countPix++;
        }

        double[][] data = new double[countPix][2];
        countPix = 0;
        for (int i = 0; i<A.length; i++) {
            if (M[i]>0) {
                data[countPix][0] = A[i];
                countPix++;
            }
        }
        countPix = 0;
        for (int i = 0; i<B.length; i++) {
            if (M[i]>0) {
                data[countPix][1] = B[i];
                countPix++;
            }
        }

        double spearmanrank_coeff = corr.calculateSpearmanRank(data);

        IJ.log("\nSpearman's Rank Coefficient:\nr="+round(spearmanrank_coeff,3));
        rt.addValue("Spearman's Rank Coefficient", spearmanrank_coeff);
    }
    
    public void Overlap(){
    	double num=0;
        double numThr=0;
        double den1=0;
        double den1Thr=0;
        double den2=0;
        double den2Thr=0;
        
        for (int i=0; i<this.length; i++){
            num+=this.A[i]*this.B[i];
            den1+=Math.pow(this.A[i], 2);
            den2+=Math.pow(this.B[i], 2);
            if (this.A[i]>thrA && this.B[i]>thrB){
                numThr+=this.A[i]*this.B[i];
                den1Thr+=Math.pow(this.A[i], 2);
                den2Thr+=Math.pow(this.B[i], 2);
            }
        }
        
        double OverlapCoeff=num/(Math.sqrt(den1*den2));
        IJ.log("\nOverlap Coefficient:\nr="+round(OverlapCoeff,3));
        IJ.log("\nr^2=k1xk2:\nk1="+round(num/den1,3)+"\nk2="+round(num/den2,3));
        double OverlapCoeffThr=numThr/(Math.sqrt(den1Thr*den2Thr));
        IJ.log("\n\nUsing thresholds (thrA="+thrA+" and thrB="+thrB+")");
        IJ.log("\nOverlap Coefficient:\nr="+round(OverlapCoeffThr,3));
        IJ.log("\nr^2=k1xk2:\nk1="+round(numThr/den1Thr,3)+"\nk2="+round(numThr/den2Thr,3));
        
        // OB: Add Overlap Coefficients to Results Table
        rt.addValue("Overlap Coefficient",OverlapCoeff);
        rt.addValue("k1",num/den1);
        rt.addValue("k2",num/den2);

        rt.addValue("Threshold A", thrA);
        rt.addValue("Threshold B", thrB);
        rt.addValue("Thresholded Overlap Coefficient",OverlapCoeffThr);
        rt.addValue("Thresholded k1", numThr/den1Thr);
        rt.addValue("Thresholded k2", numThr/den2Thr);

    }

    /**
     * Area overlap measurements
     * Measures the area above the thresholds on both channels as well as the overlap
     */
    public void Areas() {

    	double AA=0, AB=0, AAB = 0;
    	double totalPix = 0;
    	
    	for (int i=0; i<this.length; i++){
    		if (this.A[i]>thrA) { AA++; }
    		if (this.B[i]>thrB) { AB++;	}
    		if (this.A[i]>thrA && this.B[i]>thrB) { AAB++;}
    		if (this.M[i]>0) {totalPix++;}
    	}
        
    	// Calibrate
        totalPix*=Math.pow(impA.getCalibration().pixelWidth,2);
    	AA*=Math.pow(impA.getCalibration().pixelWidth,2);
    	AB*=Math.pow(impA.getCalibration().pixelWidth,2);
    	AAB*=Math.pow(impA.getCalibration().pixelWidth,2);
    	
        IJ.log(
                "\nArea ROI ("+impA.getCalibration().getXUnit()+") Area tot = "+round(totalPix,0)+
                "\nArea Measurements ("+impA.getCalibration().getXUnit()+")\n Area A="+round(AA,0)+"Area B="+round(AB,0)+
        		"\n Area Overlap="+round(AAB,0));
        rt.addValue("Area tot", totalPix);
        rt.addValue("Area A", AA);
        rt.addValue("Area B", AB);
        rt.addValue("Area Overlap", AAB);
    }
    
    /**
     * Manders Coefficients
     * Calculation logic was untouched by OB
     */
    public void MM(){
    	double sumAcoloc=0;
        double sumAcolocThr=0;
        double sumA=0;
        double sumAThr=0;
        double sumBcoloc=0;
        double sumBcolocThr=0;
        double sumB=0;
        double sumBThr=0;
        
        for (int i=0; i<this.length; i++){
            if (this.B[i]>0){sumB+=this.B[i]; if(this.A[i]>0) sumAcoloc+=this.A[i]; }
            if (this.B[i]>thrB){sumBThr+=this.B[i]; if(this.A[i]>thrA) sumAcolocThr+=this.A[i];}
            if (this.A[i]>0) {sumA+=this.A[i]; if(this.B[i]>0) sumBcoloc+=this.B[i];}
            if (this.A[i]>thrA){sumAThr+=this.A[i]; if(this.B[i]>thrB) sumBcolocThr+=this.B[i]; }
        }
                
        double M1=sumAcoloc/sumA;
        double M1Thr=sumAcolocThr/sumAThr;
        double M2=sumBcoloc/sumB;
        double M2Thr=sumBcolocThr/sumBThr;

        IJ.log("\nManders' Coefficients (original):\nM1="+round(M1,3)+" (fraction of A overlapping B)\nM2="+round(M2,3)+" (fraction of B overlapping A)");
        IJ.log("\nManders' Coefficients (using threshold value of "+thrA+" for imgA and "+thrB+" for imgB):\nM1="+round(M1Thr,3)+" (fraction of A overlapping B)\nM2="+round(M2Thr,3)+" (fraction of B overlapping A)");
        
        // OB: Add Manders to Results Table
        rt.addValue("M1", M1);
        rt.addValue("M2", M2);
        
        rt.addValue("Threshold A", thrA);
        rt.addValue("Threshold B", thrB);
        
        rt.addValue("Thresholded M1", M1Thr);
        rt.addValue("Thresholded M2", M2Thr);
    
    }

    ResultHandler<IntType> resultHandler = new ResultHandler<IntType>() {
        @Override
        public void handleImage(RandomAccessibleInterval<IntType> image, String name) {

        }

        @Override
        public void handleHistogram(Histogram2D<IntType> histogram, String name) {

        }

        @Override
        public void handleWarning(Warning warning) {
            IJ.log(warning.getShortMessage());
        }

        @Override
        public void handleValue(String name, String value) {

        }

        @Override
        public void handleValue(String name, double value) {

        }

        @Override
        public void handleValue(String name, double value, int decimals) {

        }

        @Override
        public void process() {

        }
    };
    
    public void CostesAutoThr() {
        // See https://imagej.net/media/costes-etalcoloc.pdf
        // Coloc 2 implementation: https://github.com/fiji/Colocalisation_Analysis/blob/master/src/main/java/sc/fiji/coloc/algorithms/AutoThresholdRegression.java

        ArrayImg<IntType, IntArray> imgA = ArrayImgs.ints(A, A.length);
        ArrayImg<IntType, IntArray>  imgB = ArrayImgs.ints(B, B.length);
        ArrayImg<IntType, IntArray>  imgMask = ArrayImgs.ints(M, M.length);

        int CostesThrA, CostesThrB;

        try {
            DataContainer<IntType> dc = new DataContainer<>(imgA, imgB, 0,1,"image A", "image B", imgMask, new long[]{0}, new long[] {M.length});
            PearsonsCorrelation<IntType> pc = new PearsonsCorrelation<>();
            AutoThresholdRegression<IntType> ar = new AutoThresholdRegression<>(pc);
            ar.execute(dc);
            ar.processResults(resultHandler);
            CostesThrA = (int) ar.getCh1MaxThreshold().getRealDouble();
            CostesThrB = (int) ar.getCh2MaxThreshold().getRealDouble();
            IJ.log("\nCostes' automatic threshold set to "+CostesThrA+" for imgA & "+CostesThrB+" for imgB");

            rt.addValue("Threshold A", CostesThrA);
            rt.addValue("Threshold B", CostesThrB);

            // Seeing as it was set, set the thresholds here
            this.thrA = CostesThrA;
            this.thrB = CostesThrB;

            setThresholds(thrA, thrB);

        } catch (MissingPreconditionException e) {
            throw new RuntimeException(e);
        }

    }
    
    public void CCF(int CCFx){
        double meanA;
        double meanB;
        double nPoints;
        double num;
        double den1;
        double den2;
        double CCFmin=0;
        int lmin=-CCFx;
        double CCFmax=0;
        int lmax=-CCFx;

        double [] CCFarray=new double[2*CCFx+1];
        double [] x=new double[2*CCFx+1];
        
        int count=0;
        
        IJ.log("\nVan Steensel's Cross-correlation Coefficient between "+titleA+" and "+titleB+":");
        for (int l=-CCFx; l<=CCFx; l++){
            IJ.showStatus("CCF calculation in progress: "+(count+1)+"/"+(2*CCFx+1));
            IJ.showProgress(count+1, 2*CCFx+1);
            
            if (IJ.escapePressed()) {
                IJ.showStatus("Task canceled by user");
                IJ.showProgress(2,1);
                return;
            }

            meanA=0;
            meanB=0;
            nPoints=0;

            for (int k=1; k<=this.nbSlices; k++){
                for (int j=0; j<this.height; j++){
                    for (int i=0; i<this.width; i++){
                        if (i+l>=0 && i+l<this.width){
                            int coord=offset(i,j,k);
                            int coordShift=offset(i+l,j,k);

                            meanA+=this.A[coord];
                            meanB+=this.B[coordShift];
                            nPoints++;
                        }
                    }
                }
            }

            meanA/=nPoints;
            meanB/=nPoints;

            num=0;
            den1=0;
            den2=0;
            
            for (int k=1; k<=this.nbSlices; k++){
                for (int j=0; j<this.height; j++){
                    for (int i=0; i<this.width; i++){
                        if (i+l>=0 && i+l<this.width){
                            int coord=offset(i,j,k);
                            int coordShift=offset(i+l,j,k);
                            
                            num+=(this.A[coord]-meanA)*(this.B[coordShift]-meanB);
                            den1+=Math.pow((this.A[coord]-meanA), 2);
                            den2+=Math.pow((this.B[coordShift]-meanB), 2);
                        }
                    }
                }
            }

            double CCF=num/(Math.sqrt(den1*den2));
            
            if (l==-CCFx){
                CCFmin=CCF;
                CCFmax=CCF;
            }else{
                if (CCF<CCFmin){
                    CCFmin=CCF;
                    lmin=l;
                }
                if (CCF>CCFmax){
                    CCFmax=CCF;
                    lmax=l;
                }
            }
            x[count]=l;
            CCFarray[count]=CCF;
            count++;
        }
        IJ.log ("CCF min.: "+round(CCFmin,3)+" (obtained for dx="+lmin+") CCF max.: "+round(CCFmax,3)+" (obtained for dx="+lmax+")");

        cff_plot=new Plot("Van Steensel's CCF between "+this.titleA+" and "+this.titleB,"dx", "CCF",x,CCFarray);
        cff_plot.setLimits(-CCFx, CCFx, CCFmin-(CCFmax-CCFmin)*0.05, CCFmax+(CCFmax-CCFmin)*0.05);
        cff_plot.setColor(Color.white);
        cff_plot.draw();
        
        //Previous plot is white, just to get values inserted into the plot list, the problem being that the plot is as default a line plot... Following line plots same values as circles.
        cff_plot.setColor(Color.black);
        cff_plot.addPoints(x, CCFarray, Plot.CIRCLE);
        
        double[] xline={0,0};
        double[] yline={CCFmin-(CCFmax-CCFmin)*0.05,CCFmax+(CCFmax-CCFmin)*0.05};
        cff_plot.setColor(Color.red);
        cff_plot.addPoints(xline, yline, 2);
        
        CurveFitter cf=new CurveFitter(x, CCFarray);
        double[] param={CCFmin, CCFmax, lmax, (double) CCFx}; 
        cf.setInitialParameters(param);
        cf.doFit(CurveFitter.GAUSSIAN);
        param=cf.getParams();
        IJ.log("\nResults for fitting CCF on a Gaussian (CCF=a+(b-a)exp(-(xshift-c)^2/(2d^2))):"+cf.getResultString()+"\nFWHM="+Math.abs(round(2*Math.sqrt(2*Math.log(2))*param[3], 3))+" pixels");
        for (int i=0; i<x.length; i++) CCFarray[i]=CurveFitter.f(CurveFitter.GAUSSIAN, param, x[i]);
        cff_plot.setColor(Color.BLUE);
        cff_plot.addPoints(x, CCFarray, 2);
        
        IJ.showStatus("");
        IJ.showProgress(2,1);

        //plot.show();
        
        // OB: Add CFF to results table
        rt.addValue("CCF Min", CCFmin);
        rt.addValue("CCF Min DX Location", lmin);
        
        
        rt.addValue("CCF Max", CCFmax);
        rt.addValue("CCF Max DX Location", lmax);

        rt.addValue("CCF Fit FWHM", Math.abs(2*Math.sqrt(2*Math.log(2))*param[3]));
        

        
    }
    
    public void CytoFluo(){
        
        //Plot only accepts double array: convert A & B to Adb & Bdb    
        double[] Adb=int2double(this.A);    
        double[] Bdb=int2double(this.B);    
            
        fluorogram_plot = new Plot("Cytofluorogram between "+this.titleA+" and "+this.titleB, this.titleA, this.titleB, Adb, Bdb);
        double limHigh=Math.max(this.Amax, this.Bmax);
        double limLow=Math.min(this.Amin, this.Bmin);
        fluorogram_plot.setLimits(this.Amin, this.Amax, this.Bmin, this.Bmax);
        fluorogram_plot.setColor(Color.white);
        
        this.doThat=true;
        double[] tmp=linreg(this.A, this.B, 0, 0);
        double a=tmp[0];
        double b=tmp[1];
        double CoeffCorr=tmp[2];
        fluorogram_plot.draw();
        fluorogram_plot.setColor(Color.black);
        fluorogram_plot.addPoints(Adb, Bdb, 6);
        
        double[] xline={limLow,limHigh};
        double[] yline={a*limLow+b,a*limHigh+b};
        fluorogram_plot.setColor(Color.red);
        fluorogram_plot.addPoints(xline, yline, 2);

        IJ.log("\nCytofluorogram's parameters:\na: "+round(a,3)+"\nb: "+round(b,3)+"\nCorrelation coefficient: "+round(CoeffCorr,3));

    }

    public void ICA(){
        double[] Anorm=new double[this.length];
        double[] Bnorm=new double[this.length];
        double AnormMean=0;
        double BnormMean=0;
        double prodMin=0;
        double prodMax=0;
        double lim;
        double[] x= new double[this.length];
        double ICQ=0;

        //Intensities are normalized to range from 0 to 1
        for (int i=0; i<this.length; i++){
            Anorm[i]=(double) (this.A[i]-this.Amin)/this.Amax;
            Bnorm[i]=(double) (this.B[i]-this.Bmin)/this.Bmax;

            AnormMean+=Anorm[i];
            BnormMean+=Bnorm[i];
        }
        AnormMean=AnormMean/this.length;
        BnormMean=BnormMean/this.length;

        for (int i=0; i<this.length;i++){
            x[i]=(Anorm[i]-AnormMean)*(Bnorm[i]-BnormMean);
            if (x[i]>prodMax) prodMax=x[i];
            if (x[i]<prodMin) prodMin=x[i];
            if (x[i]>0) ICQ++;
       }

       lim = Math.max(Math.abs(prodMin), Math.abs(prodMax));

       ICQ=ICQ/this.length-0.5;

       ica_a_plot = new Plot("ICA A ("+this.titleA+")", "(Ai-a)(Bi-b)", this.titleA, new double[]{0, 0}, new double[]{0, 0});
       ica_a_plot.setColor(Color.white);
       ica_a_plot.setLimits(-lim, lim, 0, 1);
       ica_a_plot.draw();
       ica_a_plot.setColor(Color.black);
       ica_a_plot.addPoints(x, Anorm, Plot.DOT);
       ica_a_plot.draw();

       ica_a_plot.setColor(Color.red);
       ica_a_plot.drawLine(0, 0, 0, 1);

       ica_b_plot = new Plot("ICA B ("+this.titleB+")", "(Ai-a)(Bi-b)", titleB, new double[]{0, 0}, new double[]{0, 0});
       ica_b_plot.setColor(Color.white);
       ica_b_plot.setLimits(-lim, lim, 0, 1);
       ica_b_plot.draw();
       ica_b_plot.setColor(Color.black);
       ica_b_plot.addPoints(x, Bnorm, Plot.DOT);

       ica_b_plot.setColor(Color.red);
       ica_b_plot.drawLine(0, 0, 0, 1);

       IJ.log("\nLi's Intensity correlation coefficient:\nICQ: "+ICQ);

       // OB: Add CFF to results table
       rt.addValue("ICQ", ICQ);

    }

    public void RandomCostes2D(boolean binarize, int squareSize, int nShuffling, boolean showPlot, boolean showShuffledImage, boolean costesGraphBoundsUserSet, double xminCostesGraph, double xmaxCostesGraph) {
        randomCostes2D(impA,impB,squareSize,nShuffling,binarize, showPlot, showShuffledImage, costesGraphBoundsUserSet, xminCostesGraph, xmaxCostesGraph);
    }

    public ImagePlus randomCostesMaskPlot;
    public ImagePlus randomCostesPlot;

    public ImagePlus randomCostesMaskExampleShuffledImg;
    public ImagePlus randomCostesExampleShuffledImg;

    /**
     * Selects Blocks within ROI of size square
     * Performs Pearson correlation coefficient with these blocks, unshifted, and then by performing random shifts of the blocks (nShuffling)
     * Performing random shifts yields a distribution which is fitted to a gaussian
     * The observed pcc (with unshifted blocks) is then normalized to this gaussian, to provide the normalized Pcc
     * For instance pcc_normalized = -2 means the two channels are anti correlated, two sigmas away from randomness
     * For instance pcc_normalized = +2 means the two channels are correlated, two sigmas away from randomness
     * Warns if not enough shuffling possibilities are available
     * <br>
     * Threshold means the randomization is done with the thresholded images or not
     */
    public RandomCostes randomCostes2D(ImagePlus imgA, ImagePlus imgB, int squareSize, int nShuffling, boolean threshold, boolean showPlot, boolean showSampleImage,
                                       boolean costesGraphBoundsUserSet,
                                               double xminCostesGraph,
                                               double xmaxCostesGraph) {
        RandomCostes rc = new RandomCostes(imgA,imgB, squareSize, nShuffling,threshold, thrA, thrB);
        rc.compute();
        if (threshold) {
            randomCostesMaskExampleShuffledImg = rc.getExampleShuffleImage();
        } else {
            randomCostesExampleShuffledImg = rc.getExampleShuffleImage();
        }
        if (showPlot) {
            if (threshold) {
                randomCostesMaskPlot = rc.getPearsonDistributionGraph(costesGraphBoundsUserSet,
                        xminCostesGraph,
                        xmaxCostesGraph).getImagePlus();
            } else {
                randomCostesPlot = rc.getPearsonDistributionGraph(costesGraphBoundsUserSet,
                        xminCostesGraph,
                        xmaxCostesGraph).getImagePlus();
            }
        }
        if (threshold) {

            rt.addValue("Random Pearson Costes 2D (Mask)", rc.pearsonNormalized);

            rt.addValue("Random Pearson Costes 2D (Mask) pValueCorrelated", rc.pValueIsCorrelated);

            rt.addValue("Random Pearson Costes 2D (Mask) pValueAntiCorrelated", rc.pValueIsAntiCorrelated);
            IJ.log("\nNormalized random Pearson Costes 2D (Mask) nShuffle = "+nShuffling+" block Size = "+squareSize+" \nnpc="+round(rc.pearsonNormalized,3)+" \npValueCorrelated="+rc.pValueIsCorrelated+" \npValueAntiCorrelated="+rc.pValueIsAntiCorrelated);
        } else {
            rt.addValue("Random Pearson Costes 2D", rc.pearsonNormalized);

            rt.addValue("Random Pearson Costes 2D pValueCorrelated", rc.pValueIsCorrelated);

            rt.addValue("Random Pearson Costes 2D pValueAntiCorrelated", rc.pValueIsAntiCorrelated);

            IJ.log("\nNormalized random Pearson Costes 2D nShuffle = "+nShuffling+" block Size = "+squareSize+" \nnpc="+round(rc.pearsonNormalized,3)+" \npValueCorrelated="+rc.pValueIsCorrelated+" \npValueAntiCorrelated="+rc.pValueIsAntiCorrelated);
        }
        return rc;
    }
    
    //----------------------------------------------------------------------------------------------------------------------------------------------
    private void buildArray(ImagePlus imgA, ImagePlus imgB, ImageProcessor imgMask){
        int index=0;
        this.Amin=(int) Math.pow(2, this.depth);
        this.Amax=0;
        this.Amean=0;
        this.Bmin=this.Amin;
        this.Bmax=0;
        this.Bmean=0;
        
        for (int z=1; z<=this.nbSlices; z++){
            imgA.setSlice(z);
            imgB.setSlice(z);
            
            ImageStatistics stA=imgA.getStatistics();
            ImageStatistics stB=imgB.getStatistics();
            
            this.Amin=Math.min(this.Amin, (int) stA.min);
            this.Bmin=Math.min(this.Bmin, (int) stB.min);
            this.Amax=Math.max(this.Amax, (int) stA.max);
            this.Bmax=Math.max(this.Bmax, (int) stB.max);
                    
            this.Amean+=stA.pixelCount*stA.mean;
            this.Bmean+=stB.pixelCount*stB.mean;
            
            for (int y=0; y<this.height; y++){
                for (int x=0; x<this.width; x++){
                    this.A[index]=imgA.getProcessor().getPixel(x,y);
                    this.B[index]=imgB.getProcessor().getPixel(x,y);
                    this.M[index]=imgMask.getPixel(x,y);
                    index++;
                }
            }
            
            this.Amean/=this.length;
            this.Bmean/=this.length;
        }
    }

    public double[] linreg(int[] Aarray, int[] Barray, int TA, int TB){
         double num=0;
         double den1=0;
         double den2=0;
         double[] coeff=new double[6];
         int count=0;
         
         if (doThat){
             sumA=0;
             sumB=0;
             sumAB=0;
             sumsqrA=0;
             Aarraymean=0;
             Barraymean=0;
             for (int m=0; m<Aarray.length; m++){
                if ((Aarray[m]>=TA && Barray[m]>=TB)&&(M[m]>0)){
                    sumA+=Aarray[m];
                    sumB+=Barray[m];
                    sumAB+=Aarray[m]*Barray[m];
                    sumsqrA+=Math.pow(Aarray[m],2);
                    count++;
                }
             }
             Aarraymean=sumA/count;
             Barraymean=sumB/count;
         }
         
         for (int m=0; m<Aarray.length; m++){
            if ((Aarray[m]>=TA && Barray[m]>=TB)&&(M[m]>0)){
                num+=(Aarray[m]-Aarraymean)*(Barray[m]-Barraymean);
                den1+=Math.pow((Aarray[m]-Aarraymean), 2);
                den2+=Math.pow((Barray[m]-Barraymean), 2);
            }
         }
        
        //0:a, 1:b, 2:corr coeff, 3: num, 4: den1, 5: den2
        coeff[0]=(count*sumAB-sumA*sumB)/(count*sumsqrA-Math.pow(sumA,2));
        coeff[1]=(sumsqrA*sumB-sumA*sumAB)/(count*sumsqrA-Math.pow(sumA,2));
        coeff[2]=num/(Math.sqrt(den1*den2));
        coeff[3]=num;
        coeff[4]=den1;
        coeff[5]=den2;
        return coeff;
     }

    private double[] int2double(int[] input){
        double[] output=new double[input.length];
        for (int i=0; i<input.length; i++) output[i]=input[i];
        return output;
    }
    
    /** Returns the index where to find the information corresponding to pixel (x, y, z).
     * @param x coordinate of the pixel.
     * @param y coordinate of the pixel.
     * @param z coordinate of the pixel.
     * @return the index where to find the information corresponding to pixel (x, y, z).
     */
    private int offset(int x, int y, int z){
        if (x+y*this.width+(z-1)*this.width*this.height>=this.width*this.height*this.nbSlices){
            return this.width*this.height*this.nbSlices-1;
        }else{
            if (x+y*this.width+(z-1)*this.width*this.height<0){
                return 0;
            }else{
                return x+y*this.width+(z-1)*this.width*this.height;
            }
        }
    }
    
    public int offsetCostes(int m,int n,int o){
        if (m+n*this.widthCostes+(o-1)*this.widthCostes*this.heightCostes>=this.widthCostes*this.heightCostes*this.nbsliceCostes){
            return this.widthCostes*this.heightCostes*this.nbsliceCostes-1;
        }else{
            if (m+n*this.widthCostes+(o-1)*this.widthCostes*this.heightCostes<0){
                return 0;
            }else{
                return m+n*this.widthCostes+(o-1)*this.widthCostes*this.heightCostes;
            }
        }
    }
    
    public double round(double y, int z){
         //Special tip to round numbers to 10^-2
         y*=Math.pow(10,z);
         y=(int) y;
         y/=Math.pow(10,z);
         return y;
    }
    
	public ImagePlus getFluorogramImage() {
		float[] minmax = getFluorogramMinMax();
		return getFluorogramImage(256, (int)Math.floor(minmax[0]), (int)Math.floor(minmax[1]));	
	}
	
	public float[] getFluorogramMinMax() {
		Plot fp = getFluorogram();
		float[] valA = fp.getXValues();
		float[] valB = fp.getYValues();
		
		
		float minA=Float.MAX_VALUE;
		float maxA=0;
		
		float minB=Float.MAX_VALUE;
		float maxB=0;
		
		for(int i=0; i<valA.length; i++) {
			
			if(minA > valA[i]) {
				minA=valA[i];
			}
			if(minB > valB[i]) {
				minB=valB[i];
			}
			if(maxA < valA[i]) {
				maxA=valA[i];
			}
			if(maxB < valB[i]) {
				maxB=valB[i];
			}
		}
		float min, max;
		min = Math.min(minA, minB);
		max = Math.max(maxA, maxB);

        return new float[]{min, max};
	}

    /*
     * Builds a pretty fluorogram for display
     */
	public ImagePlus getFluorogramImage(int nbins, int min, int max) {

		Plot fp = getFluorogram();

		float[] valA = fp.getXValues();
		float[] valB = fp.getYValues();
		
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
		
		// Build the grayLevels between the min and the max
		int binWidth= (int) Math.round((double)(max-min) / (double)(nbins-1));
		for(int i=0; i<nbins; i++) {
			for(int j=0; j<lut_size; j++) {
				int val =  (i+1)*binWidth;
				grad_A.set(i,j,val);
				grad_B.set(j,nbins-i-1,val);
			}
		}
		
		grad_A.setLut(impA.getLuts()[0]);
		grad_B.setLut(impB.getLuts()[0]);
		// Somehow load the Fire LUT
		fluo_ip.setLut(fireLUT());
		fluo_ip.setMinAndMax(0, 6);
		
		fluo_ip.convertToRGB();
		
		// Add threshold lines
		// Normalize Threshold to min/max
		int thrAb = (int) Math.round( (((double) thrA - (double) min) / (double) max * (double) (nbins-1)));
		int thrBb = (int) Math.round( (((double) thrB - (double) min) / (double) max * (double) (nbins-1)));
		
		fluo_ip.setColor(new Color(255,255,255));
		fluo_ip.drawLine(thrAb, nbins, thrAb, 0);

		//fluo_ip.setColor(impB.getLuts()[0].getColorModel().);
		fluo_ip.drawLine(0, nbins-thrBb-1, nbins, nbins-thrBb-1);
		
		// Build fluorogram
		fluo.copyBits(grad_A.convertToRGB(), lut_size+1, nbins+1, Blitter.ADD);
		fluo.copyBits(grad_B.convertToRGB(), 0, 0, Blitter.ADD);
		fluo.copyBits(fluo_ip, lut_size+1, 0, Blitter.ADD);
		
		
		return new ImagePlus(fp.getTitle(), fluo);
	}
    
    // OB: Extra functions to return the results
    public void showResults() {
    	rt.show("Results");
    }

	public Plot getCFFPlot() {
		return cff_plot;
	}

	public Plot getICAaPlot() {
		return ica_a_plot;
	}

	public Plot getICAbPlot() {
		return ica_b_plot;
	}
	
	public Plot getFluorogram() {
		return fluorogram_plot;
	}
	
	public void setThresholds(int thrA, int thrB) {
		
		this.thrA = thrA;
		this.thrB = thrB;

        if(roi!=null) {
            for(int i=0; i<impA.getStackSize(); i++) {
                impA.getStack().getProcessor(i+1).setMask(mask);
                impB.getStack().getProcessor(i+1).setMask(mask);
            }
            impA.setRoi(roi);
            impB.setRoi(roi);
        }

		impA.getProcessor().setThreshold(this.thrA , impA.getProcessor().getMaxThreshold(), ImageProcessor.NO_LUT_UPDATE);
		impB.getProcessor().setThreshold(this.thrB , impB.getProcessor().getMaxThreshold(), ImageProcessor.NO_LUT_UPDATE);
	}
	
	public int getThresholdA() {
		return this.thrA;
	}
	
	public int getThresholdB() {
		return this.thrB;
	}
	
	
	public void setThresholds(String thrMetA, String thrMetB) {
		
		// ON 16 BIT IMAGES, the DISPLAY RANGE IS RESET!!!!!!!!!
		double[][] rangeA = Utils.getDisplayRange(impA);
		double[][] rangeB = Utils.getDisplayRange(impB);
  		
		if(roi!=null) {
			for(int i=0; i<impA.getStackSize(); i++) {
				impA.getStack().getProcessor(i+1).setMask(mask);
				impB.getStack().getProcessor(i+1).setMask(mask);
			}
		}
		
		// Allow for thresholds to be either manual or automatic
		if(thrMetA.matches("\\d*")) {
			this.thrA = Integer.parseInt(thrMetA);
			impA.getProcessor().setThreshold(this.thrA , impA.getProcessor().getMaxThreshold(), ImageProcessor.RED_LUT);
            impA.setRoi(roi);

		} else {
		
			// Add the ROIs, to ensure proper thresholds
			
			impA.getProcessor().resetThreshold();
			impA.setRoi(roi);
			
			// Bug with setAutoThreshold, check whether to use stack or not by hand...
			// Otherwise we get an incorrect value
			if(impA.getNSlices() > 1) {
				IJ.setAutoThreshold(impA, thrMetA+" dark stack");
			} else {
				IJ.setAutoThreshold(impA, thrMetA+" dark");
			}
			this.thrA = (int) impA.getProcessor().getMinThreshold();
		}
		
		// Allow for thresholds to be either manual or automatic
		if(thrMetB.matches("\\d*")) {
			//IJ.log("Threshold "+thrMetB+" matches regex");
			this.thrB = Integer.parseInt(thrMetB);
			impB.getProcessor().setThreshold(this.thrB , impB.getProcessor().getMaxThreshold(), ImageProcessor.RED_LUT);
            impB.setRoi(roi);
			
		} else {
		
			impB.getProcessor().resetThreshold();
			impB.setRoi(roi);
			
			if(impB.getNSlices() > 1) {
				IJ.setAutoThreshold(impB, thrMetB+" dark stack");
			} else {
				IJ.setAutoThreshold(impB, thrMetB+" dark");
			}
			
			this.thrB = (int) impB.getProcessor().getMinThreshold();
		}
		
		rt.addValue("Auto Threshold A", thrMetA);
		rt.addValue("Auto Threshold B", thrMetB);
		
		Utils.setDisplayRange(impA, rangeA);
		Utils.setDisplayRange(impB, rangeB);

				
	}
	
	public ImagePlus getImageA() {
		return impA;
	}

	public ImagePlus getImageB() {
		return impB;
	}
	
	public ImagePlus getRGBImage(ImagePlus imp, Boolean is_zProject, float strokeWidth) {
		
		imp.killRoi();
		ImagePlus impr = imp.duplicate();
		impr.setTitle(imp.getTitle());
			
		if(is_zProject) {
			ZProjector zp = new ZProjector(impr);
			zp.setMethod(ZProjector.MAX_METHOD);
			zp.doHyperStackProjection(true);
			impr = zp.getProjection();
		}
		
		return flattenRoi(impr, strokeWidth);
	}

    public ImagePlus getRandomCostesMaskPlot(Boolean is_Z, float strokeWidth) {
        return getRGBImage(impA, is_Z, strokeWidth);
    }

	public ImagePlus getRGBImageA(Boolean is_Z, float strokeWidth) {
		return getRGBImage(impA, is_Z, strokeWidth);
	}
	
	public ImagePlus getRGBImageB(Boolean is_Z, float strokeWidth) {
		return getRGBImage(impB, is_Z, strokeWidth);
	}
	
	public ImagePlus getRGBColocImage(float strokeWidth) {
		// This should return a rgb image of the composite of the two channels, with the pixels as a mask
		// Make a composite of the two images
		ImagePlus[] images = {impA, impB};
		
		ImagePlus comp = RGBStackMerge.mergeChannels(images, true);
		return flattenRoi(comp, strokeWidth);
		
	}
	
	public ImagePlus getMaskA() {
		return Utils.binarize(impA, thrA);
		
	}
	
	public ImagePlus getMaskB() {
		return Utils.binarize(impB, thrB);
	}
	
	public ImagePlus getRGBMaskA(float strokeWidth) {
		return flattenRoi(Utils.binarize(impA, thrA), strokeWidth);
		
	}
	
	public ImagePlus getRGBMaskB(float strokeWidth) {
		return flattenRoi(Utils.binarize(impB, thrB), strokeWidth);
		
	}
	
	public ImagePlus getANDMask() {
		ImagePlus impMA = getMaskA();
		ImagePlus impMB = getMaskB();
		
		ImagePlus impAND = new ImagePlus(impA.getTitle()+" AND "+impB.getTitle(), impMA.getStack());
		// DO an AND
		int nSlices = impMA.getStackSize();
		if(nSlices > 1) {
			ImageStack stackA = impMA.getStack();
			ImageStack stackB = impMB.getStack();
			ImageStack andStack = impMA.createEmptyStack();
			ImageProcessor ip;
			for(int i=1; i<=nSlices; i++) {
				ip = stackA.getProcessor(i).duplicate();
				ip.copyBits(stackB.getProcessor(i), 0, 0, Blitter.AND);
				andStack.addSlice(ip);	
			}
			impAND.setStack(andStack);
			return impAND;
		} else {
			ImageProcessor ip = impMA.getProcessor().duplicate();
			ip.copyBits(impMB.getProcessor(), 0, 0, Blitter.AND);
			impAND.setProcessor(ip);
			return impAND;
		}
	}
	
	public ImagePlus getRGBANDMask(float strokeWidth) {
		return flattenRoi(getANDMask(), strokeWidth);
	}
	
	
	
	// Proper flattening of the ROI onto the image or stack
	private ImagePlus flattenRoi(ImagePlus imp, float strokeWidth) {
		if(roi != null) {
            if (strokeWidth>0) {
                roi.setStrokeColor(Color.white);
                roi.setStrokeWidth(strokeWidth);
                imp.setOverlay(new Overlay(roi));
                imp.setHideOverlay(false);
            }
		} else {
			// Make an empty overlay
			imp.setOverlay(new Roi(1, 1, 1, 1), new Color(0,0,0,0), 1, new Color(0,0,0,0));
			imp.setHideOverlay(false);
			
		}
		if (imp.getStackSize() >1) {
			imp.flattenStack();
		} else {
			imp = imp.flatten();
		}
		return imp;
		
	}
	

	
	public void addResult(String name, int value) {
		rt.addValue(name, value);
	}
    public void addResult(String name, double value) {
        rt.addValue(name, value);
    }
	
	public void addResult(String name, String value) {
		rt.addValue(name, value);
	}
	
	public void removeLastRow() {
		rt.deleteRow(rt.getCounter()-1);
	}

    // Fire LUT

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
