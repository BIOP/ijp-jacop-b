package ch.epfl.biop.coloc.utils;

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

import ij.*;
import ij.gui.*;
import ij.ImagePlus.*;
import ij.measure.*;
import ij.plugin.ChannelSplitter;
import ij.plugin.RGBStackMerge;
import ij.plugin.Thresholder;
import ij.plugin.ZProjector;
import ij.process.*;

import java.awt.*;

import fiji.threshold.Auto_Threshold;

/**
 *
 * @author Fabrice Cordelieres
 */
public class ImageColocalizer {
    int width, height, nbSlices, depth, length, widthCostes, heightCostes, nbsliceCostes, lengthCostes;
    String titleA, titleB;
    int[] A, B;
    int Amin, Amax, Bmin, Bmax;
    double Amean, Bmean;
    Calibration cal, micronCal;
    Counter3D countA, countB;
    ResultsTable rt; 
    
    
    //Values for stats
    boolean doThat;
    double sumA, sumB, sumAB, sumsqrA, Aarraymean, Barraymean;
    
    // OB: Add output for fluorogram ICQs and the likes to be separated from the actuall processing (avoid .show())
    private Plot fluorogram_plot, icq_plot, costes_plot, cff_plot;
    
    // OB: Add output for Costes Mask
    private ImagePlus costesMask;
	private Plot ica_a_plot;
	private Plot ica_b_plot;
	private ImagePlus costes_rand;
	private Plot costes_rand_plot;
	private ImagePlus centers_img;
	private ImagePlus coinc_img;
    
	// OB: ROI Handling
	private Roi roi=null;
    private Color roiColor = Color.white;
    String roiName="";
	// OB: Get the images
	private ImagePlus impA, impB;
	
	// OB: Keep the Thresholds
	private int thrA, thrB;
	
    /** Creates a new instance of ImageColocalizer */
    
    /** Creates a new instance of ImageColocalizer */
    public ImageColocalizer(ImagePlus ipA, ImagePlus ipB) {
        setup(ipA, ipB, new Calibration());
    }
    
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
    
    
    
	public void setup(ImagePlus oipA, ImagePlus oipB, Calibration cal) {
    	// Add functionality to handle ROIs
    	// If either image has a ROI, start by clearing the pixels outside of the ROI, and use the ROI name in the image names
    	if (oipA.getRoi() != null || oipB.getRoi() != null) {
    		roi = (oipA.getRoi() != null ? oipA.getRoi() : oipB.getRoi());
    		
    		roiName = roi.getName();
    		if (roiName == null) {
    			roiName = "ROI";
    		}
    		oipA.deleteRoi();
    		oipB.deleteRoi();
    		
    		impA = oipA.duplicate();
    		impB = oipB.duplicate();

    		
    		oipA.setRoi(roi);
    		oipB.setRoi(roi);
    		
    		// Clear outside for both
    		for(int i = 1; i<=impA.getNSlices(); i++) {
    			impA.getStack().getProcessor(i).fillOutside(roi);
    			impB.getStack().getProcessor(i).fillOutside(roi);
    		}
    		// Set the title of the images
    		impA.setTitle(oipA.getTitle());
    		impB.setTitle(oipB.getTitle());
    		
    	} else {
    		impA = oipA;
    		impB = oipB;
    	}
    	
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
        this.titleA=impA.getTitle();
        this.titleB=impB.getTitle();
        this.cal=cal;
        this.micronCal=(Calibration) cal.clone();
        
        this.micronCal.pixelDepth/=1000;
        this.micronCal.pixelHeight/=1000;
        this.micronCal.pixelWidth/=1000;
        this.micronCal.setUnit("um");
        
        buildArray(impA, impB);
        
        IJ.log("**************************************************\nImage A: "+this.titleA+"\nImage B: "+this.titleB);
        
        // OB: Build results Table
        rt = ResultsTable.getResultsTable();
        if(rt == null) { rt = new ResultsTable(); }
        
        // Increment the counter already
        rt.incrementCounter();
        rt.addValue("Image A", this.titleA);
        rt.addValue("Image B", this.titleB);
        if (roi != null)
        rt.addValue("ROI", roiName);
        
     }
    
   
    public void Pearson() {
        this.doThat=true;
        // OB: Add Pearsons Coefficient to Results Table
        double pearsons_coeff = linreg(A,B,0,0)[2];
        IJ.log("\nPearson's Coefficient:\nr="+round(pearsons_coeff,3));
        rt.addValue("Pearson's Coefficient", pearsons_coeff);
    }
    
    /*public void Pearson(int TA, int TB) {
        doThat=true;
        IJ.log("\nPearson's Coefficient using thesholds:\nr="+round(linreg(A,B,TA,TB)[2],3));
    }*/
    
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
    /*
     * Area overlap measurements
     * Measures the area above the thresholds on both channels as well as the overlap
     */
    public void Areas() {
    	double AA=0, AB=0, AAB = 0;
    	
    	for (int i=0; i<this.length; i++){
    		if (this.A[i]>thrA) { AA++; }
    		if (this.B[i]>thrB) { AB++;	}
    		if (this.A[i]>thrA && this.B[i]>thrB) { AAB++;}
    	}
        
    	// Calibrate
    	AA*=Math.pow(impA.getCalibration().pixelWidth,2);
    	AB*=Math.pow(impA.getCalibration().pixelWidth,2);
    	AAB*=Math.pow(impA.getCalibration().pixelWidth,2);
    	
        IJ.log("\nArea Measurements ("+impA.getCalibration().getXUnit()+")\n Area A="+round(AA,0)+"Area B="+round(AB,0)+
        		"\n Area Overlap="+round(AAB,0));

    	rt.addValue("Area A", AA);
        rt.addValue("Area B", AB);
        rt.addValue("Area Overlap", AAB);
    }
    
    /*
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
    
    public void CostesAutoThr() {
        int CostesThrA=this.Amax;
        int CostesThrB=this.Bmax;
        double CostesSumAThr=0;
        double CostesSumA=0;
        double CostesSumBThr=0;
        double CostesSumB=0;
        double CostesPearson=1;
        double [] rx= new double[this.Amax-this.Amin+1];
        double [] ry= new double[this.Amax-this.Amin+1];
        double rmax=0;
        double rmin=1;
        this.doThat=true;
        int count=0;
        
        //First Step: define line equation
        this.doThat=true;
        double[] tmp=linreg(this.A, this.B, 0, 0);
        double a=tmp[0];
        double b=tmp[1];
        double CoeffCorr=tmp[2];
        this.doThat=false;
        
        int LoopMin= (int) Math.max(this.Amin, (this.Bmin-b)/a);
        int LoopMax= (int) Math.min(this.Amax, (this.Bmax-b)/a);
        
        
        //Minimize r of points below (thrA,a*thrA+b)
        for (int i=LoopMax;i>=LoopMin;i--){
            IJ.showStatus("Costes' threshold calculation in progress : "+(int)(100*(LoopMax-i)/(LoopMax-LoopMin))+"% done");
            IJ.showProgress(LoopMax-i, LoopMax-LoopMin);
            
            if (IJ.escapePressed()) {
                IJ.showStatus("Task canceled by user");
                IJ.showProgress(2,1);
                return;
            }
            
            CostesPearson=linregCostes(this.A, this.B, i, (int) (a*i+b))[2];
            
            rx[count]=i;
            ry[count]=CostesPearson;
            if (((Double) CostesPearson).isNaN()){
                if (count!=LoopMax){
                    ry[count]=ry[count-1];
                }else{
                    ry[count]=1;
                }
            }
            
            if (CostesPearson<=rmin && i!=LoopMax){
                CostesThrA=i;
                CostesThrB=(int)(a*i+b);
                //i=Amin-1;
            }
            
            rmax=Math.max(rmax,ry[count]);
            rmin=Math.min(rmin,ry[count]);
            count++;
            
            
        }
        
        
        for (int i=0; i<this.length; i++){
            CostesSumA+=this.A[i];
            if (this.A[i]>CostesThrA) CostesSumAThr+=this.A[i];
            CostesSumB+=this.B[i];
            if (this.B[i]>CostesThrB) CostesSumBThr+=this.B[i];
        }
        
        costes_plot=new Plot("Costes' threshold "+this.titleA+" and "+this.titleB,"ThrA", "Pearson's coefficient below",rx,ry);
        costes_plot.setLimits(LoopMin, LoopMax, rmin, rmax);
        costes_plot.setColor(Color.black);
        costes_plot.draw();
        
        //Draw the zero line
        double[] xline={CostesThrA, CostesThrA};
        double[] yline={rmin, rmax};
        costes_plot.setColor(Color.red);
        costes_plot.addPoints(xline, yline, 2);
        
        // OB: Removed output so that we can display the plots programmatically 
        //plot.show();
        
        // OB: Created local variable to be able to request the image as needed.
        costesMask=NewImage.createRGBImage("Costes' mask",this.width,this.height,this.nbSlices,0);
        costesMask.getProcessor().setValue(Math.pow(2, this.depth));
        for (int k=1; k<=this.nbSlices; k++){
        	costesMask.setSlice(k);
            for (int j=0; j<this.height; j++){
                for (int i=0; i<this.width; i++){
                    int position=offset(i,j,k);
                    int [] color=new int[3];
                    color[0]=this.A[position];
                    color[1]=this.B[position];
                    color[2]=0;
                    if (color[0]>CostesThrA && color[1]>CostesThrB){
                        //CostesMask.getProcessor().setValue(((A[position]-CostesThrA)/(LoopMax-CostesThrA))*Math.pow(2, depthA));
                        //CostesMask.getProcessor().drawPixel(i,j);
                        for (int l=0; l<=2; l++) color[l]=255;
                    }
                    costesMask.getProcessor().putPixel(i,j,color);
                }
            }
        }
        costesMask.setCalibration(this.cal);
        costesMask.setSlice(1);
        //costesMask.show();
               
        
        IJ.showStatus("");
        IJ.showProgress(2,1);
        
        this.doThat=true;
        
        IJ.log("\nCostes' automatic threshold set to "+CostesThrA+" for imgA & "+CostesThrB+" for imgB");
        IJ.log("Pearson's Coefficient:\nr="+round(linreg(this.A, this.B,CostesThrA,CostesThrB)[2],3)+" ("+round(CostesPearson,3)+" below thresholds)");
        IJ.log("M1="+round(CostesSumAThr/CostesSumA,3)+" & M2="+round(CostesSumBThr/CostesSumB,3));
        
        // OB: Add Costes to results table 
        rt.addValue("Threshold A", CostesThrA);
        rt.addValue("Threshold B", CostesThrB);
        
        rt.addValue("Thresholded M1", CostesSumAThr/CostesSumA);
        rt.addValue("Thresholded M2", CostesSumBThr/CostesSumB);
        
        // Seeing as it was set, set the thresholds here
        this.thrA = CostesThrA;
        this.thrB = CostesThrB;
        
        
    }
    
    public void CCF(int CCFx){
        double meanA;
        double meanB;
        double nPoints;
        double num;
        double den1;
        double den2;
        double CCF0=0;
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
                CCF0=CCF;
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
        for (int i=0; i<x.length; i++) CCFarray[i]=cf.f(CurveFitter.GAUSSIAN, param, x[i]);
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
                
        //cyto.show();
        //plot.show();
        IJ.log("\nCytofluorogram's parameters:\na: "+round(a,3)+"\nb: "+round(b,3)+"\nCorrelation coefficient: "+round(CoeffCorr,3));

    }
    
    public void ICA(){
        double[] Anorm=new double[this.length];
        double[] Bnorm=new double[this.length];
        double AnormMean=0;
        double BnormMean=0;
        double prodMin=0;
        double prodMax=0;
        double lim=0;
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
       
       if (Math.abs(prodMin)>Math.abs(prodMax)){
           lim=Math.abs(prodMin);
       }else{
           lim=Math.abs(prodMax);
       }
       
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
       //plotA.show();
       
       
       /*double[] xline={0,0};
       double[] yline={0,1};
       plotA.setColor(Color.red);
       plotA.addPoints(xline, yline, Plot.LINE);
        
       plotA.show();*/
       
       ica_b_plot = new Plot("ICA B ("+this.titleB+")", "(Ai-a)(Bi-b)", titleB, new double[]{0, 0}, new double[]{0, 0});
       ica_b_plot.setColor(Color.white);
       ica_b_plot.setLimits(-lim, lim, 0, 1);
       ica_b_plot.draw();
       ica_b_plot.setColor(Color.black);
       ica_b_plot.addPoints(x, Bnorm, Plot.DOT);
       
       
       ica_b_plot.setColor(Color.red);
       ica_b_plot.drawLine(0, 0, 0, 1);
       //plotB.addPoints(xline, yline, Plot.LINE);
        
       //plotB.show();
       
       IJ.log("\nLi's Intensity correlation coefficient:\nICQ: "+ICQ);
       
       // OB: Add CFF to results table
       rt.addValue("ICQ", ICQ);
        
     }
    
    public void CostesRand(int xyBlock, int zBlock, int nbRand, double binWidth, int fillMeth, boolean xyRand, boolean zRand, boolean showRand){
        int[] ACostes, BCostes, BRandCostes;
        
        if (fillMeth==0){
            this.widthCostes=((int)(this.width/xyBlock))*xyBlock;
            this.heightCostes=((int)(this.height/xyBlock))*xyBlock;
        }else{
            this.widthCostes=(((int)(this.width/xyBlock))+1)*xyBlock;
            this.heightCostes=(((int)(this.height/xyBlock))+1)*xyBlock;
        }

        if (zRand){
            if (fillMeth==0){
                this.nbsliceCostes=((int)(this.nbSlices/zBlock))*zBlock;
            }else{
                this.nbsliceCostes=(((int)(this.nbSlices/zBlock))+1)*zBlock;
            }
            if (this.nbSlices==1) nbsliceCostes=1;
            
        }else{
            this.nbsliceCostes=this.nbSlices;
        }
        
        this.lengthCostes=this.widthCostes*this.heightCostes*this.nbsliceCostes;
        ACostes=new int[this.lengthCostes];       
        BCostes=new int[this.lengthCostes];
        BRandCostes=new int[this.lengthCostes];
        
        int index=0;
        for (int k=1; k<=this.nbsliceCostes; k++){
            for (int j=0; j<this.heightCostes; j++){
                for (int i=0; i<this.widthCostes; i++){
                    int offset=offset(i, j, k);
                    ACostes[index]=A[offset];
                    BCostes[index]=B[offset];
                    index++;
                }
            }
        }
        
        
        
         double direction;
         int shift;
         int newposition;
         if (xyRand || this.nbsliceCostes==1){
             //If slices independent 2D there is no need to take into account the z thickness and ranndomization along z axis should not be done
             zBlock=1;
             zRand=false;
         }
         this.doThat=true;
         double r2test=linreg(ACostes, BCostes, 0, 0)[2];
         this.doThat=false;
         double[] arrayR= new double[nbRand];
         double mean=0;
         double SD=0;
         double Pval=0;
         double[] arrayDistribR= new double[(int)(2/binWidth+1)];
         double[] x= new double[arrayDistribR.length];
         
         
         for (int f=0; f<nbRand; f++){
             
             //Randomization by shifting along x axis
             for (int e=1; e<=this.nbsliceCostes-zBlock+1; e+=zBlock){
                 for (int d=0; d<this.heightCostes-xyBlock+1; d+=xyBlock){
                     
                     //Randomization of the shift's direction
                     direction=1;
                     if(Math.random()<0.5) direction=-1;
                     //Randomization of the shift: should be a multiple of the xy block size
                     shift=((int) (direction*Math.random()*this.widthCostes/xyBlock))*xyBlock;
                     
                     for (int a=0; a<this.widthCostes; a++){
                        for (int b=d; b<d+xyBlock; b++){
                            for (int c=e; c<e+zBlock; c++){
                                newposition=a+shift;
                                if (newposition>=this.widthCostes) newposition-=this.widthCostes;
                                if (newposition<0) newposition+=this.widthCostes;
                                BRandCostes[offsetCostes(newposition,b,c)]=BCostes[offsetCostes(a,b,c)];
                            }
                        }
                     }
                 }
             }
             for (int i=0; i<BCostes.length; i++) BCostes[i]=BRandCostes[i];
             
             //Randomization by shifting along y axis
             for (int e=1; e<=this.nbsliceCostes-zBlock+1; e+=zBlock){
                 for (int d=0; d<this.widthCostes-xyBlock+1; d+=xyBlock){
                     
                     //Randomization of the shift's direction
                     direction=1;
                     if(Math.random()<0.5) direction=-1;
                     //Randomization of the shift: should be a multiple of the xy block size
                     shift=((int) (direction*Math.random()*this.heightCostes/xyBlock))*xyBlock;
                     
                     for (int a=0; a<this.heightCostes; a++){
                        for (int b=d; b<d+xyBlock; b++){
                            for (int c=e; c<e+zBlock; c++){
                                newposition=a+shift;
                                if (newposition>=this.heightCostes) newposition-=this.heightCostes;
                                if (newposition<0) newposition+=this.heightCostes;
                                BRandCostes[offsetCostes(b,newposition,c)]=BCostes[offsetCostes(b,a,c)];
                            }
                        }
                     }
                 }
             }
             for (int i=0; i<BCostes.length; i++) BCostes[i]=BRandCostes[i];
             
             if (zRand){
                 //Randomization by shifting along z axis
                 for (int e=0; e<this.heightCostes-xyBlock+1; e+=xyBlock){
                     for (int d=0; d<this.widthCostes-xyBlock+1; d+=xyBlock){

                         //Randomization of the shift's direction
                         direction=1;
                         if(Math.random()<0.5) direction=-1;
                         //Randomization of the shift: should be a multiple of the z block size
                         shift=((int) (direction*Math.random()*this.nbsliceCostes/zBlock))*zBlock;

                         for (int a=1; a<=this.nbsliceCostes; a++){
                            for (int b=d; b<d+xyBlock; b++){
                                for (int c=e; c<e+xyBlock; c++){
                                    newposition=a+shift;
                                    if (newposition>this.nbsliceCostes) newposition-=this.nbsliceCostes;
                                    if (newposition<1) newposition+=this.nbsliceCostes;
                                    BRandCostes[offsetCostes(b,c,newposition)]=BCostes[offsetCostes(b,c,a)];
                                }
                            }
                         }
                     }
                 }
                 for (int i=0; i<BCostes.length; i++) BCostes[i]=BRandCostes[i];
             }
         arrayR[f]=linreg(ACostes, BCostes, 0, 0)[2];
         //if (arrayR[f]<r2test) Pval++;
         mean+=arrayR[f];
         arrayDistribR[(int)((arrayR[f]+1)/binWidth)]++;
         x[(int)((arrayR[f]+1)/binWidth)]+=arrayR[f];
         IJ.showStatus("Costes' randomization loop n�"+f+"/"+nbRand);
         }
         
         //Draw the last randomized image, if required
         if (showRand){
             costes_rand = NewImage.createImage("Randomized images of "+this.titleB,this.widthCostes,this.heightCostes,this.nbsliceCostes,this.depth, 1);
             
             index=0;
             for (int k=1; k<=this.nbsliceCostes; k++){
            	 costes_rand.setSlice(k);
                 for (int j=0;j<this.heightCostes; j++){
                     for (int i=0; i<this.widthCostes;i++){
                    	 costes_rand.getProcessor().putPixel(i, j, BRandCostes[index]);
                         index++;
                     }
                 }
             }
             costes_rand.setCalibration(this.cal);
             costes_rand.setSlice(1);
             //Rand.show();
             //IJ.setMinAndMax(this.Bmin,this.Bmax);
         }
         
         //Plots the r probability distribution
         double minx=-1;
         double maxx=1;
         double maxy=0;
         for (int i=0; i<arrayDistribR.length;i++) x[i]=arrayDistribR[i]==0?i*binWidth-1+binWidth/2:x[i]/arrayDistribR[i];
         for (int i=0; i<arrayDistribR.length;i++) arrayDistribR[i]/=nbRand;
             
         for (int i=0; i<arrayDistribR.length;i++){
             //x[i]=i*binWidth-1+binWidth/2;
             if (minx==-1 && arrayDistribR[i]!=0) minx=x[i];
             if (maxy<arrayDistribR[i]) maxy=arrayDistribR[i];
         }
         minx=Math.min(minx,r2test);
         
         int i=arrayDistribR.length-1;
         while (arrayDistribR[i]==0) {
             maxx=x[i];
             i--;
         }
         
         maxx=Math.max(maxx,r2test);
         
         //Remove from arraDistribR all values equals to zero.
         int newLength=0;
         for (i=0; i<arrayDistribR.length; i++) if(arrayDistribR[i]!=0) newLength++;
         double[] xNew=new double[newLength], arrayNew=new double[newLength];
         newLength=0;
         for (i=0; i<arrayDistribR.length; i++) if(arrayDistribR[i]!=0){ xNew[newLength]=x[i]; arrayNew[newLength++]=arrayDistribR[i];}
         x=xNew;
         arrayDistribR=arrayNew;
         
         
         costes_rand_plot = new Plot("Costes' method ("+this.titleA+" & "+this.titleB+")", "r", "Probability density of r", x, arrayDistribR);
         costes_rand_plot.setLimits(minx-10*binWidth, maxx+10*binWidth, 0, maxy*1.05);
         costes_rand_plot.setColor(Color.white);
         costes_rand_plot.draw();
        
         //Previous plot is white, just to get values inserted into the plot list, the problem being that the plot is as default a line plot... Following line plots same values as circles.
         costes_rand_plot.setColor(Color.black);
         costes_rand_plot.addPoints(x, arrayDistribR, Plot.CIRCLE);
         
        //Draw the r line
         double[] xline={r2test,r2test};
         double[] yline={0,maxy*1.05};
         costes_rand_plot.setColor(Color.red);
         costes_rand_plot.addPoints(xline, yline, 2);
         
         
         //Retrieves the mean, SD and P-value of the r distribution
         for (i=1; i<nbRand; i++) SD+=Math.pow(arrayR[i]-mean,2);
         mean/=nbRand;
         SD=Math.sqrt(SD/(nbRand-1));
         //Pval/=nbRand;
         
         
         IJ.log("\nCostes' randomization based colocalization:\nParameters: Nb of randomization rounds: "+nbRand+", Resolution (bin width): "+binWidth);
         
         
         CurveFitter cf=new CurveFitter(x, arrayDistribR);
         double[] param={0, maxy, mean, SD}; 
         cf.setInitialParameters(param);
         cf.doFit(CurveFitter.GAUSSIAN);
         param=cf.getParams();
         mean=param[2];
         SD=param[3];
         
         //Algorithm 26.2.17 from Abromowitz and Stegun, Handbook of Mathematical Functions for approximation of the cumulative density function (max. error=7.5e^-8). 
         double[] b={0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
         double p=0.2316419;
         double z=(1/Math.sqrt(2*Math.PI))*Math.exp(-Math.pow((r2test-mean)/SD, 2)/2);
         double t=1/(1+p*Math.abs((r2test-mean)/SD));
        
         if(r2test>=0){
             Pval=1-z*t*(t*(t*(t*(t*b[4]+b[3])+b[2])+b[1])+b[0]);
         }else {
             Pval= z*t*(t*(t*(t*(t*b[4]+b[3])+b[2])+b[1])+b[0]);
         }
        
         IJ.log("r (original)="+round(r2test,3)+"\nr (randomized)="+round(mean,3)+"�"+round(SD,3)+" (calculated from the fitted data)\nP-value="+round(Pval*100,2)+"% (calculated from the fitted data)");
         
         IJ.log("\nResults for fitting the probability density function on a Gaussian (Probability=a+(b-a)exp(-(R-c)^2/(2d^2))):"+cf.getResultString()+"\nFWHM="+Math.abs(round(2*Math.sqrt(2*Math.log(2))*param[3], 3)));
         for (i=0; i<x.length; i++) arrayDistribR[i]=cf.f(CurveFitter.GAUSSIAN, param, x[i]);
         costes_rand_plot.setColor(Color.BLUE);
         costes_rand_plot.addPoints(x, arrayDistribR, 2);
         //plot.show();
         
     }
    
    public void distBetweenCentres(int minSize, int maxSize, double limXY, double limZ, boolean cMass, boolean fullList, boolean showImage){
        if (this.countA==null) this.countA=new Counter3D(this.A, this.titleA, this.width, this.height, this.nbSlices, thrA, minSize, maxSize, this.cal);
        if (this.countB==null) this.countB=new Counter3D(this.B, this.titleB, this.width, this.height, this.nbSlices, thrB, minSize, maxSize, this.cal);
        
        double[][] cenA, cenB;
        boolean[] cenAbool, cenBbool;
        int nbColocA=0, nbColocB=0;
        
        if (cMass){
            cenA=this.countA.getCentreOfMassList();
            cenB=this.countB.getCentreOfMassList();
        }else{
            cenA=this.countA.getCentroidList();
            cenB=this.countB.getCentroidList();
        }
        
        String[] header={"Centre_A_ni", "Centre_B_um", "d(A-B)", "reference_dist", "phi", "theta", "XA", "YA", "ZA", "XB", "YB", "ZB"};
        ResultsTable rt=new ResultsTable();
        for (int i=0; i<header.length; i++) rt.setHeading(i, header[i]);
        int index=0;
        
        cenAbool=new boolean[cenA.length];
        cenBbool=new boolean[cenB.length];
        
        for (int i=0; i<cenAbool.length; i++) cenAbool[i]=false;
        for (int i=0; i<cenBbool.length; i++) cenBbool[i]=false;
        
        for (int i=0; i<cenA.length; i++){
            for (int j=0; j<cenB.length; j++){
                double x=(cenB[j][0]-cenA[i][0])*this.cal.pixelWidth;
                double y=(cenB[j][1]-cenA[i][1])*this.cal.pixelHeight;
                double z=(cenB[j][2]-cenA[i][2])*this.cal.pixelDepth;
                
                double distXY=Math.sqrt(x*x+y*y);
                double distXYZ=Math.sqrt(distXY*distXY+z*z);
                
                /*
                 *The first Airy disc in 3D is not a sphere but rather egg shaped. Therefore, while the maximimum ditance between two colocalising spots in 2D is equal to the xy optical resolution 
                 *its hard to figure it out along a xz section as the cross section is an ellipse rather than a sphere. What if this section is not coincident with the equatorial plane ?!!!
                 *The only mean is to calculate the distance on the Airy "egg shape"...
                 *First, we convert the system: centre A becomes the origin of the spherical space. Then we calculate the two coordinates of B into the new space (phi, theta) ie angles in reference
                 *to axis Z and X.
                */
                
                double theta=0;
                if (distXYZ!=0) theta=Math.acos(z/distXYZ);
                
                double phi=Math.PI/2;
                if (distXY!=0) phi=Math.acos(x/distXY);
                
                /*
                 *Second, we use the two angles in the equation of the "egg shape" to estimate the coordinates of the pixel on its border. Then, we calculate the distance between the origin and this
                 *pixel: it will be used as the reference distance...
                */
                
                double xRef=limXY*Math.sin(theta)*Math.cos(phi);
                double yRef=limXY*Math.sin(theta)*Math.sin(phi);
                double zRef=limZ*Math.cos(theta);
                
                double distRef=Math.sqrt(xRef*xRef+yRef*yRef+zRef*zRef);
                
                if (distXYZ<=distRef || fullList){
                    if (distXYZ<=distRef){
                        cenAbool[i]=true;
                        cenBbool[j]=true;
                    }
                    
                    rt.incrementCounter();
                    rt.setValue("Centre_A_ni", index, i+1);
                    rt.setValue("Centre_B_ni", index, j+1);
                    if (fullList) rt.setLabel(distXYZ<=distRef?"Colocalization":"No colocalization", index);
                    rt.setValue("d(A-B)", index, distXYZ);
                    rt.setValue("reference_dist", index, distRef);
                    rt.setValue("theta", index, theta);
                    rt.setValue("phi", index, phi);
                    rt.setValue("XA", index, cenA[i][0]);
                    rt.setValue("YA", index, cenA[i][1]);
                    if (this.nbSlices>1) rt.setValue("ZA", index, cenA[i][2]);
                    rt.setValue("XB", index, cenB[j][0]);
                    rt.setValue("YB", index, cenB[j][1]);
                    if (this.nbSlices>1) rt.setValue("ZB", index, cenB[j][2]);
                    index++;
                }
            }
        }
        
        if (rt.getCounter()==0){
            rt.incrementCounter();
            rt.addLabel("Result", "No colocalization found");
        }
        
        for (int i=0; i<cenAbool.length; i++) if (cenAbool[i]) nbColocA++;
        for (int i=0; i<cenBbool.length; i++) if (cenBbool[i]) nbColocB++;
        
        String title="Distance based colocalization between "+this.titleA+ " and "+this.titleB+(cMass?" (centres of mass)":" (geometrical centres)");
        
        if (showImage){
            centers_img  = NewImage.createImage(title, this.width, this.height, this.nbSlices, 24, 1);
            for (int i=0; i<cenA.length; i++){
                if (cenAbool[i] || fullList){
                	centers_img.setSlice((int) cenA[i][2]);
                    int[] val={255, 0, 0};
                    if (cenAbool[i]) val[2]=255;
                    centers_img.getProcessor().putPixel((int) cenA[i][0], (int) cenA[i][1], val);
                }
            }
            
            //centers_img.show();
            
            for (int i=0; i<cenB.length; i++){
                if (cenBbool[i] || fullList){
                	centers_img.setSlice((int) cenB[i][2]);
                    int[] val=centers_img.getPixel((int) cenB[i][0], (int) cenB[i][1]);
                    val[1]=255;
                    if (cenBbool[i]) val[2]=255;
                    if (val[0]==255 && val[1]==255) val[2]=0;
                    centers_img.getProcessor().putPixel((int) cenB[i][0], (int) cenB[i][1], val);
                }
            }
            centers_img.setCalibration(this.micronCal);
            centers_img.getProcessor().resetMinAndMax();
            //img.show();
           //centers_img.updateAndDraw();
        }
        
        rt.show(title);
        IJ.log("\nColocalization based on distance between "+(cMass?"centres of mass":"geometrical centres"));
        IJ.log("Threshold for Image A="+thrA+"; Image B="+thrB);
        IJ.log("Particles size between "+minSize+" & "+maxSize);
        IJ.log("Image A: "+nbColocA+" centre(s) colocalizing out of "+cenA.length);
        IJ.log("Image B: "+nbColocB+" centre(s) colocalizing out of "+cenB.length);
        
    }
    
    public void coincidenceCentreParticle(int minSize, int maxSize, boolean cMass, boolean fullList, boolean showImage){
        if (this.countA==null) this.countA=new Counter3D(this.A, this.titleA, this.width, this.height, this.nbSlices, thrA, minSize, maxSize, this.cal);
        if (this.countB==null) this.countB=new Counter3D(this.B, this.titleB, this.width, this.height, this.nbSlices, thrB, minSize, maxSize, this.cal);
        
        String title=(cMass?" (Centres of mass)":" (Geometrical centres)")+" of "+this.titleA+"-Particles of "+this.titleB+" based colocalization";
        ResultsTable rt=new ResultsTable();
        String[] header={"Centre_"+this.titleA+"_ni", "Particle_"+this.titleB+"_ni", "X", "Y", "Z"};
        for (int i=0; i<header.length; i++) rt.setHeading(i, header[i]);
        
        Object3D[] objCentres=this.countA.getObjectsList();
        Object3D[] objParticles=this.countB.getObjectsList();
        
        boolean[] centBool=new boolean[objCentres.length];
        boolean[] partBool=new boolean[objParticles.length];
        int index=0;
        int nbColocA=0, nbColocB=0;
        
        for (int i=0; i<objCentres.length; i++){
            double[] currCent=cMass?(objCentres[i].c_mass):(objCentres[i].centroid);
            for (int j=0; j<objParticles.length; j++){
                
                int[] boundTL=objParticles[j].bound_cube_TL;
                int[] boundBR=objParticles[j].bound_cube_BR;
                boolean isColoc=false;
                 
                boolean isProbablyColoc=currCent[0]>=boundTL[0] && currCent[0]<=boundBR[0] && currCent[1]>=boundTL[1] && currCent[1]<=boundBR[1] && currCent[2]>=boundTL[2] && currCent[2]<=boundBR[2];
                
                if (isProbablyColoc){
                    for (int k=0; k<objParticles[j].size; k++){
                        isColoc=(int) currCent[0]==objParticles[j].obj_voxels[k][0] && (int) currCent[1]==objParticles[j].obj_voxels[k][1] && (int) currCent[2]==objParticles[j].obj_voxels[k][2];
                        if (isColoc) k=objParticles[j].size;
                    }
                }
                
                if (isColoc || fullList){
                    if (isColoc){
                        centBool[i]=true;
                        partBool[j]=true;
                    }
                    
                    rt.incrementCounter();
                    rt.setValue("Centre_"+this.titleA+"_ni", index, i+1);
                    rt.setValue("Particle_"+this.titleB+"_ni", index, j+1);
                    if (fullList) rt.setLabel(isColoc?"Colocalization":"No colocalization", index);
                    rt.setValue("X", index, currCent[0]);
                    rt.setValue("Y", index, currCent[1]);
                    rt.setValue("Z", index, currCent[2]);
                    index++;
                }
            }
        }
        
        if (rt.getCounter()==0){
            rt.incrementCounter();
            rt.addLabel("Result", "No colocalization found");
        }
        
        for (int i=0; i<centBool.length; i++) if (centBool[i]) nbColocA++;
        for (int i=0; i<partBool.length; i++) if (partBool[i]) nbColocB++;
        
        if (showImage){
            coinc_img=NewImage.createImage(title, this.width, this.height, this.nbSlices, 24, 1);
            for (int i=0; i<objParticles.length; i++){
                if (partBool[i] || fullList){
                    for (int j=0; j<objParticles[i].size; j++){
                    	coinc_img.setSlice(objParticles[i].obj_voxels[j][2]);
                        int[] val={255, 0, 0};
                        coinc_img.getProcessor().putPixel(objParticles[i].obj_voxels[j][0], objParticles[i].obj_voxels[j][1], val);
                    }
                }
            }
            
            //coinc_img.show();
            
            for (int i=0; i<objCentres.length; i++){
                if (centBool[i] || fullList){
                    double[] currCent=cMass?(objCentres[i].c_mass):(objCentres[i].centroid);
                    coinc_img.setSlice((int) currCent[2]);
                    int[] val=coinc_img.getPixel((int) currCent[0], (int) currCent[1]);
                    val[1]=255;
                    coinc_img.getProcessor().putPixel((int) currCent[0], (int) currCent[1], val);
                }
            }
            
            coinc_img.setCalibration(this.micronCal);
            coinc_img.getProcessor().resetMinAndMax();
            //coinc_img.updateAndDraw();
        }
        
        rt.show(title);
        IJ.log("\nColocalization based on "+(cMass?"centres of mass":"geometrical centres")+"-particles coincidence");
        IJ.log("Threshold for Image A="+thrA+"; Image B="+thrB);
        IJ.log("Particles size between "+minSize+" & "+maxSize);
        IJ.log("Image A: "+nbColocA+" centre(s) colocalizing out of "+this.countA.nbObj);
        IJ.log("Image B: "+nbColocB+" centre(s) colocalizing out of "+this.countB.nbObj);
    }
    
    //----------------------------------------------------------------------------------------------------------------------------------------------
    private void buildArray(ImagePlus imgA, ImagePlus imgB){
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
                    index++;
                }
            }
            
            this.Amean/=this.length;
            this.Bmean/=this.length;
        }
    }
    
    /** Generates the ImagePlus base on the input array and title.
     * @param array containing the pixels intensities (integer array).
     * @param title to attribute to the ImagePlus (string).
     */
    private ImagePlus buildImg(int[] array, String title){
        int index=0;
        double min=array[0];
        double max=array[0];
        ImagePlus img=NewImage.createImage(title, this.width, this.height, this.nbSlices, this.depth, 1);
        
        for (int z=1; z<=this.nbSlices; z++){
            IJ.showStatus("Creating the image...");
            img.setSlice(z);
            for (int y=0; y<this.height; y++){
                for (int x=0; x<this.width; x++){
                    int currVal=array[index];
                    min=Math.min(min, currVal);
                    max=Math.max(max, currVal);
                    img.getProcessor().putPixel(x,y, currVal);
                    index++;
                }
            }
        }
        IJ.showStatus("");
        img.setCalibration(this.micronCal);
        img.getProcessor().setMinAndMax(min, max);
        return img;
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
                if (Aarray[m]>=TA && Barray[m]>=TB){
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
            if (Aarray[m]>=TA && Barray[m]>=TB){
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
    
    public double[] linregCostes(int[] Aarray, int[] Barray, int TA, int TB){
         double num=0;
         double den1=0;
         double den2=0;
         double[] coeff=new double[3];
         int count=0;
         
         sumA=0;
         sumB=0;
         sumAB=0;
         sumsqrA=0;
         Aarraymean=0;
         Barraymean=0;
         
         for (int m=0; m<Aarray.length; m++){
            if (Aarray[m]<TA && Barray[m]<TB){
                sumA+=Aarray[m];
                sumB+=Barray[m];
                sumAB+=Aarray[m]*Barray[m];
                sumsqrA+=Math.pow(Aarray[m],2);
                count++;
            }
        }

             Aarraymean=sumA/count;
             Barraymean=sumB/count;
                  
         
         for (int m=0; m<Aarray.length; m++){
            if (Aarray[m]<TA && Barray[m]<TB){
                num+=(Aarray[m]-Aarraymean)*(Barray[m]-Barraymean);
                den1+=Math.pow((Aarray[m]-Aarraymean), 2);
                den2+=Math.pow((Barray[m]-Barraymean), 2);
            }
         }
        
        coeff[0]=(count*sumAB-sumA*sumB)/(count*sumsqrA-Math.pow(sumA,2));
        coeff[1]=(sumsqrA*sumB-sumA*sumAB)/(count*sumsqrA-Math.pow(sumA,2));
        coeff[2]=num/(Math.sqrt(den1*den2));
        return coeff;
     }
    
    private double[] int2double(int[] input){
        double[] output=new double[input.length];
        for (int i=0; i<input.length; i++) output[i]=input[i];
        return output;
    }
    
    /** Returns the index where to find the informations corresponding to pixel (x, y, z).
     * @param x coordinate of the pixel.
     * @param y coordinate of the pixel.
     * @param z coordinate of the pixel.
     * @return the index where to find the informations corresponding to pixel (x, y, z).
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
    
    
	public ImagePlus getFluorogramImage(int min, int max) {

		Plot fp = getFluorogram();

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
		fluo_ip.setLut(ColocOutput.fireLUT());
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
    
    // OB: Extra functions to return the results
    public void showResults() {
    	rt.show("Results");
    }

	public Plot getICAPlot() {
		return icq_plot;
	}

	public Plot getCostesPlot() {
		return costes_plot;
	}

	public Plot getCFFPlot() {
		return cff_plot;
	}

	public ImagePlus getCostesMask() {
		return costesMask;
	}

	public Plot getICAaPlot() {
		return ica_a_plot;
	}

	public Plot getICAbPlot() {
		return ica_b_plot;
	}

	public ImagePlus getCostesRandImg() {
		return costes_rand;
	}

	public Plot getCostes_rand_plot() {
		return costes_rand_plot;
	}

	public ImagePlus getCentersImg() {
		return centers_img;
	}

	public ImagePlus getCoincidentImg() {
		return coinc_img;
	}
	
	public Plot getFluorogram() {
		return fluorogram_plot;
	}
	
	public void setThresholds(int thrA, int thrB) {
		
		this.thrA = thrA;
		this.thrB = thrB;
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
		double[][] rangeA = getDisplayRange(impA);
		double[][] rangeB = getDisplayRange(impB);
  		showDisplayRange(impA);
  		
		// Allow for thresholds to be either manual or automatic
		if(thrMetA.matches("\\d*")) {
			this.thrA = Integer.valueOf(thrMetA);
			impA.getProcessor().setThreshold(this.thrA , impA.getProcessor().getMaxThreshold(), ImageProcessor.NO_LUT_UPDATE);
		
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
			IJ.log("Threshold "+thrMetB+" matches regex");
			this.thrB = Integer.valueOf(thrMetB);
			impB.getProcessor().setThreshold(this.thrB , impB.getProcessor().getMaxThreshold(), ImageProcessor.NO_LUT_UPDATE);
			
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
				
		//impA.killRoi();
		//impB.killRoi();
		
		rt.addValue("Auto Threshold A", thrMetA);
		rt.addValue("Auto Threshold B", thrMetB);
		
		setDisplayRange(impA, rangeA);
		setDisplayRange(impB, rangeB);

				
	}
	
	public ImagePlus getImageA() {
		return impA;
	}

	public ImagePlus getImageB() {
		return impB;
	}
	
	public ImagePlus getRGBImage(ImagePlus imp, Boolean is_zProject) {
		
		imp.killRoi();
		ImagePlus impr = imp.duplicate();
		impr.setTitle(imp.getTitle());
			
		if(is_zProject) {
			ZProjector zp = new ZProjector(impr);
			zp.setMethod(ZProjector.MAX_METHOD);
			zp.doHyperStackProjection(true);
			impr = zp.getProjection();
		}
		
		return flattenRoi(impr);
	}

	public ImagePlus getRGBImageA(Boolean is_Z) {
		return getRGBImage(impA, is_Z);
	}
	
	public ImagePlus getRGBImageB(Boolean is_Z) {
		return getRGBImage(impB, is_Z);
	}
	
	public ImagePlus getRGBColocImage() {
		// This should return an rgb image of the composite of the two channels, with the pixels as a mask
		// Make a composite of the two images
		ImagePlus[] images = {impA, impB};
		
		ImagePlus comp = RGBStackMerge.mergeChannels(images, true);
		return flattenRoi(comp);
		
	}
	
	public ImagePlus getMaskA() {
		return binarize(impA, thrA);
		
	}
	
	public ImagePlus getMaskB() {
		return binarize(impB, thrB);
	}
	
	public ImagePlus getRGBMaskA() {
		return flattenRoi(binarize(impA, thrA));
		
	}
	
	public ImagePlus getRGBMaskB() {
		return flattenRoi(binarize(impB, thrB));
		
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
	
	public ImagePlus getRGBANDMask() {
		return flattenRoi(getANDMask());
	}
	
	
	
	// Proper flattening of the ROI onto the image or stack
	private ImagePlus flattenRoi(ImagePlus imp) {
		if(roi != null) {
			roi.setStrokeWidth(2);
			roi.setStrokeColor(roiColor);
			imp.setOverlay(roi, roiColor, 2, null);
			imp.setHideOverlay(false);
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
	
	
	// Seeing as binarizing images properly is a bitch, let's do it ourselves...
	private ImagePlus binarize(ImagePlus imp, int lowerThr) { 
		
		if (imp.getStackSize() == 1) {
			ImagePlus imp2 = new ImagePlus("Mask "+imp.getTitle(), binarize(imp.getProcessor(), lowerThr));
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
			stack2.addSlice(label, binarize(ip, lowerThr));
		}
		
		ImagePlus imp2 = new ImagePlus("Mask "+imp.getTitle(), stack2);
		imp2.setCalibration(imp.getCalibration()); //update calibration
		
		return imp2;
	}
	
	// Binarize image processor
	private ImageProcessor binarize(ImageProcessor ip, int lowerThr) {
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
	
	public void addResult(String name, int value) {
		rt.addValue(name, value);
	}
	
	public void addResult(String name, String value) {
		rt.addValue(name, value);
	}
	
	public void removeLastRow() {
		rt.deleteRow(rt.getCounter()-1);
	}

	private double[][] getDisplayRange(ImagePlus imp) {
		double[][] drange = new double[imp.getNChannels()][];
		
		for(int i=0;i<drange.length;i++) {
			imp.setC(i+1);
			drange[i] = new double[2];
			drange[i][0] = imp.getDisplayRangeMin();
			drange[i][1] = imp.getDisplayRangeMax();
		}
		return drange;
	}
	
	private void setDisplayRange(ImagePlus impr, double[][] range) {
		// TODO Auto-generated method stub
		for(int i=0; i<range.length;i++) {
			impr.setC(i+1);
			impr.setDisplayRange(range[i][0], range[i][1]);
		}
	}
	private void showDisplayRange(ImagePlus imp) {
		for(int i=0; i<imp.getNChannels();i++)
		{
			imp.setC(i+1);
			IJ.log("Channel "+(i+1)+" MIN: "+imp.getDisplayRangeMin()+", MAX: "+imp.getDisplayRangeMax()+" For Image "+imp.getTitle());
		}
	}

	
}