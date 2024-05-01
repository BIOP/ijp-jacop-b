/*-
 * #%L
 * An update to the JACoP Plugin that helps in the management of ROIs, Z sections and helps generate cleaner reports
 * %%
 * Copyright (C) 2009 - 2023 Susanne Bolte, Fabrice P. Cordelières, Olivier Burri
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

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.CombinatoricsUtils;

/**
 * @author Nicolas Chiaruttini
 */

public class RandomCostes {

    public ImagePlus imp_orig;
    public ImagePlus imgA;
    public ImagePlus imgB;
    public int squareSize;
    public int nShuffling;
    public boolean binarize;
    public int thrA, thrB;
    public Roi roi;

    public double pearsonNormalized = Double.NaN;
    public double pearson = Double.NaN;

    public double pValueIsCorrelated = Double.NaN;
    public double pValueIsAntiCorrelated = Double.NaN;

    public RandomCostes(ImagePlus imgA, ImagePlus imgB, int squareSize, int nShuffling, boolean binarize, int thrA, int thrB) {

        this.imgA=imgA.duplicate();
        this.imgB=imgB.duplicate();

        roi = imgA.getRoi();
        if (roi==null) {
            roi = new Roi(0,0,imgA.getWidth(), imgA.getHeight());
        }
        this.squareSize = squareSize;
        this.nShuffling=nShuffling;
        this.binarize=binarize;

        if (binarize) {
            imgA.getProcessor().setThreshold(thrA,Double.MAX_VALUE,ImageProcessor.NO_LUT_UPDATE);
            ImagePlus imgTempA = new ImagePlus();
            imgTempA.setProcessor(imgA.getProcessor().createMask());
            this.imgA = imgTempA;

            imgB.getProcessor().setThreshold(thrB,Double.MAX_VALUE,ImageProcessor.NO_LUT_UPDATE);
            ImagePlus imgTempB = new ImagePlus();
            imgTempB.setProcessor(imgB.getProcessor().createMask());
            this.imgB = imgTempB;
        }

        imp_orig = RGBStackMerge.mergeChannels(new ImagePlus[]{imgA,imgB},true);

        this.thrA=thrA;
        this.thrB=thrB;
    }

    public void compute() {
        // Copy image with roi and clear what's outside of the ROI
        imp_orig.setDisplayMode(IJ.GRAYSCALE);
        //Roi roi = imp_orig.getRoi();
        imp_orig.setRoi((Roi) null);
        ImagePlus imp = imp_orig.duplicate();
        imp.setRoi(roi);

        ShapeRoi sroi = new ShapeRoi(roi);

        Shape shape = sroi.getShape();

        Rectangle bounds = shape.getBounds();

        int nBlocks=0;

        ImagePlus impCH1 = IJ.createImage("Costes Block CH1", "32-bit black", squareSize, squareSize, 1);
        ImagePlus impCH2 = IJ.createImage("Costes Block CH2", "32-bit black", squareSize, squareSize, 1);

        for (int x=bounds.x;x+squareSize<bounds.x+bounds.width;x+=squareSize) {
            for (int y=bounds.y;y+squareSize<bounds.y+bounds.height;y+=squareSize) {
                Rectangle r = new Rectangle(x,y,squareSize, squareSize);
                if (sroi.getShape().contains(r)) {

                    // Copy block to stack
                    imgA.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(imgA, "Copy", "");
                    IJ.run(impCH1, "Add Slice", "");
                    impCH1.setSlice(nBlocks);
                    IJ.run(impCH1, "Paste", "");

                    imgB.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(imgB, "Copy", "");
                    IJ.run(impCH2, "Add Slice", "");
                    impCH2.setSlice(nBlocks);
                    IJ.run(impCH2, "Paste", "");

                    nBlocks++;
                }
            }
        }

        double pearson = getPearson(impCH1, impCH2);

        // Normal indexes
        ArrayList<Integer> blockIndexes = new ArrayList<>();

        for (int i = 0;i<nBlocks ; i++) {
            blockIndexes.add(i);
        }

        imp.setRoi(roi);

        double maxShuffling = CombinatoricsUtils.factorialDouble(nBlocks);
        //System.out.println("Maximal number of shuffling = "+maxShuffling)

        if (nShuffling<100) {
            System.out.println("Take care! Low number of shuffling...");
        }

        if (nShuffling>2*maxShuffling) {
            System.out.println("Take care! You ask for "+nShuffling+" shufflings, but the maximal number of shuffling is "+maxShuffling);
            System.out.println("You are sampling many times the same shuffling.");
        }

        sd = new StandardDeviation(false);
        double sum = 0;
        int nAbove = 0;
        int nBelow = 0;
        valuesShuffling = new double[nShuffling];

        ImagePlus impCH2Shuffled = IJ.createImage("Costes Block CH2 Shuffled", "32-bit black", squareSize, squareSize, nBlocks);

        for (int iShuffle=0;iShuffle<nShuffling;iShuffle++) {
            Collections.shuffle(blockIndexes);
            for (int iSlice=0;iSlice<nBlocks;iSlice++) {
                impCH2Shuffled.getStack().setProcessor(impCH2.getStack().getProcessor(blockIndexes.get(iSlice)+1), iSlice+1);
            }

            double value = getPearson(impCH1, impCH2Shuffled);

            if (value>pearson) {
                nAbove=nAbove+1;
            }
            if (value<pearson) {
                nBelow=nBelow+1;
            }
            sum+=value;
            sd.increment(value);
            valuesShuffling[iShuffle] = value;
        }

        double shufflingMean = (sum/(double)nShuffling);
        double shufflingStd = sd.getResult();
        this.pearsonNormalized = (pearson-shufflingMean) / shufflingStd;
        this.pearson = pearson;
        this.pValueIsCorrelated = Erf.erfc(-this.pearsonNormalized)/2.0;
        this.pValueIsAntiCorrelated = Erf.erfc(this.pearsonNormalized)/2.0;

        imp.changes=false;

        ImagePlus impSampleShuffleCH1 = IJ.createImage("Sample CH1", "32-bit black", imp.getWidth(), imp.getHeight(), 1);
        ImagePlus impSampleShuffleCH2 = IJ.createImage("Sample Shuffle CH2", "32-bit black", imp.getWidth(), imp.getHeight(), 1);

        int iSlice=0;

        for (int x=bounds.x;x+squareSize<bounds.x+bounds.width;x+=squareSize) {
            for (int y=bounds.y;y+squareSize<bounds.y+bounds.height;y+=squareSize) {
                Rectangle r = new Rectangle(x,y,squareSize, squareSize);
                if (sroi.getShape().contains(r)) {

                    // Copy block to stack
                    imgA.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(imgA, "Copy", "");
                    impSampleShuffleCH1.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(impSampleShuffleCH1, "Paste", "");

                    iSlice++;
                    impCH2Shuffled.setSlice(iSlice);
                    IJ.run(impCH2Shuffled, "Copy", "");
                    impSampleShuffleCH2.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(impSampleShuffleCH2, "Paste", "");

                    nBlocks++;
                }
            }
        }

        impSampleShuffleCH1.setLut(imgA.getLuts()[0]);

        impSampleShuffleCH2.setLut(imgB.getLuts()[0]);



        imp.close();

        impSampleShuffle = RGBStackMerge.mergeChannels(new ImagePlus[]{impSampleShuffleCH1,impSampleShuffleCH2},true);

    }
    ImagePlus impSampleShuffle;

    StandardDeviation sd;
    double[] valuesShuffling;

    public Plot getPearsonDistributionGraph(boolean userDefinedLimits, double userMin, double userMax) {

        double binRes = 0.25;

        int minBin = (int) (-6/binRes);
        int maxBin = (int) (+6/binRes);

        if (Double.isNaN(pearsonNormalized)) {
        } else if (pearsonNormalized == Double.NEGATIVE_INFINITY) {
        } else if (pearsonNormalized == Double.POSITIVE_INFINITY) {
        } else if (pearsonNormalized>0) {
            maxBin = Math.max(maxBin, (int) (pearsonNormalized*1.5/binRes+1) );
        } else if (pearsonNormalized<0) {
            minBin = Math.min(minBin, (int) (pearsonNormalized*1.5/binRes-1) );
        }

        double[] bins = new double[maxBin-minBin];
        double[] xs = new double[maxBin-minBin];
        double[] gaussFit = new double[maxBin-minBin];

        for (int iV=0;iV<valuesShuffling.length;iV++) {
            int iBin = (int) (  ((valuesShuffling[iV]/(sd.getResult()))/binRes-minBin)+0.5);
            if (iBin<0) {
                bins[0]++;
            } else if (iBin>=bins.length) {
                bins[bins.length-1]++;
            } else {
                bins[iBin]++;
            }
        }

        double normFactor = 1/Math.sqrt(2*Math.PI);
        for (int i=0;i<xs.length;i++) {
            xs[i]=(minBin+i)*binRes*sd.getResult();
            bins[i]/=(valuesShuffling.length);
            bins[i]/=binRes;
            gaussFit[i]= Math.exp(-xs[i]*xs[i]/(2*sd.getResult()*sd.getResult()))*normFactor;
        }

        String titlePlot = "Random Costes";
        if (binarize) {titlePlot+=" Mask";}

        Plot plot = new Plot(titlePlot,"Pearson","Probability");

        plot.setColor("black");
        plot.setLineWidth(3);
        plot.addPoints(xs,gaussFit,Plot.LINE);
        plot.addPoints(xs,bins, Plot.CIRCLE);

        plot.setColor("red");
        plot.addPoints(new double[]{pearson,pearson},new double[]{0,0.4},Plot.LINE );
        //plot.setLimits(plot.getLimits()[0], plot.getLimits()[1], 0,0.4);
        if (userDefinedLimits) {
            plot.setLimits(userMin, userMax, 0,0.4);
        }
        return plot;
    }

    public ImagePlus getExampleShuffleImage() {
        return impSampleShuffle;
    }

    // Deeply inspired from https://github.com/fiji/Colocalisation_Analysis/blob/master/src/main/java/sc/fiji/coloc/Colocalisation_Threshold.java

    double getPearson(ImagePlus img1, ImagePlus img2){

        float ch1, ch2;
        int nSlices = img1.getNSlices();
        int rheight = img1.getHeight();
        int rwidth = img1.getWidth();

        double N = nSlices*rheight*rwidth;

        double sumX=0, sumXX=0, sumXY=0, sumYY=0, sumY=0;

        //calulate variances
        for (int s=1; s<nSlices; s++) {
            ImageProcessor ip1 = img1.getStack().getProcessor(s);
            ImageProcessor ip2 = img2.getStack().getProcessor(s);
            for (int y=0; y<rheight; y++) {
                for (int x=0; x<rwidth; x++) {
                    ch1 = ip1.getPixelValue(x,y);
                    ch2 = ip2.getPixelValue(x,y);
                    sumX = sumX+ch1;
                    sumXY = sumXY + (ch1 * ch2);
                    sumXX = sumXX + (ch1 * ch1);
                    sumYY = sumYY + (ch2 * ch2);
                    sumY = sumY + ch2;
                }
            }
        }

        double pearsons1 = sumXY - (sumX*sumY/N);
        double pearsons2 = sumXX - (sumX*sumX/N);
        double pearsons3 = sumYY - (sumY*sumY/N);
        double rTotal = pearsons1/(Math.sqrt(pearsons2*pearsons3));

        return rTotal;
    }

}
