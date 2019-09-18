package ch.epfl.biop.coloc.utils;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.awt.Rectangle;
import java.awt.Shape;
import java.util.Collections;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.CombinatoricsUtils;

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

    public RandomCostes(ImagePlus imgA, ImagePlus imgB, int squareSize, int nShuffling, boolean binarize, int thrA, int thrB) {

        this.imgA=imgA.duplicate();
        this.imgB=imgB.duplicate();



        roi = imgA.getRoi();
        this.squareSize = squareSize;
        this.nShuffling=nShuffling;
        this.binarize=binarize;

        if (binarize) {
            //System.err.println("2D BIOP Random Costes with Binarization do not work yet");
            imgA.getProcessor().setThreshold(thrA,Double.MAX_VALUE,ImageProcessor.NO_LUT_UPDATE);
            ImagePlus imgTempA = new ImagePlus();
            imgTempA.setProcessor(imgA.getProcessor().createMask());
            imgA = imgTempA;

            imgB.getProcessor().setThreshold(thrB,Double.MAX_VALUE,ImageProcessor.NO_LUT_UPDATE);
            ImagePlus imgTempB = new ImagePlus();
            imgTempB.setProcessor(imgB.getProcessor().createMask());
            imgB = imgTempB;


        }
        //imgA.show();
        //imgB.show();

        imp_orig = RGBStackMerge.mergeChannels(new ImagePlus[]{imgA,imgB},true);
        //imp_orig.show();

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

                    //rm.addRoi(new Roi(x+sroi.getXBase(),y+sroi.getYBase(),squareSize, squareSize)); // Uncomment to vizualize blocks
                    // Copy block to stack
                    //imp.setC(0);
                    imgA.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(imgA, "Copy", "");
                    IJ.run(impCH1, "Add Slice", "");
                    impCH1.setSlice(nBlocks);
                    IJ.run(impCH1, "Paste", "");

                    //imp.setC(1);

                    imgB.setRoi((int) (x+sroi.getXBase()),(int) (y+sroi.getYBase()),squareSize, squareSize);
                    IJ.run(imgB, "Copy", "");
                    IJ.run(impCH2, "Add Slice", "");
                    impCH2.setSlice(nBlocks);
                    IJ.run(impCH2, "Paste", "");

                    nBlocks++;
                }
            }
        }

        double im1Avg = getMean(impCH1);
        double im2Avg = getMean(impCH2);
        double im3Avg = im1Avg+im2Avg;

        double pearson = getPearson(impCH1, impCH2, im1Avg, im2Avg, im3Avg);

        //System.out.println("pearson= "+pearson );

        // Normal indexes
        ArrayList<Integer> blockIndexes = new ArrayList<>();

        for (int i = 0;i<nBlocks ; i++) {
            blockIndexes.add(i);
        }

        imp.setRoi(roi);

        //if (printWarning) {
            double maxShuffling = CombinatoricsUtils.factorialDouble(nBlocks);
            //System.out.println("Maximal number of shuffling = "+maxShuffling)

            if (nShuffling<100) {
                System.out.println("Take care! Low number of shuffling...");
            }

            if (nShuffling>2*maxShuffling) {
                System.out.println("Take care! You ask for "+nShuffling+" shufflings, but the maximal number of shuffling is "+maxShuffling);
                System.out.println("You are probably sampling many times the same shuffling");
            }
        //}

        StandardDeviation sd = new StandardDeviation(false);
        double sum = 0;
        int nAbove = 0;
        int nBelow = 0;

        ImagePlus impCH2Shuffled = IJ.createImage("Costes Block CH2 Shuffled", "32-bit black", squareSize, squareSize, nBlocks);

        im1Avg = getMean(impCH1);
        im2Avg = getMean(impCH2);
        im3Avg = im1Avg+im2Avg;

        for (int iShuffle=0;iShuffle<nShuffling;iShuffle++) {
            Collections.shuffle(blockIndexes);
            for (int iSlice=0;iSlice<nBlocks;iSlice++) {
                impCH2Shuffled.getStack().setProcessor(impCH2.getStack().getProcessor(blockIndexes.get(iSlice)+1), iSlice+1);
            }

            double value = getPearson(impCH1, impCH2Shuffled, im1Avg, im2Avg, im3Avg);
            //System.out.println(value);

            if (value>pearson) {
                nAbove=nAbove+1;
            }
            if (value<pearson) {
                nBelow=nBelow+1;
            }
            sum+=value;
            sd.increment(value);
        }

        double shufflingMean = (sum/(double)nShuffling);
        double shufflingStd = sd.getResult();
        pearsonNormalized = (pearson-shufflingMean) / shufflingStd;

        imp.changes=false;
        imp.close();

    }

    // Deeply inspired from https://github.com/fiji/Colocalisation_Analysis/blob/master/src/main/java/sc/fiji/coloc/Colocalisation_Threshold.java

    public static double getMean(ImagePlus img) {
        int nSlices = img.getNSlices();
        int rheight = img.getHeight();
        int rwidth = img.getWidth();

        double N = nSlices*rheight*rwidth;
        double chSum = 0;

        for (int s=1; s<=nSlices; s++) {
            ImageProcessor ip = img.getStack().getProcessor(s);
            for (int y=0; y<rheight; y++) {
                for (int x=0; x<rwidth; x++) {
                    double ch = ip.getPixelValue(x,y);
                    chSum+=ch;
                }
            }
        }
        return chSum/N;
    }


    double getPearson(ImagePlus img1, ImagePlus img2, double ch1Mean, double ch2Mean, double ch3Mean ) {

        float ch1 = 0, ch2 = 0, ch3 = 0;
        double ch1Sum = 0, ch2Sum = 0, ch3Sum = 0;

        int nSlices = img1.getNSlices();
        int rheight = img1.getHeight();
        int rwidth = img1.getWidth();

        double N = nSlices*rheight*rwidth;

	/*for (int s=1; s<=nSlices; s++) {
		ip1 = img1.getStack().getProcessor(s);
		ip2 = img2.getStack().getProcessor(s);
		for (int y=0; y<rheight; y++) {
			for (int x=0; x<rwidth; x++) {
				ch1 = ip1.getPixelValue(x,y);
				ch2 = ip2.getPixelValue(x,y);
				ch3 = ch1+ch2;
				ch1Sum+=ch1;
				ch2Sum+=ch2;
				ch3Sum+=ch3;
			}
		}
	}

	def double ch1Mean = ch1Sum/N;
	def double ch2Mean = ch2Sum/N;
	def double ch3Mean = ch3Sum/N;*/

        double ch1mch1MeanSqSum = 0, ch2mch2MeanSqSum = 0, ch3mch3MeanSqSum = 0;

        double sumX=0, sumXX=0, sumXY=0, sumYY=0, sumY=0;

        //calulate variances
        for (int s=1; s<nSlices; s++) {
            ImageProcessor ip1 = img1.getStack().getProcessor(s);
            ImageProcessor ip2 = img2.getStack().getProcessor(s);
            for (int y=0; y<rheight; y++) {
                for (int x=0; x<rwidth; x++) {
                    ch1 = ip1.getPixelValue(x,y);
                    ch2 = ip2.getPixelValue(x,y);
                    ch3 = ch1+ch2;
                    ch1mch1MeanSqSum+= (ch1-ch1Mean)*(ch1-ch1Mean);
                    ch2mch2MeanSqSum+= (ch2-ch2Mean)*(ch2-ch2Mean);
                    ch3mch3MeanSqSum+= (ch3-ch3Mean)*(ch3-ch3Mean);
                    //calc pearsons for original image
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
