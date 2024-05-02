[![Build Status](https://github.com/BIOP/ijp-jacop-b/actions/workflows/build.yml/badge.svg)](https://github.com/BIOP/ijp-jacop-b/actions/workflows/build.yml)

# ijp-jacop-b
[Jacop by F. Cordelières](https://imagej.net/plugins/jacop) revamped by the BIOP

An update to the JACoP Plugin that helps in the management of ROIs, Z sections and helps generate cleaner reports

Based on [JACoP](https://imagej.net/plugins/jacop) and [Coloc2](https://imagej.net/plugins/coloc-2)

# Installation

Please use our PTBIOP update site

## Available methods/metrics

* [Pearson's correlation Coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
* [Mander's coefficient Article](https://imagej.net/_images/2/24/Manders.pdf)
* [Costes Article](https://imagej.net/_images/e/e2/Costes_etalColoc.pdf)
* [Li Article](https://imagej.net/_images/9/9e/LietAlColoc.pdf)

## Input
To run BIOP JACoP you will need
- A multi-channel image ( no need to split into independent channels )
- To define channels to use 
- To define Thresholds, either using a User-selected manually fixed value (left) or an [Automatic Thresholding Method](https://imagej.net/Auto_Threshold) (right)

![manual or auto threshold](https://github.com/BIOP/ijp-jacop-b/assets/20223054/26fedd0a-3230-40b4-ba50-30d21f98ba86)

We recommend to use the Automatic threshold methods. They will compensate for subtle changes in illumination (microscopes are not perfect), or (in case of variable expression of fluorescent proteins) the intensity from one cell to another (you would need to define a ROI for each cell first). 
Nevertheless, it might fail if your images are very heterogenous (number of cells, intensities range, areas of signal...)
If you can't find any auto-threshold method that gives satisfying results, then you should consider using a manually defined fixed value (based on controls).

(IMPORTANT) We urge you to have controls (mono-stained samples, acquired the same way as the test ones, same number of channels, etc...) to verify that the defined thresholds are above the signal one can observe in an unstained sample. You will always have some crosstalk/bleedthrough, always!
 
## Outputs

- Results Table
- Output image with thresholded mask, fluorogram (optional), example randomized image.

![jacop output](https://github.com/BIOP/ijp-jacop-b/assets/20223054/a7ca4c93-18c2-42ee-9ba1-fb23d3037e3f)

## Some functionalities

### Region of Interest (ROIs)

With many ROIs in the ROI manager (defining cells for example), you will get an analysis per ROI.

- **Crop ROIs** , generates a cropped image for each ROI.

### Z-stack
-  **Consider Z slices Separately** , output a single Z-stack image **BUT** the analysis is performed on each individual slice.
We recommend always starting your pilot experiment with z-stacks and assess (using this option) if the result depends on doing the analysis on 2D or 3D image and continue your acquisition campaign accordingly.

-  **Set Auto Thresholds On Stack Histogram** , when using one of the [automatic threshold methods](https://imagej.net/Auto_Threshold) you can apply the auto-threshold either on individual slices or on the entire stack.
> [!WARNING]  
> This checkbox is only valid if **Consider Z slices Separately** is checked????

![histogram or not](https://github.com/BIOP/ijp-jacop-b/assets/20223054/0bfb815a-8284-4759-9c00-64b13cd84039)

### Advanced parameters

You can fine tune the histogram size and the costes plot. This is useful essentially in order to keep the same bounds for different images or for the all slices of a single image.

![advanced parameters](https://github.com/BIOP/ijp-jacop-b/assets/20223054/008ac0a7-f843-463f-bcd6-e97517b34a15)




