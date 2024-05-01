//Delete log, delete results and close all windows
print("\\Clear");
run("Clear Results");
close("*");
close("Results");
close("Log");

//Load data and crop it

run("Neuron (5 channels)");
makeRectangle(139, 166, 128, 128);
run("Duplicate...", "title=Rat_Hippocampal_Neuron-2C.tif duplicate channels=1-2");
selectImage("Rat_Hippocampal_Neuron.tif");
// See https://forum.image.sc/t/inconsistent-costes-threshold-results-across-different-plugins/94558

close();

//Run BIOP-JACOP on 16-bit data and repeat flipping channels

selectImage("Rat_Hippocampal_Neuron-2C.tif");
run("BIOP JACoP", "channel_a=1 channel_b=2 threshold_for_channel_a=[Costes Auto-Threshold] threshold_for_channel_b=[Costes Auto-Threshold] manual_threshold_a=954 manual_threshold_b=846 get_manders costes_block_size=5 costes_number_of_shuffling=100");

selectImage("Rat_Hippocampal_Neuron-2C.tif");
run("BIOP JACoP", "channel_a=2 channel_b=1 threshold_for_channel_a=[Costes Auto-Threshold] threshold_for_channel_b=[Costes Auto-Threshold] manual_threshold_a=954 manual_threshold_b=846 get_manders costes_block_size=5 costes_number_of_shuffling=100");

//Make an 8-bit version of the data

selectImage("Rat_Hippocampal_Neuron-2C.tif");
run("Duplicate...", "title=Rat_Hippocampal_Neuron-2C-8bit.tif duplicate");
setOption("ScaleConversions", true);
run("8-bit");

//Run BIOP-JACOP on 8-bit data and repeat flipping channels
/*
selectImage("Rat_Hippocampal_Neuron-2C-8bit.tif");
run("BIOP JACoP", "channel_a=1 channel_b=2 threshold_for_channel_a=[Costes Auto-Threshold] threshold_for_channel_b=[Costes Auto-Threshold] manual_threshold_a=954 manual_threshold_b=846 get_manders costes_block_size=5 costes_number_of_shuffling=100");

selectImage("Rat_Hippocampal_Neuron-2C-8bit.tif");
run("BIOP JACoP", "channel_a=2 channel_b=1 threshold_for_channel_a=[Costes Auto-Threshold] threshold_for_channel_b=[Costes Auto-Threshold] manual_threshold_a=954 manual_threshold_b=846 get_manders costes_block_size=5 costes_number_of_shuffling=100");
*/
//Split stacks into 1 image per channel for other colocalization plugins

/*
selectImage("Rat_Hippocampal_Neuron-2C.tif");
run("Split Channels");

selectImage("Rat_Hippocampal_Neuron-2C-8bit.tif");
run("Split Channels");

//Run Coloc2 on 16-bit and 8-bit data flipping order of the channels in each case

run("Coloc 2", "channel_1=C1-Rat_Hippocampal_Neuron-2C.tif channel_2=C2-Rat_Hippocampal_Neuron-2C.tif roi_or_mask=<None> threshold_regression=Costes manders'_correlation psf=3 costes_randomisations=10");

run("Coloc 2", "channel_1=C2-Rat_Hippocampal_Neuron-2C.tif channel_2=C1-Rat_Hippocampal_Neuron-2C.tif roi_or_mask=<None> threshold_regression=Costes manders'_correlation psf=3 costes_randomisations=10");

run("Coloc 2", "channel_1=C1-Rat_Hippocampal_Neuron-2C-8bit.tif channel_2=C2-Rat_Hippocampal_Neuron-2C-8bit.tif roi_or_mask=<None> threshold_regression=Costes manders'_correlation psf=3 costes_randomisations=10");

run("Coloc 2", "channel_1=C2-Rat_Hippocampal_Neuron-2C-8bit.tif channel_2=C1-Rat_Hippocampal_Neuron-2C-8bit.tif roi_or_mask=<None> threshold_regression=Costes manders'_correlation psf=3 costes_randomisations=10");

//Run Colocalization threshold on 16-bit and 8-bit data flipping order of the channels in each case

run("Colocalization Threshold", "channel_1=C1-Rat_Hippocampal_Neuron-2C.tif channel_2=C2-Rat_Hippocampal_Neuron-2C.tif use=None channel=[Red : Green] show_0 set show show_0 mander's_0");

run("Colocalization Threshold", "channel_1=C2-Rat_Hippocampal_Neuron-2C.tif channel_2=C1-Rat_Hippocampal_Neuron-2C.tif use=None channel=[Red : Green] show_0");

run("Colocalization Threshold", "channel_1=C1-Rat_Hippocampal_Neuron-2C-8bit.tif channel_2=C2-Rat_Hippocampal_Neuron-2C-8bit.tif use=None channel=[Red : Green] show_0");

run("Colocalization Threshold", "channel_1=C2-Rat_Hippocampal_Neuron-2C-8bit.tif channel_2=C1-Rat_Hippocampal_Neuron-2C-8bit.tif use=None channel=[Red : Green] show_0");

*/


