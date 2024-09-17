setBatchMode(true);

rename("input");
run("CLIJ2 Macro Extensions", "cl_device=[NVIDIA RTX A5000]");

// 1. MEDIAN FILTER
image14 = "input";
Ext.CLIJ2_push(image14);
image15 = "median1";
radius_x = 4.0;
radius_y = 4.0;
radius_z = 1.0;
Ext.CLIJ2_median3DSphere(image14, image15, radius_x, radius_y, radius_z);
Ext.CLIJ2_pull(image15);

// 2. DOG FILTER 1
// 2.a. Difference of gaussian
Ext.CLIJ2_push(image15);
image16 = "difference_of_gaussian1";
sigma1x = 4.0;
sigma1y = 4.0;
sigma2x = 8.0;
sigma2y = 8.0;
Ext.CLIJ2_differenceOfGaussian2D(image15, image16, sigma1x, sigma1y, sigma2x, sigma2y);
Ext.CLIJ2_pull(image16);

// 2.b. Threshold
Ext.CLIJ2_push(image16);
image17 = "threshold1";
Ext.CLIJ2_thresholdTriangle(image16, image17);
Ext.CLIJ2_pull(image17);

// 3. DOG FILTER 2 (OPTIONAL)
if(false) { // Set the value to true if you want to apply it, or false if you don't.
	// 3.a. Difference of gaussian
	Ext.CLIJ2_push(image15);
	image18 = "difference_of_gaussian2";
	sigma1x = 7.0;
	sigma1y = 7.0;
	sigma2x = 14.0;
	sigma2y = 14.0;
	Ext.CLIJ2_differenceOfGaussian2D(image15, image18, sigma1x, sigma1y, sigma2x, sigma2y);
	Ext.CLIJ2_pull(image18);
	
	// 3.b. Threshold
	Ext.CLIJ2_push(image18);
	image19 = "threshold2";
	Ext.CLIJ2_thresholdTriangle(image18, image19);
	Ext.CLIJ2_pull(image19);
	
	// 3.c. Max of the two stacks
	imageCalculator("Max stack", "threshold1","threshold2");
}

// 4. CLOSING FILTER
// 4.a. Maximum
Ext.CLIJ2_push(image17);
image20 = "maximum";
radius_x = 8.0;
radius_y = 8.0;
radius_z = 1.0;
Ext.CLIJ2_maximum3DSphere(image17, image20, radius_x, radius_y, radius_z);
Ext.CLIJ2_pull(image20);

// 4.b. Minimum
Ext.CLIJ2_push(image20);
image21 = "minimum";
radius_x = 8.0;
radius_y = 8.0;
radius_z = 1.0;
Ext.CLIJ2_minimum3DSphere(image20, image21, radius_x, radius_y, radius_z);
Ext.CLIJ2_pull(image21);

// 5. MEDIAN FILTER
Ext.CLIJ2_push(image21);
image22 = "result";
radius_x = 1.0;
radius_y = 1.0;
radius_z = 1.0;
Ext.CLIJ2_median3DSphere(image21, image22, radius_x, radius_y, radius_z);
Ext.CLIJ2_pull(image22);

close("median1");
close("difference_of_gaussian1");
close("threshold1");
close("difference_of_gaussian2");
close("threshold2");
close("sum");
close("maximum");
close("minimum");
Ext.CLIJ2_clear();

setBatchMode(false);
