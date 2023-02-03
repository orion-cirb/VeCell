# Sox_10

* **Developed for:** Anne-Cécile
* **Team:** Cohen-Salmon
* **Date:** June 2022
* **Software:** Fiji


### Images description

3D images taken with a x2.3 objective on an AxioZoom

3 channels:
  1. *EGFP:* Vessels
  2. *DsRed:* Sox9/Sox10 cells
  3. *Cy5:* Vessels
  
A *.roi* or *.zip* file containing ROI(s) must be provided with each image.

### Plugin description

In each ROI:
* Detect cells with Stardist or DoG filtering + thresholding
* Compute the distance between each cell and its neighbors
* Compute the G-function related spatial distribution index of the population of cells
* Detect vessels with a Median filtering + DoG filtering + thresholding
* Compute the vessels 2D skeleton, 3D distance map and 3D inverse distance map
* Get the distance between each cell and its nearest vessel using the inverse distance map
* Compute the radius of the corresponding vessels using the vessels skeleton and distance map
* Give cells volume and intensity

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Stardist** conda environment + *StandardFluo.zip* model

### Version history

Version 1 released on June 17, 2022.
