# Sox_10

* **Developed for:** Naomie
* **Team:** Cohen-Salmon
* **Date:** April 2023
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
* Detect cells with Cellpose
* Compute the distance between each cell and its neighbors
* Compute the G-function related spatial distribution index of the population of cells
* Detect vessels with:
  * Median filtering + DoG filtering + thresholding + Median filtering
  * Background subtraction + tubeness filtering + thresholding + Median filtering 
* Compute the vessels 3D skeleton, 3D distance map and 3D inverse distance map
* Get the distance between each cell and its nearest vessel using the inverse distance map
* Compute the radius of the corresponding vessels using the vessels skeleton and distance map
* Give cells volume and intensity

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto* model

### Version history

Version 3 released on April 26, 2023.

Small improvements compared to version 2:
* Compute the difference between the area of the cells G-function curve and the area of the random organisation G-function curve
* Better segment vessels (the user can choose between DoG or tubeness filtering)
