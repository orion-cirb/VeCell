# Sox_10

* **Developed for:** Naomie
* **Team:** Cohen-Salmon
* **Date:** April 2023
* **Software:** Fiji


### Images description

3D images taken with a 1.6x objective on an AxioZoom microscope.

3 channels:
  1. *EGFP:* Vessels
  2. *DsRed:* Sox9 cells
  3. *Cy5:* Vessels (optional)
  
A *.roi* or *.zip* file containing ROI(s) must be provided with each image.

### Plugin description

In each ROI:
* Detect cells with Cellpose
* Compute the distance between each cell and its neighbors
* Compute the G-function related spatial distribution index of the population of cells
* Detect vessels with median filtering + DoG filtering + thresholding + median filtering
* Compute the vessels 3D skeleton, 3D distance map and 3D inverse distance map
* Get the distance between each cell and its nearest vessel using the inverse distance map
* Compute the radius of the corresponding vessels using the vessels skeleton and distance map
* Give vessels volume and cells volume and intensity

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto* model

### Version history

Version 4 released on October 31, 2023.

Small improvements compared to version 3:
* Adapted cells and vessels detection to smaller resolution of the images (1.6x objective instead of 2.3x)
* Added vessels total volume in each ROI in globalResults.csv file
* Results for each ROI drawn in same image
