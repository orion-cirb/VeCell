# Sox_10

* **Developed for:** Katia
* **Team:** Cohen-Salmon
* **Date:** February 2023
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
* Compute the distance between each cell and its neighboring cells
* Compute the G-function related spatial distribution index of the population of cells
* Detect vessels with a Median filtering + DoG filtering + thresholding + Median filtering
* Compute the distance between each cell and its nearest vessel and compute the radius of the corresponding vessel
* Give cells volume and intensity

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto* and model

### Version history

Version 2 released on February 3, 2023.
