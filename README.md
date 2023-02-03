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

* Detect vessels with a Median filtering + DoG filtering + thresholding + Median filtering
* Detect cells with Cellpose
* Compute each cell distance to nearest vessel
Keep Th cells colocalizing with a nucleus only
* Measure ORF1p intensity in the nucleus and the cytoplasm of each Th cell
* Measure ORF1p intensity in Th-negative nuclei


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto* and model

### Version history

Version 2 released on February 3, 2023.
