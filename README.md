# VeCell

* **Developed for:** Naomie
* **Team:** Cohen-Salmon
* **Date:** September 2025
* **Software:** Fiji


### Images description

3D images taken on an Axio Zoom microscope with Apotome module:
* Zoom: 160x 
* File format: .czi 
* Voxel size: 0.4063 µm (XY), 4.00 µm (Z)

2 channels:
  1. *DsRed:* Astrocytes (mandatory)
  2. *Cy5:* Vessels (optional)

For each image, a corresponding `imageName.zip` file containing several ROIs must be provided, as analysis is performed within each ROI. If no ROI file is supplied, the image is not analyzed. ROIs should be named `position layerName` (e.g. `0095-0429 l2`). If ROIs are drawn on a lower-resolution image, scaling factor must be specified in the dialog box.

### Plugin description

In each ROI:
* Detect astrocytes with Cellpose
* Compute distance between each cell and its nearest neighbors
* Compute G-function related spatial distribution index of the population of cells
* If vessels channel provided:
  * Detect vessels with Quantile Based Normalization + median filtering + DoG filtering + thresholding + closing filtering + median filtering + fill holes
  * Compute vessels skeleton and provide vessels volume, diameter, length, branches number, and junctions number
  * Compute distance between each cell and its nearest vessel

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + fine-tuned model `cyto2_sox9_p5-15-60_27-11-24`
  
The Fiji plugins, Cellpose environment, and associated model are installed on the workstations in the ORION image analysis room.

The dataset used for Cellpose model fine-tuning is stored in the ORION storage space on the ISIS_PROD_NAS server.

### Version history

Version 5.1 released on September 22, 2025.

Improvements compared to version 4:
* Plugin renamed
* Code cleaned
* Dialog box changed
* Vessels segmentation improved:
  * Channel quantile-based normalization
  * 2 DoG filters can be applied if thin and thick vessels appear in the image
  * Fill 2D holes with areas below specified max value
* Vessels skeleton branches filtered out by length
* Vessels skeleton length, branches number, and junctions number provided
