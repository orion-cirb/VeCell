import VeCell_Tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


/**
 * @author ORION-CIRB
 */   
public class VeCell implements PlugIn {
    
    private Tools tools = new Tools();
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            }
            
            String imageDir = IJ.getDirectory("Choose images folder");
            if (imageDir == null) {
                return;
            }
            
            // Find extension of first image in input folder
            String fileExt = tools.findImageType(new File(imageDir));
            // Find all images with corresponding extension in folder
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("ERROR", "No image found with " + fileExt + " extension in " + imageDir + " folder");
                return;
            }
            
            // Instantiate metadata and reader
            IMetadata meta = MetadataTools.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);

            // Find channels name
            String[] chMeta = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] chOrder = tools.dialog(chMeta);
            if (chOrder == null) {
                return;
            } else if(chOrder[0].equals("None")) {
                IJ.showMessage("ERROR", "Astrocytes channel not defined");
                return;
            }
            
            // Create output folder
            String outDir = imageDir + File.separator + "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write header in results file
            tools.writeHeaders(outDir, !chOrder[1].equals("None"));
            
            IJ.setForegroundColor(255, 255, 255);
            IJ.setBackgroundColor(0, 0, 0);
                 
            for (String f: imageFiles) {
                reader.setId(f);
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ---");
                
                tools.print("- Loading ROIs -");
                List<Roi> rawRois = tools.loadRois(imageDir, rootName);
                if (rawRois == null) continue;
                
                List<Roi> rois = tools.scaleRois(rawRois, tools.roiScaling);
                Region bBox = tools.getBoundingBox(rois);
                tools.translateRois(rois, bBox);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setCrop(true);
                options.setCropRegion(0, bBox);
                options.doCrop();
                
                 // Astrocytes channel
                tools.print("- Opening astrocytes channel -");
                int chIndex = ArrayUtils.indexOf(chMeta, chOrder[0]);                
                options.setCBegin(tools.imgSeries, chIndex);
                options.setCEnd(tools.imgSeries, chIndex);
                ImagePlus imgCell = BF.openImagePlus(options)[0];
                
                tools.print("- Detecting astrocytes -");
                ImagePlus imgCellMask = tools.cellposeDetection(imgCell);
                
                // Vessels channel
                ImagePlus imgVessel = null, imgVesselMask = null, imgVesselSkel = null;
                ImageFloat imgVesselDistMap = null, imgVesselDistMapInv = null;
                if (!chOrder[1].equals("None")){
                    tools.print("- Opening vessels channels -");
                    chIndex = ArrayUtils.indexOf(chMeta, chOrder[1]);
                    options.setCBegin(tools.imgSeries, chIndex);
                    options.setCEnd(tools.imgSeries, chIndex);
                    imgVessel = BF.openImagePlus(options)[0];
                    if (!chOrder[2].equals("None")) {
                        chIndex = ArrayUtils.indexOf(chMeta, chOrder[2]);
                        options.setCBegin(tools.imgSeries, chIndex);
                        options.setCEnd(tools.imgSeries, chIndex);
                        ImagePlus imgVessel2 = BF.openImagePlus(options)[0];
                        // Add two vessels channels
                        new ImageCalculator().run("Add stack", imgVessel, imgVessel2);
                        tools.closeImage(imgVessel2);
                    }
                    
                    tools.print("- Detecting vessels -");
                    imgVesselMask = tools.vesselsSegmentation(imgVessel);
                    
                    tools.print("- Computing vessels distance maps -");
                    imgVesselDistMap = tools.distanceMap3D(imgVesselMask, false);
                    imgVesselDistMapInv = tools.distanceMap3D(imgVesselMask, true);
                    
                    tools.print("- Computing vessels skeleton -");
                    imgVesselSkel = tools.skeletonize3D(imgVesselMask);
                    // Prune vessels skeleton small branches
                    imgVesselSkel = tools.pruneSkeleton(imgVesselSkel);
                }
                
                // For each roi, open cropped image
                ImageHandler imhCell = ImageHandler.wrap(imgCell).createSameDimensions();
                ImageHandler imhCellDist = imhCell.createSameDimensions();
                ImageHandler imhVessel = imhCell.createSameDimensions();
                ImageHandler imhSkel = imhCell.createSameDimensions();
                for (Roi roi : rois) {
                    tools.print("- Computing parameters and saving results for ROI " + roi.getName() + " -");

                    // Cells channel
                    Objects3DIntPopulation popCellRoi = tools.getPopInRoi(imgCellMask, imgCell, roi);
                    
                    // Vessels channel
                    double vesselVolRoi = 0;
                    Object3DInt objVesselRoiDil = null, objSkelRoi = null, objSkelRoiDil = null;
                    if (imgVessel != null){
                        vesselVolRoi = new MeasureVolume(tools.getObjInRoi(imgVesselMask, roi, 0)).getVolumeUnit();
                        objVesselRoiDil = tools.getObjInRoi(imgVesselMask, roi, tools.roiDilation);
                        objSkelRoi = tools.getObjInRoi(imgVesselSkel, roi, 0);
                        objSkelRoiDil = tools.getObjInRoi(imgVesselSkel, roi, tools.roiDilation);

                        tools.print("Computing distance between each cell and its closest vessel...");
                        tools.computeCellsVesselsDists(popCellRoi, imgVesselDistMapInv);
                    }
                    
                    tools.print("Drawing results...");
                    tools.drawResultsInRoi(popCellRoi, objVesselRoiDil, objSkelRoi, imhCell, imhCellDist, imhVessel, imhSkel);   
                    
                    tools.print("Writing results...");
                    tools.writeResultsInRoi(popCellRoi, objSkelRoiDil, objSkelRoi, imgVesselDistMap, imgCell, roi, vesselVolRoi, outDir, rootName);
                }
                
                tools.saveCloseDrawings(imgCell, imgVessel, imhCell, imhCellDist, imhVessel, imhSkel, outDir, rootName);
                
                tools.closeImage(imgCell);
                tools.closeImage(imgCellMask);
                if (imgVessel != null){
                    tools.closeImage(imgVessel);
                    tools.closeImage(imgVesselMask);
                    tools.closeImage(imgVesselDistMap.getImagePlus());
                    tools.closeImage(imgVesselDistMapInv.getImagePlus());
                    tools.closeImage(imgVesselSkel);
                }
            }
            tools.closeResults();
            tools.print("--- All done! ---");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(VeCell.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

