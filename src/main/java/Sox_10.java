import Sox10_Tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;


/**
 * @author ORION-CIRB
 */   
public class Sox_10 implements PlugIn {
    
    private Tools tools = new Tools();
    private String imageDir;
    private static String outDirResults;
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            }
            
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with fileExt extension
            String fileExt = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header in results file
            tools.writeHeaders(outDirResults);
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);

            // Find channels name
            String[] chsName = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            int[] channelIndex = tools.dialog(chsName);
            if (channelIndex == null) {
                IJ.showMessage("Error", "Plugin canceled");
                return;
            }
                 
            for (String f : imageFiles) {
                reader.setId(f);
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                
                // Find ROI(s)
                String roiFile = imageDir+rootName+".roi";
                if (!new File(roiFile).exists())
                    roiFile = imageDir+rootName+".zip";
                if (!new File(roiFile).exists()) {
                    tools.print("ERROR: No ROI file found for image " + rootName);
                    IJ.showMessage("Error", "No ROI file found for image " + rootName);
                    continue;
                }
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roiFile);
                List<Roi> rois = Arrays.asList(rm.getRoisAsArray());
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setCrop(true);
                
                // Compute pyramidal factor
                int series = 0;
                int pyramidalFactor = tools.getPyramidalFactor(reader);
                
                // Vessels channel
                ImagePlus imgVessel = null;
                ImagePlus vesselsDetection = null;
                ImagePlus vesselsDistMap = null;
                ImagePlus vesselsDistMapInv = null;
                ImagePlus vesselsSkel = null;
                if (tools.vessel){
                    tools.print("Opening vessels channels...");
                    options.setCBegin(series, channelIndex[0]);
                    options.setCEnd(series, channelIndex[0]); 
                    ImagePlus imgVessel1 = BF.openImagePlus(options)[0];
                    
                    options.setCBegin(series, channelIndex[1]);
                    options.setCEnd(series, channelIndex[1]);
                    ImagePlus imgVessel2 = BF.openImagePlus(options)[0];
                    // Add two vessels chanels
                    imgVessel =  new ImageCalculator().run("add stack create", imgVessel1, imgVessel2);
                    tools.closeImage(imgVessel1);
                    tools.closeImage(imgVessel2);
                    
                    tools.print("Detecting vessels...");
                    vesselsDetection = tools.vesselsDetection(imgVessel, pyramidalFactor);
                    tools.print("Computing vessels skeleton...");
                    vesselsSkel = tools.vesselsSkeletonization(vesselsDetection);
                    tools.print("Computing vessels distance maps...");
                    vesselsDistMap = tools.localThickness3D(vesselsDetection, false);
                    vesselsDistMapInv = tools.localThickness3D(vesselsDetection, true);
                }
                 // Cells channel
                tools.print("Opening cells channel...");
                options.setCBegin(series, channelIndex[2]);
                options.setCEnd(series, channelIndex[2]);
                ImagePlus imgCellsOrg = BF.openImagePlus(options)[0];
                // if vessels exist remove vessel in image cells 
                ImagePlus imgCells = (tools.vessel) ?  new ImageCalculator().run("subtract stack create", imgCellsOrg, imgVessel) : imgCellsOrg;
                tools.closeImage(imgCellsOrg);
                tools.print("Detecting cells...");
                ImagePlus cellsDetection = tools.cellposeDetection(imgCells);
                
                // For each roi, open cropped image
                for (Roi roi : rois) {
                    String roiName = roi.getName();
                    tools.print("- Analyzing ROI " + roiName + " -");
                    Roi scaledRoi = tools.scaleRoi(roi, pyramidalFactor);

                    // Cells channel
                    Objects3DIntPopulation cellPop = tools.getCellsInRoi(cellsDetection, imgCells, scaledRoi);
                    System.out.println(cellPop.getNbObjects() + " cells found in ROI "+scaledRoi.getName());
                    tools.closeImage(cellsDetection);
                    
                    // Vessels channel
                    Objects3DIntPopulation vesselPop = new Objects3DIntPopulation();
                    ArrayList<Double> dist = new ArrayList<>();
                    ArrayList<Double> radius = new ArrayList<>();
                    if (tools.vessel){
                        vesselPop = tools.getVesselsInRoi(vesselsDetection, scaledRoi);
                        Object3DInt vesselsSkelObj = tools.getVesselsSkelInRoi(vesselsSkel, scaledRoi);
                        
                        tools.print("Computing cells distance to vessels...");
                        dist = tools.findCellVesselDist(cellPop, vesselsDistMapInv);
                        
                        tools.print("Computing vessels radius...");
                        radius = tools.findVesselRadius(cellPop, vesselsDistMap, vesselsSkelObj);
                    }
                    
                    tools.print("Drawing results...");
                    tools.drawResults(cellPop, vesselPop, imgCells, imgVessel, dist, outDirResults+rootName+"_"+roi.getName());        
                    tools.print("Writing results...");
                    tools.writeResults(cellPop, vesselPop, dist, radius, imgCells, roi.getName(), scaledRoi, rootName, outDirResults);
                    tools.closeImage(imgVessel);
                }
                
                tools.closeImage(imgCells);
                
                if (tools.vessel){
                    tools.closeImage(imgVessel);
                    tools.closeImage(vesselsDetection);
                    tools.closeImage(vesselsDistMap);
                    tools.closeImage(vesselsDistMapInv);
                    tools.closeImage(vesselsSkel);
                }
            }
            
            tools.closeResults();
            tools.print("--- All done! ---");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Sox_10.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

