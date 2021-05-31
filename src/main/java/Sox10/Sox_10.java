package Sox10;




import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.frame.RoiManager;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import org.apache.commons.io.FilenameUtils;


/**
 *
 * @author phm
 */
        
        
public class Sox_10 implements PlugIn {
    
    private boolean canceled;
    private String imageDir;
    private static String outDirResults;
    private BufferedWriter outPutResults;
    public static Calibration cal = new Calibration();
    // min and max cells filter volume
    private final double minCells = 10;
    private final double maxCells = 5000;
    
    private Sox_10_Tools sox = new Sox_10_Tools();
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            if (!sox.checkInstalledModules()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Images folder");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String fileExt = sox.findImageType(inDir);
            ArrayList<String> imageFiles = sox.findImages(imageDir, fileExt);
            if (imageFiles == null) {
                return;
            }
            
            // create output folder
            outDirResults = imageDir + "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }

            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            int nucIndex = 0;
            // Find channel names , calibration
            reader.setId(imageFiles.get(0));
            cal = sox.findImageCalib(meta);
            String[] chsName = sox.findChannels(imageFiles.get(0), meta, reader);
            
            int[] channelIndex = sox.dialog(chsName);
            cal = sox.getCalib();
            if (channelIndex == null)
                return;
            
            // Write headers results for results file{
            FileWriter fileResults = new FileWriter(outDirResults +"Results.xls", false);
            outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("Image name\tRoi name\tRoi volume\tNb Cell\tCell mean intensity\tCell sd intensity\tCell mean volume\t Cell sd volume\ttotalVolume\tAverageDistanceToClosestNeighbor\tSDDistanceToClosestNeighbor\tDistributionAleatoireStat \n");
            outPutResults.flush();
                    
                
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                int series = reader.getSeries();
                int width = meta.getPixelsSizeX(series).getValue();
                int height = meta.getPixelsSizeY(series).getValue();
                reader.setSeries(series);
                boolean roiFile = true;
                // Test if ROI file exist
                String roi_file  = (new File(imageDir+rootName+".zip").exists()) ? imageDir+rootName+".zip" :  ((new File(imageDir+rootName+".roi").exists()) ? imageDir+rootName+".roi" : null);
                if (roi_file == null) {
                    IJ.showStatus("No ROI file found !") ;
                    roiFile = false;
                }
                List<Roi> rois = new ArrayList<>();
                if (roiFile) {
                    // find rois
                    RoiManager rm = new RoiManager(false);
                    rm.runCommand("Open", roi_file);
                    rois = Arrays.asList(rm.getRoisAsArray());
                }
                // define roi as all image
                else {
                    //rois.add(new Roi(0, 0, width, height));
                    IJ.log("No ROIs found ! byebye \n");
                    return;
                }
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setStitchTiles(true);
                options.setOpenAllSeries(true);
               
                // open Cell channel, for all series
                int nseries = reader.getSeriesCount();
                for (int s=0; s<nseries;s++) {
                options.setCBegin(s, channelIndex[1]);
                options.setCEnd(s, channelIndex[1]);      
                }
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setQuiet(true);
                ImagePlus[] allImages = BF.openImagePlus(options);
                // get the stitched image
                ImagePlus wholeImage = allImages[allImages.length-1];
                // close the other images
                for (int k=0; k<(allImages.length-1); k++){
                    sox.closeImages(allImages[k]);
                }
                
                // For each roi open cropped image
                for (Roi roi : rois) {
                    nucIndex++;
                    
                    Rectangle rectRoi = roi.getBounds();
                    wholeImage.setRoi(rectRoi);
                    ImagePlus imgCells = new Duplicator().run(wholeImage);
                    Objects3DPopulation cellPop = sox.findCellsDoG(imgCells);
                    System.out.println(cellPop.getNbObjects()+" cells found");
                    // find parameters
                    sox.computeNucParameters(cellPop, imgCells, roi.getName(), rootName, outDirResults, outPutResults);
                    sox.saveCellsImage(cellPop, imgCells, outDirResults+rootName+"_"+roi.getName()+".tif");
                    sox.closeImages(imgCells);
                }
                sox.closeImages(wholeImage);
            }
            outPutResults.close();
            IJ.showStatus("Processing done....");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Sox_10.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

