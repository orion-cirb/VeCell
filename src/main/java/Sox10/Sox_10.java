package Sox10;




import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.measure.Calibration;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.RoiScaler;
import java.awt.Rectangle;
import java.io.File;
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
            
            // write headers for results tables
            sox.writeHeaders(outDirResults);
                        
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                // Test if ROI file exist
                String roi_file  = (new File(imageDir+rootName+".zip").exists()) ? imageDir+rootName+".zip" :  ((new File(imageDir+rootName+".roi").exists()) ? imageDir+rootName+".roi" : null);
                if (roi_file == null) {
                    IJ.showStatus("No ROI file found !") ;
                    IJ.log("No ROIs found ! bye bye \n");
                    return;
                }
                List<Roi> rois = new ArrayList<>();
                    // find rois
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roi_file);
                rois = Arrays.asList(rm.getRoisAsArray());
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                // open Cell channel
                // case Zeiss 
                // take first serie in image pyramidal
                int nSeries = reader.getSeriesCount();
                int series = (fileExt.equals("czi")) ? 0 : nSeries -1;
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setQuiet(true);
                options.setCrop(true);
                
                // For each roi open cropped image
                for (Roi roi : rois) {
                    nucIndex++;
                    Rectangle rectRoi = roi.getBounds();
                    options.setCropRegion(series, new Region(rectRoi.x, rectRoi.y, rectRoi.width, rectRoi.height));
                    options.setCBegin(series, channelIndex[2]);
                    options.setCEnd(series, channelIndex[2]); 
                    // cells image
                    ImagePlus imgCells = BF.openImagePlus(options)[0];
                    Objects3DPopulation cellPop = new Objects3DPopulation();
                    if (sox.cellsDetection.equals("DOG"))
                        cellPop = sox.findCellsDoG(imgCells, roi);
                    else
                        cellPop = sox.stardistNucleiPop(imgCells, roi);
                    System.out.println(cellPop.getNbObjects()+" cells found");
                    
                    // Vessel channel
                    Objects3DPopulation vesselPop = new Objects3DPopulation();
                    ArrayList<Double> dist = new ArrayList<>();
                    ArrayList<Double> diam = new ArrayList<>();
                    // vessel image
                    if (sox.vessel){
                        options.setCBegin(series, channelIndex[0]);
                        options.setCEnd(series, channelIndex[0]); 
                        ImagePlus imgVessel1 = BF.openImagePlus(options)[0];
                        options.setCBegin(series, channelIndex[1]);
                        options.setCEnd(series, channelIndex[1]);
                        ImagePlus imgVessel2 = BF.openImagePlus(options)[0];
                        ImagePlus imgVessel =  new ImageCalculator().run("add stack create", imgVessel1, imgVessel2);
                        sox.closeImages(imgVessel1);
                        sox.closeImages(imgVessel2);
                        ImagePlus imgVesselTube = sox.tubeness(imgVessel);
                        vesselPop = sox.findVessel(imgVesselTube, roi);
                        dist = sox.findCellVesselDist(cellPop, vesselPop);
                        ImagePlus imgVesselMap = sox.localThickness3D(imgVesselTube);
                        diam = sox.findVesselDiameter(cellPop, vesselPop, imgVesselMap);
                        sox.closeImages(imgVesselTube);
                        sox.closeImages(imgVesselMap);
                    }
                    
                    sox.saveCellsImage(cellPop, vesselPop, imgCells, dist, outDirResults+rootName+"_"+roi.getName()+".tif");
                                        
                    // find parameters
                    sox.computeNucParameters(cellPop, vesselPop, dist, diam, imgCells, roi.getName(), roi, rootName, outDirResults);
                    sox.closeImages(imgCells);
                }
            }
            sox.closeResults();
            IJ.showStatus("Processing done....");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Sox_10.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

