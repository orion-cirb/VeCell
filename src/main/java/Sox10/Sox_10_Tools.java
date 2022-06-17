package Sox10;


import Sox10_StardistOrion.StarDist2D;
import features.TubenessProcessor;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.CalibrationBar;
import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.RoiScaler;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;
import java.awt.Rectangle;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import static mcib3d.geom.Object3D_IJUtils.createObject3DVoxels;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;
import mcib3d.utils.ThreadUtil;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;
import sc.fiji.localThickness.Clean_Up_Local_Thickness;
import sc.fiji.localThickness.EDT_S1D;
import sc.fiji.localThickness.Local_Thickness_Parallel;
import mcib_plugins.analysis.SpatialAnalysis;       
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class Sox_10_Tools {
   
    // min volume in µm^3 for cells
    public double minCell = 20;
    // max volume in µm^3 for cells
    public double maxCell = 1000;
    // sigmas for DoG
    private double sigma1 = 20;
    private double sigma2 = 40;
    // Dog parameters
    private Calibration cal = new Calibration(); 
    public boolean vessel = false;
    private boolean doF = false;
    private double radiusNei = 50; //neighboring radius
    private int nbNei = 10; // K neighborg
        
    private BufferedWriter outPutResults;
    private BufferedWriter outPutDistances;
    
    // Stardist
    private Object syncObject = new Object();
    private final double stardistPercentileBottom = 0.2;
    private final double stardistPercentileTop = 99.8;
    private final double stardistProbThreshNuc = 0.80;
    private final double stardistOverlayThreshNuc = 0.25;
    private File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    private String stardistOutput = "Label Image"; 
    public String stardistCellModel = "";
    
    public String[] cellsDetections = {"Stardist", "DOG"};
    public String cellsDetection = "";
    // dots threshold method
    private String thMet = "Moments";
    private String vesselThMet = "Moments";
    
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
       
     /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        if (chs == 3)
            vessel = true;
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal = new Calibration();  
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        return(cal);
    }
    
    public Calibration getCalib()
    {
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return cal;
    }
    
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
     /**
     * return objects population in an ClearBuffer image
     * @param imgCL
     * @return pop objects population
     */

    public Objects3DPopulation getPopFromClearBuffer(ClearCLBuffer imgBin, double min, double max) {
        ClearCLBuffer labelsCL = clij2.create(imgBin);
        clij2.connectedComponentsLabelingBox(imgBin, labelsCL);
        ClearCLBuffer labelsSizeFilter = clij2.create(imgBin);
        // filter size
        clij2.excludeLabelsOutsideSizeRange(labelsCL, labelsSizeFilter, min/(cal.pixelWidth*cal.pixelWidth*cal.pixelDepth),
                max/(cal.pixelWidth*cal.pixelWidth*cal.pixelDepth));
        clij2.release(labelsCL);
        ImagePlus img = clij2.pull(labelsSizeFilter);
        clij2.release(labelsSizeFilter);
        ImageHandler imh = ImageHandler.wrap(img);
        closeImages(img);
        Objects3DPopulation pop = new Objects3DPopulation(imh);
        pop.setCalibration(cal.pixelWidth, cal.pixelDepth, cal.getUnit());
        imh.closeImagePlus();
        return(pop);
    }  
    
    
   /**
     * compute local thickness
     * @param img
     * @return imgMap
    **/
    public ImageFloat localThickness3D (ImagePlus img, boolean inverse) {
        IJ.showStatus("Computing distance map...");
        img.setCalibration(cal);
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0,inverse, ThreadUtil.getNbCpus());
        return(edt);
    }
    
    
    /**
     * Find cell distance to nearest vessel with inverse distance map
     * @param cellPop
     * @param vesselPop
     * @return 
     */
    public ArrayList<Double> findCellVesselDist(Objects3DPopulation cellPop, ImagePlus imgVesselBin) {
        ArrayList<Double> cellDist = new ArrayList<>();
        ImageFloat imgMapInv = localThickness3D(imgVesselBin, true);
        int cellNb = cellPop.getNbObjects();
        for (Object3D cellObj : cellPop.getObjectsList()) {
            IJ.showStatus("Finding vessel distance for cell "+cellObj.getValue()+"/"+cellNb);
            Point3D pt = cellObj.getCenterAsPoint();
            double dist = imgMapInv.getPixel(pt);
            cellDist.add(dist);	
        }
        imgMapInv.closeImagePlus();
        return(cellDist);
    }
    
    
    
    /**
     * Tubeness
     * @param img
     * @return 
     */
    public ImagePlus tubeness(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img); 
        ClearCLBuffer imgCLDOG = DOG(imgCL, 5, 10);
        clij2.release(imgCL);        
        ImagePlus imgTube = clij2.pull(imgCLDOG); 
        clij2.release(imgCLDOG);
        return(imgTube);
    }
    
     /**
     * Find vessel diameter to nearest astrocyte
     * @param astroPop
     * @param imgVesselBin
     * @return 
     */
    public ArrayList<Double> findVesselDiameter(Objects3DPopulation astroPop, ImagePlus imgVesselBin) {
        ArrayList<Double> vesselDiam = new ArrayList<>();
        // get skeleton
        ImagePlus imgSkel = new Duplicator().run(imgVesselBin);
        IJ.setAutoThreshold(imgSkel, "Default dark");
        Prefs.blackBackground = false;
        IJ.run(imgSkel, "Convert to Mask", "method=Default background=Dark");
        IJ.run(imgSkel, "Skeletonize", "stack");
        ClearCLBuffer skelCL = clij2.push(imgSkel);
        Objects3DPopulation vesselSkelPop = getPopFromClearBuffer(skelCL, 0, Double.MAX_VALUE);
        clij2.release(skelCL);
        ImageFloat imgVesselmap = localThickness3D(imgVesselBin, false);
        int astroNb = astroPop.getNbObjects();
        for (Object3D astroObj : astroPop.getObjectsList()) {
            IJ.showStatus("Finding vessel diameter for cell "+astroObj.getValue()+"/"+astroNb);
            Object3D vesselObj = vesselSkelPop.closestBorder(astroObj);
            Voxel3D[] ptsBorder = astroObj.VoxelsBorderBorder(vesselObj);
            Point3D ptBorder = ptsBorder[1];
            double diam = imgVesselmap.getPixel(ptBorder);            
            vesselDiam.add(diam);
        }
        closeImages(imgSkel);
        return(vesselDiam);
    }
    
    public ImagePlus thresholdVessel(ImagePlus img, Roi roi) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer threshold = threshold(imgCL, vesselThMet); 
        clij2.release(imgCL);
        ImagePlus imgTh = clij2.pull(threshold);
        roi.setLocation(0, 0);
        IJ.setBackgroundColor(0, 0, 0);
        imgTh.setRoi(roi);
        IJ.run(imgTh, "Clear Outside", "stack");
        threshold = clij2.push(imgTh);
        ImagePlus imgFinal = clij2.pull(threshold);
        closeImages(imgTh);
        clij2.release(threshold);
        imgFinal.setCalibration(cal);
        return imgFinal;
    }
    
    public Objects3DPopulation findVessel(ImagePlus imgTh) {
        ClearCLBuffer threshold = clij2.push(imgTh);
        Objects3DPopulation vesselPop = getPopFromClearBuffer(threshold, 0, Double.MAX_VALUE);
        System.out.println(vesselPop.getNbObjects()+" vessels found");
        return(vesselPop);
    }
    
    
    /* Median filter 
     * Using CLIJ2
     * @param ClearCLBuffer
     * @param sizeXY
     * @param sizeZ
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer  imgCL, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        clij2.release(imgCL);
        return(imgCLMed);
    }
    
    /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double size1, double size2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        ClearCLBuffer imgCLDOGmed = median_filter(imgCLDOG, 2, 2);
        clij2.release(imgCLDOG);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOGmed, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        return(imgCLDOGmed);
    }
    
    /**
     * Fill hole
     * USING CLIJ2
     */
    private void fillHole(ClearCLBuffer imgCL) {
        long[] dims = clij2.getDimensions(imgCL);
        ClearCLBuffer slice = clij2.create(dims[0], dims[1]);
        ClearCLBuffer slice_filled = clij2.create(slice);
        for (int z = 0; z < dims[2]; z++) {
            clij2.copySlice(imgCL, slice, z);
            clij2.binaryFillHoles(slice, slice_filled);
            clij2.copySlice(slice_filled, imgCL, z);
        }
        clij2.release(slice);
        clij2.release(slice_filled);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     * @param fill 
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        return(imgCLBin);
    }
    
    
    /**
     * Find image type
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                    return fileExt;
                case "czi" :
                  return fileExt;
                case "lif"  :
                    return fileExt;
                case "isc2" :
                   return fileExt;
                default :
                   ext = fileExt;
                   break; 
            }
        }
        return(ext);
    }
    
     /*
    Find starDist models in Fiji models folder
    */
    public String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        Arrays.sort(models);
        return(models);
    } 
    
    /**
     * Write headers for results file
     * 
     * @param outDirResults
    */
    public void writeHeaders(String outDirResults) throws IOException {
        // global results
        FileWriter fileResults = new FileWriter(outDirResults +"Results.xls", false);
        outPutResults = new BufferedWriter(fileResults);
        outPutResults.write("Image name\tRoi name\tRoi volume\tNb Cell\tCell mean intensity\tCell sd intensity\t"
                + "Cell mean volume\t Cell sd volume\ttotal cell Volume\tCell average Distance To Closest Neighbor\tCell SD DistanceTo Closest Neighbor"
                + "\tCell Mean Of Average Distance of "+nbNei+" neighbors"+"\tCell SD Of Average Distance of "+nbNei+" neighbors"+"\tCell Mean Of Max Distance of "+nbNei+" neighbors"+
                "\tCell SD Of Max Distance of "+nbNei+" neighbors\tCell Area Curves\tCell Distribution Aleatoire Stat\tCell Mean distance to Vessel\tMean Vessel radius\n");
        outPutResults.flush();

        // distances results
        FileWriter fileDistances = new FileWriter(outDirResults +"Distances.xls", false);
        outPutDistances = new BufferedWriter(fileDistances);
        outPutDistances.write("Image name\tRoi name\tCell Volume\tCell Distance To Closest Neighbor\tCell Distance to Closest Vessel\tClosest Vessel radius\n");
        outPutDistances.flush();
    }
    
   /**
     * Dialog 
     * 
     * @param channels
     * @return 
     */
    public int[] dialog(String[] channels) {
        String[] models = findStardistModels();
        String[] chNames = new String[]{"Vessel1", "Vessel2", "Sox"};
        String[] methods = AutoThresholder.getMethods();
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName : chNames) {
            gd.addChoice(chName+" : ", channels, channels[index]);
            index++;
        }
        gd.addMessage("Cells parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min cell size (µm3) : ", minCell, 3);
        gd.addNumericField("Max cell size (µm3) : ", maxCell, 3);
        
        gd.addMessage("--- Stardist model ---", Font.getFont(Font.MONOSPACED), Color.blue);
        if (models.length > 0) {
            gd.addChoice("StarDist cell model :",models, models[0]);
        }
        else {
            gd.addMessage("No StarDist model found in Fiji !!", Font.getFont("Monospace"), Color.red);
            gd.addFileField("StarDist cell model :", stardistCellModel);
        }
        gd.addChoice("Cells detection method : ", cellsDetections, cellsDetection);
        gd.addMessage("--- DOG detection ---", Font.getFont(Font.MONOSPACED), Color.blue);
        gd.addChoice("Thresholding method :", methods, thMet);
        gd.addMessage("Difference of Gaussian (radius1 < radius2)", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("radius 1 (pixels) : ", sigma1, 1);
        gd.addNumericField("radius 2 (pixels) : ", sigma2, 1);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Calibration xy (µm)  :", cal.pixelWidth, 3);
        gd.addNumericField("Calibration z (µm)  :", cal.pixelDepth, 3);
        gd.addMessage("Spatial distribution", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Radius for neighboring analysis : ", radiusNei, 2);
        gd.addNumericField("Number of neighbors : ", nbNei, 2);
        gd.addCheckbox("Do comparaison to random distribution", false);
        gd.addCheckbox("Compute distance to vessel", false);
        gd.addChoice("Vessel thresholding method :", methods, vesselThMet);
        gd.showDialog();
        
        int[] chChoices = new int[chNames.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        minCell = gd.getNextNumber();
        maxCell = gd.getNextNumber();
        if (models.length > 0) {
            stardistCellModel = modelsPath+File.separator+gd.getNextChoice();
        }
        else {
            stardistCellModel = gd.getNextString();
        }
        cellsDetection = gd.getNextChoice();
        thMet = gd.getNextChoice();
        sigma1 = gd.getNextNumber();
        sigma2 = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = gd.getNextNumber();
        radiusNei = gd.getNextNumber();
        nbNei = (int)gd.getNextNumber();
        doF = gd.getNextBoolean();
        vessel = gd.getNextBoolean();
        vesselThMet = gd.getNextChoice();
        if (gd.wasCanceled())
                chChoices = null;
        return(chChoices);
    }
    
    
    /** 
     * Find cells with DOG method
     * @param img channel
     * @return cells population
     */
    public Objects3DPopulation findCellsDoG(ImagePlus img, Roi roi) {
        IJ.run(img, "Subtract Background...", "rolling=50 stack");
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = DOG(imgCL, sigma1, sigma2);
        clij2.release(imgCL);
        ClearCLBuffer imgCLBin = threshold(imgCLDOG, thMet); 
        clij2.release(imgCLDOG);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        imgBin.setRoi(roi);
        roi.setLocation(0, 0);
        IJ.setBackgroundColor(0, 0, 0);
        IJ.run(imgBin, "Clear Outside", "stack");
        Objects3DPopulation pmlPop = getPopFromClearBuffer(clij2.push(imgBin), minCell, maxCell);
        closeImages(imgBin);
        return(pmlPop);
    } 
    
      /** Look for all nuclei
    Do z slice by slice stardist 
    * return nuclei population
    */
   public Objects3DPopulation stardistNucleiPop(ImagePlus imgNuc, Roi roi) throws IOException{
       ImagePlus img =  new Duplicator().run(imgNuc);
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLM = median_filter(imgCL, 2, 2);
       clij2.release(imgCL);
       ImagePlus imgM = clij2.pull(imgCLM);
       clij2.release(imgCLM);
       closeImages(img);

       // Go StarDist
       File starDistModelFile = new File(stardistCellModel);
       StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
       star.loadInput(imgM);
       star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThreshNuc, stardistOverlayThreshNuc, stardistOutput);
       star.run();
       closeImages(imgM);
       // label in 3D
       ImagePlus nuclei = new Duplicator().run(star.getLabelImagePlus());
       nuclei.setRoi(roi);
       roi.setLocation(0, 0);
       IJ.setBackgroundColor(0, 0, 0);
       IJ.run(nuclei, "Clear Outside", "stack");
       imgCL = clij2.push(nuclei);
       Objects3DPopulation nPop = getPopFromClearBuffer(imgCL, minCell, maxCell); 
       clij2.release(imgCL);
       closeImages(nuclei);
       return(nPop);
   }
    
    
    /**
     * Create cells population image
     * @param cellPop
     * @param img
     * @param pathName

     */
    public void saveCellsImage(Objects3DPopulation cellPop, Objects3DPopulation vesselPop, ImagePlus img, ArrayList<Double> dists, String pathName) {
        ImageHandler imhCell = ImageHandler.wrap(img).createSameDimensions();
        if (vessel) {
            double maxDist = Collections.max(dists);
            for (int i = 0; i < cellPop.getNbObjects(); i++) {
                Object3D obj = cellPop.getObject(i);
                double dist = dists.get(i);
                obj.draw(imhCell, (float)dist);
            }
            vesselPop.draw(imhCell, (int)maxDist+10);
        }
        else
            cellPop.draw(imhCell);
        
        // save image for objects population
        ImagePlus imgCells = imhCell.getImagePlus();
        imgCells.setCalibration(cal);
        if (vessel) {
            IJ.run(imgCells, "Fire", "");
            IJ.run(imgCells,"Calibrate...","function=None unit="+cal.getUnit());
            IJ.run(imgCells, "Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=2 font=12 zoom=1 overlay show");
        }

        FileSaver ImgObjectsFile = new FileSaver(imgCells);
        ImgObjectsFile.saveAsTiff(pathName); 
        closeImages(imgCells);
        imhCell.closeImagePlus();
    }
    
    /** \brief bounding box of the population*/
     protected Object3D maskBounding(Objects3DPopulation pop, ImagePlus imgCells) {
         // change to convex hull ?
        ImageHandler imh = ImageHandler.wrap(imgCells).createSameDimensions();
        pop.draw(imh, 255);
        Object3D objMask = createObject3DVoxels(imh.getImagePlus(), 255);
        imh.closeImagePlus();
        int[] bbMask = objMask.getBoundingBox();
        ImagePlus imgMask = new Duplicator().run(imgCells);
        IJ.run(imgMask,"8-bit","");
        Roi roi = new Roi(bbMask[0], bbMask[2], bbMask[1]-bbMask[0], bbMask[3] - bbMask[2]);
        imgMask.setRoi(roi);
        IJ.setForegroundColor(255,255,255);
        IJ.run(imgMask, "Fill", "stack");
        IJ.setBackgroundColor(0, 0, 0);
        IJ.run(imgMask, "Clear Outside", "stack");
        imgMask.deleteRoi();
        objMask = createObject3DVoxels(imgMask, 255);
        return(objMask);
    }
     
    
    /**
    * For compute G function, with parallel version
    * 
    * 
     * @param pop
     * @param mask
     * @return G SDI
    **/ 
    public double[] processGParallel (Objects3DPopulation pop, ImagePlus imgCells, String outDirResults, String imgName, String roiName, Roi roi) {
        
        // change to convex hull ?
        Object3D objMask = maskBounding(pop, imgCells);
        String outname = outDirResults + imgName + "_Gplot_" + roiName + ".tif";
        double minDist = Math.pow(3.0/4.0/Math.PI*minCell, 1.0/3.0) * 2; // calculate minimum cell radius -> *2 min distance
        return processG(pop, objMask, 50, minDist, outname, roiName);
  
    }
    
    private double[] processG(Objects3DPopulation pop, Object3D mask, final int numRandomSamples, final double distHardCore, String filename, String roiName) {
        //final Calibration calibration = Object3D_IJUtils.getCalibration(mask);
        final double sxy = mask.getResXY();
        final double sz = mask.getResZ();
        final String unit = mask.getUnits();
        //final Calibration calibration = mask.getCalibration();
        final int nbSpots = pop.getNbObjects();

        // observed G
        ArrayUtil observedDistancesG;
        ArrayUtil observedCDG;
        observedDistancesG = pop.distancesAllClosestCenter();
        observedDistancesG.sort();
        observedCDG = CDFTools.cdf(observedDistancesG);

        // Average G
        ArrayUtil xEvalG;
        final ArrayUtil[] sampleDistancesG;
        ArrayUtil averageCDG;

        xEvalG = new ArrayUtil(numRandomSamples * nbSpots);
        sampleDistancesG = new ArrayUtil[numRandomSamples];

        // PARALLEL
        final Object3D mask2 = mask;
        final AtomicInteger ai = new AtomicInteger(0);
        final int nCpu = ThreadUtil.getNbCpus();
        Thread[] threads = ThreadUtil.createThreadArray(nCpu);
        final int dec = (int) Math.ceil((double) numRandomSamples / (double) nCpu);
        for (int iThread = 0; iThread < threads.length; iThread++) {
            threads[iThread] = new Thread() {
                @Override
                public void run() {
                    ArrayUtil distances2;
                    //image.setShowStatus(show);
                    for (int k = ai.getAndIncrement(); k < nCpu; k = ai.getAndIncrement()) {
                        for (int i = dec * k; ((i < (dec * (k + 1))) && (i < numRandomSamples)); i++) {
                            Objects3DPopulation popRandom = new Objects3DPopulation();
                            //popRandom.setCalibration(calibration);
                            popRandom.setCalibration(sxy, sz, unit);
                            popRandom.setMask(mask2);
                            popRandom.createRandomPopulation(nbSpots, distHardCore);
                            distances2 = popRandom.distancesAllClosestCenter();
                            distances2.sort();
                            sampleDistancesG[i] = distances2;
                        }
                    }
                }
            };
        }
        ThreadUtil.startAndJoin(threads);
        for (int i = 0; i < numRandomSamples; i++) {
            xEvalG.insertValues(i * nbSpots, sampleDistancesG[i]);
        }
        xEvalG.sort();
        averageCDG = CDFTools.cdfAverage(sampleDistancesG, xEvalG);

        // Envelope G
        ai.set(0);
        for (int iThread = 0; iThread < threads.length; iThread++) {
            threads[iThread] = new Thread() {
                @Override
                public void run() {
                    ArrayUtil distances2;
                    //image.setShowStatus(show);
                    for (int k = ai.getAndIncrement(); k < nCpu; k = ai.getAndIncrement()) {
                        for (int i = dec * k; ((i < (dec * (k + 1))) && (i < numRandomSamples)); i++) {
                            
                            Objects3DPopulation popRandom = new Objects3DPopulation();
                            //popRandom.setCalibration(calibration);
                            popRandom.setCalibration(sxy, sz, unit);
                            popRandom.setMask(mask2);
                            popRandom.createRandomPopulation(nbSpots, distHardCore);
                            distances2 = popRandom.distancesAllClosestCenter();
                            distances2.sort();
                            sampleDistancesG[i] = distances2;
                        }
                    }
                }
            };
        }
        ThreadUtil.startAndJoin(threads);

        double sdi_G = CDFTools.SDI(observedDistancesG, sampleDistancesG, averageCDG, xEvalG);
        double area = integrate(observedDistancesG,observedCDG)-integrate(xEvalG,averageCDG)+(xEvalG.getMaximum() - observedDistancesG.getMaximum());

        // plot
        Plot plotG = null;
        plotG = createPlot(xEvalG, sampleDistancesG, observedDistancesG, observedCDG, averageCDG, "G");
        plotG.draw();
        
        plotG.addLabel(0.1, 0.1, "sdi = " + String.format("%.2f",sdi_G));
        plotG.addLabel(0.1, 0.15, "Area diff. = " + String.format("%.3f",area));
        ImagePlus imgPlot = plotG.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(filename);
        closeImages(imgPlot); 
        //showPlotData(roiName, filename,observedDistancesG,observedCDG,xEvalG,averageCDG, plotResults);
        double[] res = {sdi_G, area}; 
        return res;
    }

    public void showPlotData(String roiName,String imgName,ArrayUtil xObs, ArrayUtil yObs, ArrayUtil xRand, ArrayUtil yRand,BufferedWriter plotResults) 
            throws IOException{
        ArrayUtil a = (xObs.size()>xRand.size()) ? xObs : xRand;
        for (int i=0; i<a.size(); i++){          
            plotResults.write(imgName+"\t"+roiName+"\t"+xObs.getValue(i)+"\t"+yObs.getValue(i)+"\t"+xRand.getValue(i)+"\t"+yRand.getValue(i)+"\n");
            plotResults.flush();
        }
    }

     private Plot createPlot(ArrayUtil xEvals, ArrayUtil[] sampleDistances, ArrayUtil observedDistances, ArrayUtil observedCD, ArrayUtil averageCD, String function) {
        Color ColorAVG = Color.red;
        Color ColorENV = Color.green;
        Color ColorOBS = Color.blue;
        int nbBins = 1000;
        double env = 0.25;
        double plotMaxX = observedDistances.getMaximum();
        double plotMaxY = observedCD.getMaximum();

        // low env
        double max = xEvals.getMaximum();
        ArrayUtil xEval0 = new ArrayUtil(nbBins);
        for (int i = 0; i < nbBins; i++) {
            xEval0.addValue(i, ((double) i) * max / ((double) nbBins));
        }
        // get the values
        ArrayUtil samplesPc5 = CDFTools.cdfPercentage(sampleDistances, xEval0, env / 2.0);
        ArrayUtil samplesPc95 = CDFTools.cdfPercentage(sampleDistances, xEval0, 1.0 - env / 2.0);
        // get the limits
        if (xEval0.getMaximum() > plotMaxX) {
            plotMaxX = xEval0.getMaximum();
        }
        if (samplesPc5.getMaximum() > plotMaxY) {
            plotMaxY = samplesPc5.getMaximum();
        }
        if (samplesPc95.getMaximum() > plotMaxY) {
            plotMaxY = samplesPc95.getMaximum();
        }
        if (xEvals.getMaximum() > plotMaxX) {
            plotMaxX = xEvals.getMaximum();
        }
        if (averageCD.getMaximum() > plotMaxY) {
            plotMaxY = averageCD.getMaximum();
        }
        if (observedCD.getMaximum() > plotMaxY) {
            plotMaxY = observedCD.getMaximum();
        }
        if (observedDistances.getMaximum() > plotMaxX) {
            plotMaxX = observedDistances.getMaximum();
        }
        // create the plot
        Plot plot = new Plot(function + "-function", "distance", "cumulated frequency");
        plot.setLimits(0, plotMaxX, 0, plotMaxY);

        // envelope  for e.g 10 % at 5 and 95 %
        plot.setColor(ColorENV);
        plot.addPoints(xEval0.getArray(), samplesPc5.getArray(), Plot.LINE);

        // envelope  for e.g 10 % at 5 and 95 %
        plot.setColor(ColorENV);
        plot.addPoints(xEval0.getArray(), samplesPc95.getArray(), Plot.LINE);

        // average
        plot.setColor(ColorAVG);
        plot.addPoints(xEvals.getArray(), averageCD.getArray(), Plot.LINE);

        // observed
        plot.setColor(ColorOBS);
        plot.addPoints(observedDistances.getArray(), observedCD.getArray(), Plot.LINE);

        return plot;
    }

    
    public void getDistNeighbors(Object3D obj, Objects3DPopulation pop, DescriptiveStatistics cellNbNeighborsDistMean, 
           DescriptiveStatistics cellNbNeighborsDistMax) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        double[] dist= pop.kClosestDistancesSquared(obj.getCenterX(), obj.getCenterY(), obj.getCenterZ(), nbNei);
        for (double d : dist) {
           stats.addValue(Math.sqrt(d)); 
        }
        cellNbNeighborsDistMax.addValue(stats.getMax());
        cellNbNeighborsDistMean.addValue(stats.getMean());
    }
    
    /**
     * 
     * @param x
     * @param y
     * @return 
     */
    private double integrate (ArrayUtil x, ArrayUtil y){
        double sum = 0.0;
        for(int i=0; i<x.size()-1; i++){
            sum += (x.getValue(i+1)-x.getValue(i))*y.getValue(i);
            sum += (x.getValue(i+1)-x.getValue(i))*(y.getValue(i+1)-y.getValue(i))/2.0;
        }
        return sum;
    }
    
    /**
    * Compute global cells parameters
    * @param cellPop cell population
    * @param imgCell read cell intensity
     * @param roiName
     * @param roi
    * @param imgName image file
    * @param outDirResults results file
     * @throws java.io.IOException
    **/
    public void computeNucParameters(Objects3DPopulation cellPop, Objects3DPopulation vesselPop, ArrayList<Double> dist, ArrayList<Double> diam, ImagePlus imgCell, String roiName, Roi roi,
        String imgName, String outDirResults) throws IOException {
        
        DescriptiveStatistics cellIntensity = new DescriptiveStatistics();
        DescriptiveStatistics cellVolume = new DescriptiveStatistics();
        DescriptiveStatistics cellNbNeighborsDistMean = new DescriptiveStatistics();
        DescriptiveStatistics cellNbNeighborsDistMax = new DescriptiveStatistics();
        double cellVolumeSum = 0.0;
        // do individual stats
        cellPop.createKDTreeCenters();
        ArrayUtil alldistances = cellPop.distancesAllClosestCenter();
        for (int i = 0; i < cellPop.getNbObjects(); i++) {
            IJ.showStatus("Computing cell "+i+" parameters ....");
            IJ.showProgress(i, cellPop.getNbObjects());
            Object3D cellObj = cellPop.getObject(i);
            cellIntensity.addValue(cellObj.getIntegratedDensity(ImageHandler.wrap(imgCell)));
            cellVolume.addValue(cellObj.getVolumeUnit());
            cellVolumeSum += cellObj.getVolumeUnit();
            double vesselDist = Double.NaN;
             double vesselDiam = Double.NaN;
            if (vessel) {
                Object3D vesselObj = vesselPop.closestBorder(cellObj);
                vesselDist = dist.get(i);
                vesselDiam = diam.get(i);
                outPutDistances.write(imgName+"\t"+roiName+"\t"+cellObj.getVolumeUnit()+"\t"+alldistances.getValue(i)+"\t"+vesselDist+"\t"+vesselDiam+"\n");
                outPutDistances.flush();
            }
            else {
                outPutDistances.write(imgName+"\t"+roiName+"\t"+cellObj.getVolumeUnit()+"\t"+alldistances.getValue(i)+"\t"+"-\n");
                outPutDistances.flush();
            }
            getDistNeighbors(cellObj, cellPop, cellNbNeighborsDistMean, cellNbNeighborsDistMax);
        }
        
        double sdiF = Double.NaN;
        double areaCurves = Double.NaN;
        if (doF) {
            IJ.showStatus("Computing spatial distribution G Function ...");
            double[] temp = processGParallel(cellPop, imgCell, outDirResults, imgName, roiName, roi);
            sdiF = temp[0];
            areaCurves = temp[1];
        }
        double vesselDistMean = 0;
        double vesselDiamMean = 0;
        if (vessel) {
            vesselDistMean = dist.stream().mapToDouble(val -> val).average().orElse(0.0);
            vesselDiamMean = diam.stream().mapToDouble(val -> val).average().orElse(0.0);
        }
        double minDistCenterMean = alldistances.getMean(); 
        double minDistCenterSD = alldistances.getStdDev();
        // compute statistics
        double cellIntMean = cellIntensity.getMean();
        double cellIntSD = cellIntensity.getStandardDeviation();
        double cellVolumeMean = cellVolumeSum/cellPop.getNbObjects();
        double cellVolumeSD = cellVolume.getStandardDeviation();
        double roiVol = imgCell.getWidth()*imgCell.getHeight()*imgCell.getNSlices()*imgCell.getCalibration().pixelDepth;
        double neiMeanMean = cellNbNeighborsDistMean.getMean();
        double neiMeanSD = cellNbNeighborsDistMean.getStandardDeviation();
        double neiMaxMean = cellNbNeighborsDistMax.getMean();
        double neiMaxSD = cellNbNeighborsDistMax.getStandardDeviation();
       
        outPutResults.write(imgName+"\t"+roiName+"\t"+roiVol+"\t"+cellPop.getNbObjects()+"\t"+cellIntMean+"\t"+cellIntSD+"\t"+cellVolumeMean+"\t"+
                cellVolumeSD+"\t"+cellVolumeSum+"\t"+minDistCenterMean+"\t"+minDistCenterSD+"\t"+neiMeanMean+"\t"+neiMeanSD+"\t"+neiMaxMean+"\t"+neiMaxSD
                +"\t"+areaCurves+"\t"+sdiF+"\t"+vesselDistMean+"\t"+vesselDiamMean+"\n");
        
        outPutResults.flush();
    }
    
    
    public void closeResults() throws IOException {
       outPutResults.close();
       outPutDistances.close();         
       //outPutPlot.close();
    }
    
    public int getPyramidalFactor(ImageProcessorReader reader) {
        reader.setSeries(0);
        int sizeXseries0 = reader.getSizeX();
        reader.setSeries(reader.getSeriesCount()-1);
        int sizeXseriesN = reader.getSizeX();
        return (sizeXseries0 / sizeXseriesN); 
    }
    
    public Roi scaleRoi(Roi roi, int scale) {
        Roi scaledRoi = new RoiScaler().scale(roi, scale, scale, false);
        Rectangle rect = roi.getBounds();
        scaledRoi.setLocation(rect.x*scale, rect.y*scale);
        return scaledRoi;
    }
}
