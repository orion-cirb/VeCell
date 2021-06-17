package Sox10;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
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
import mcib3d.geom.Object3DVoxels;
import static mcib3d.geom.Object3D_IJUtils.createObject3DVoxels;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.spatial.analysis.SpatialStatistics;
import mcib3d.spatial.descriptors.F_Function;
import mcib3d.spatial.descriptors.G_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;
import mcib3d.utils.ThreadUtil;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

        
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
    public double minCell = 200;
    // max volume in µm^3 for cells
    public double maxCell = 3000;
    // sigmas for DoG
    private double sigma1 = 2;
    private double sigma2 = 4;
    // Dog parameters
    private Calibration cal = new Calibration(); 
    
    private boolean doF = false;
    private double radiusNei = 50; //neighboring radius
    

    // dots threshold method
    private String thMet = "Moments";
    
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
        try {
            loader.loadClass("uk.ac.sussex.gdsc.utils.DifferenceOfGaussians_PlugIn");
        } catch (ClassNotFoundException e) {
            IJ.log("GDSC Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("TurboReg_");
        } catch (ClassNotFoundException e) {
            IJ.log("TurboReg not installed, please install from http://bigwww.epfl.ch/thevenaz/turboreg/");
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
    public static String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
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
                    if (meta.getChannelID(0, n) == null)
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
                    channels[0] = Integer.toString(n);
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
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
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
    
    public Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
   
    
    
    /* Median filter 
     * Using CLIJ2
     * @param ClearCLBuffer
     * @param sizeXY
     * @param sizeZ
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer  imgCL, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.mean3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
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
        long[] dims = clij2.getSize(imgCL);
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
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed, boolean fill) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        if (fill)
            fillHole(imgCLBin);
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
    
    
    /**
     * Dialog 
     * 
     * @param channels
     * @return 
     */
    public int[] dialog(String[] channels) {
        String[] chNames = {"Nucleus", "TF"};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName : channels) {
            gd.addChoice(chNames[index]+" : ", channels, channels[index]);
            index++;
        }
        String[] methods = AutoThresholder.getMethods();
        gd.addMessage("Cells parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min cell size (µm3) : ", minCell, 3);
        gd.addNumericField("Max cell size (µm3) : ", maxCell, 3);
        gd.addChoice("Thresholding method :", methods, thMet);
        gd.addMessage("Difference of Gaussian", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("radius 1 (pixels) : ", sigma1, 1);
        gd.addNumericField("radius 2 (pixels) : ", sigma2, 1);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Calibration xy (µm)  :", cal.pixelWidth, 3);
        if ( cal.pixelDepth == 1) cal.pixelDepth = 0.5;
        gd.addNumericField("Calibration z (µm)  :", cal.pixelDepth, 3);
        gd.addMessage("Spatial distribution", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Radius for neighboring analysis", radiusNei, 2);
        gd.addCheckbox("Do comparaison to random distribution", false);
        gd.showDialog();
        int[] chChoices = new int[channels.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        minCell = gd.getNextNumber();
        maxCell = gd.getNextNumber();
        thMet = gd.getNextChoice();
        sigma1 = gd.getNextNumber();
        sigma2 = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = gd.getNextNumber();
        radiusNei = gd.getNextNumber();
        doF = gd.getNextBoolean();
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
        ClearCLBuffer imgCLBin = threshold(imgCLDOG, thMet, false); 
        clij2.release(imgCLDOG);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        imgBin.setRoi(roi);
        roi.setLocation(0, 0);
        IJ.setBackgroundColor(0, 0, 0);
        IJ.run(imgBin, "Clear Outside", "stack");
        imgBin.setCalibration(img.getCalibration());
        Objects3DPopulation pmlPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minCell, maxCell, true));
        closeImages(imgBin);
        return(pmlPop);
    } 
    

    
    /**
     * Create cells population image
     * @param cellPop
     * @param img
     * @param pathName

     */
    public void saveCellsImage(Objects3DPopulation cellPop, ImagePlus img, String pathName) {
        ImageHandler imh = ImageHandler.wrap(img);
        cellPop.draw(imh, 255);

        // Save diffus
        FileSaver imgCells = new FileSaver(imh.getImagePlus());
        imgCells.saveAsTiff(pathName);
        imh.closeImagePlus();
    }
    
    /** \brief bounding box of the population*/
     protected void maskBounding(Objects3DPopulation pop, ImagePlus empty, Roi roi) {
         int[] res = {Integer.MAX_VALUE, 0, Integer.MAX_VALUE, 0, Integer.MAX_VALUE, 0};
        
         // get min, max z
        for (Object3D obj : pop.getObjectsList()) {
            int[] bb = obj.getBoundingBox();
            if (res[4] > bb[4]) res[4]= bb[4];
            if (res[5] < bb[5]) res[5]= bb[5];
        }
        
        for (int z=res[4]; z<=res[5]; z++)
        {
            empty.setSlice(z);
            ImageProcessor proc = empty.getProcessor();
            empty.setRoi(roi);
            IJ.setForegroundColor(255,255,255);
            IJ.run(empty, "Fill", "slice");
       }
       IJ.run(empty, "Select None", "");
    }
    
    /**
    * For compute G function, with parallel version
    * 
    * 
     * @param pop
     * @param mask
     * @return G SDI
    **/ 
    public double processGParallel (Objects3DPopulation pop, ImagePlus imgCells, String outDirResults, String imgName, String roiName, Roi roi) {
        
        // change to convex hull ?
        ImagePlus imgMask = new Duplicator().run(imgCells);
        IJ.run(imgMask, "Select All", "");
        IJ.run(imgMask, "Clear", "stack");
        IJ.run(imgMask, "Select None", "");
        maskBounding(pop, imgMask, roi);
        //imgMask.show();
        //new WaitForUserDialog("test").show();
        Object3D mask = createObject3DVoxels(imgMask, 255);
        closeImages(imgMask);
        String outname = outDirResults + imgName + "_Gplot_" + roiName + ".tif";
        double minDist = Math.pow(3.0/4.0/Math.PI*minCell, 1.0/3.0) * 2; // calculate minimum cell radius -> *2 min distance
        return processG(pop, mask, 50, minDist, outname);
  
    }
    
    private double processG(Objects3DPopulation pop, Object3D mask, final int numRandomSamples, final double distHardCore, String filename) {
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
        IJ.log("SDI G=" + sdi_G);
        // plot
        Plot plotG = null;
        plotG = createPlot(xEvalG, sampleDistancesG, observedDistancesG, observedCDG, averageCDG, "G");
        plotG.draw();
        
        plotG.addLabel(0.1, 0.1, "sdi = " + String.valueOf(sdi_G));
        ImagePlus imgPlot = plotG.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(filename);
        closeImages(imgPlot);

        return sdi_G;
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

    
    /**
     * Return dilated object restriced to image borders
     * @param img
     * @param obj
     * @return
     */
    private double getVolumeObjInside(ImagePlus img, Object3D obj) {
        double x = obj.getCenterX();
        double y = obj.getCenterY();
        double z = obj.getCenterZ();
        Object3DVoxels sphere = new Object3DVoxels();
        sphere.createSphereUnit((float)x, (float)y, (float)z, (float)radiusNei);
        // check if object go outside image
        if (sphere.getXmin() < 0 || sphere.getXmax() > img.getWidth() || sphere.getYmin() < 0 || sphere.getYmax() > img.getHeight()
                || sphere.getZmin() < 0 || sphere.getZmax() > img.getNSlices()) {
            Object3DVoxels voxObj = new Object3DVoxels(sphere.listVoxels(ImageHandler.wrap(img)));
            return(voxObj.getVolumeUnit());
        }
        else
            return(sphere.getVolumeUnit());
    }
    
    public double getNbNeighbors(Object3D obj, Objects3DPopulation pop, ImagePlus img) {
        double vol = getVolumeObjInside(img, obj);
        ArrayList<Object3D> objs = pop.getObjectsWithinDistanceCenter(obj, radiusNei);
        int nobj = objs.size();
        return nobj/vol;
    }
    
    /**
    * Compute global cells parameters
    * @param cellPop cell population
    * @param imgCell read cell intensity
    * @param imgName image file
    * @param outDirResults results file
     * @param results buffer
    **/
    public void computeNucParameters(Objects3DPopulation cellPop, ImagePlus imgCell, String roiName, Roi roi,
          String imgName, String outDirResults, BufferedWriter results, BufferedWriter distances) throws IOException {
        
        DescriptiveStatistics cellIntensity = new DescriptiveStatistics();
        DescriptiveStatistics cellVolume = new DescriptiveStatistics();
        DescriptiveStatistics cellNbNeighbors = new DescriptiveStatistics();
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
            cellNbNeighbors.addValue( getNbNeighbors(cellObj, cellPop, imgCell) );
            distances.write(imgName+"\t"+roiName+"\t"+cellObj.getVolumeUnit()+"\t"+alldistances.getValue(i)+"\n");
        }
        distances.flush();
        double sdiF = Double.NaN;
        if (doF) {
            IJ.showStatus("Computing spatial distribution G Function ...");
            sdiF = processGParallel(cellPop, imgCell, outDirResults, imgName, roiName, roi);
        }
        
        
        double minDistCenterMean = alldistances.getMean(); 
        double minDistCenterSD = alldistances.getStdDev();
        // compute statistics
        double cellIntMean = cellIntensity.getMean();
        double cellIntSD = cellIntensity.getStandardDeviation();
        double cellVolumeMean = cellVolumeSum/cellPop.getNbObjects();
        double cellVolumeSD = cellVolume.getStandardDeviation();
        double roiVol = imgCell.getWidth()*imgCell.getHeight()*imgCell.getNSlices()*imgCell.getCalibration().pixelDepth;
        double neiMean = cellNbNeighbors.getMean();
        double neiSD = cellNbNeighbors.getStandardDeviation();
        
        results.write(imgName+"\t"+roiName+"\t"+roiVol+"\t"+cellPop.getNbObjects()+"\t"+cellIntMean+"\t"+cellIntSD+"\t"+cellVolumeMean+"\t"+
                cellVolumeSD+"\t"+cellVolumeSum+"\t"+minDistCenterMean+"\t"+minDistCenterSD+"\t"+neiMean+"\t"+neiSD+"\t"+sdiF+"\n");
        
        results.flush();
    }
}