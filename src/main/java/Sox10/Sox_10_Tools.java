package Sox10;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
import mcib3d.spatial.analysis.SpatialStatistics;
import mcib3d.spatial.descriptors.F_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
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
    public double minCell = 0.02;
    // max volume in µm^3 for cells
    public double maxCell = 60;
    // Dog parameters
    private double radius = 2;
    private Calibration cal = new Calibration(); 
    

    // dots threshold method
    private String thMet = "Otsu";
    
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
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        return(imgCLDOG);
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
        gd.addMessage("Cells parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min cell size (µm3) : ", minCell, 3);
        gd.addNumericField("Max cell size (µm3) : ", maxCell, 3);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Calibration xy (µm)  :", cal.pixelWidth, 3);
        if ( cal.pixelDepth == 1) cal.pixelDepth = 0.5;
        gd.addNumericField("Calibration z (µm)  :", cal.pixelDepth, 3);
        gd.showDialog();
        int[] chChoices = new int[channels.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        minCell = gd.getNextNumber();
        maxCell = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = gd.getNextNumber();
        if (gd.wasCanceled())
                chChoices = null;
        return(chChoices);
    }
    

    
    /** 
     * Find cells with DOG method
     * @param img channel
     * @return cells population
     */
    public Objects3DPopulation findCellsDoG(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        double sig1 = radius/(3*img.getCalibration().pixelWidth);
        double sig2 = radius/img.getCalibration().pixelWidth;
        ClearCLBuffer imgCLDOG = DOG(imgCL, sig1, sig2);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(threshold(imgCLDOG, thMet, false));
        clij2.release(imgCLDOG);
        imgBin.setCalibration(img.getCalibration());
        //imgBin.show();
        //new WaitForUserDialog("test").show();
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
        
    /**
    * For compute F function
     * @param pop
     * @param mask
     * @return F SDI
    **/ 
    public double processF (Objects3DPopulation pop, ImagePlus imgCells, String outDirResults, String imgName, String roiName) {
        Object3D mask = createObject3DVoxels(imgCells, 0);
        
        // define spatial descriptor, model
        SpatialDescriptor spatialDesc = new F_Function(pop.getNbObjects(), mask);
        pop.setMask(mask);
        SpatialModel spatialModel = new SpatialRandomHardCore(pop.getNbObjects(), 0.8, mask);
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, 50, pop);
        spatialStatistics.setEnvelope(0.25);
        spatialStatistics.setVerbose(false);
        Plot fPlot = spatialStatistics.getPlot();
        fPlot.draw();
        fPlot.addLabel(0.1, 0.1, "p = " + String.valueOf(spatialStatistics.getSdi()));
        ImagePlus imgPlot = fPlot.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(outDirResults + imgName + "_Fplot_" + roiName + ".tif");
        closeImages(imgPlot);
        System.out.println("Sdi = " + spatialStatistics.getSdi());
        return(spatialStatistics.getSdi());
    }
    
    
    
    /**
    * Compute global cells parameters
    * @param cellPop cell population
    * @param imgCell read cell intensity
    * @param imgName image file
    * @param outDirResults results file
     * @param results buffer
    **/
    public void computeNucParameters(Objects3DPopulation cellPop, ImagePlus imgCell, String roiName,
            String imgName, String outDirResults, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing cells parameters ....");
        
        DescriptiveStatistics cellIntensity = new DescriptiveStatistics();
        DescriptiveStatistics cellVolume = new DescriptiveStatistics();
        double minDistCenterMean = Double.NaN;
        double minDistCenterSD = Double.NaN;
        double sdiF = Double.NaN;
        for (int i = 0; i < cellPop.getNbObjects(); i++) {
            Object3D cellObj = cellPop.getObject(i);
            cellIntensity.addValue(cellObj.getIntegratedDensity(ImageHandler.wrap(imgCell)));
            cellVolume.addValue(cellObj.getVolumeUnit());
        }
        sdiF = processF(cellPop, imgCell, outDirResults, imgName, roiName);

        minDistCenterMean = cellPop.distancesAllClosestCenter().getMean(); 
        minDistCenterSD = cellPop.distancesAllClosestCenter().getStdDev();
        // compute statistics
        double cellIntMean = cellIntensity.getMean();
        double cellIntSD = cellIntensity.getStandardDeviation();
        double cellVolumeMean = cellVolume.getMean();
        double cellVolumeSD = cellVolume.getStandardDeviation();
        double cellVolumeSum = cellVolume.getSum();
        double roiVol = imgCell.getWidth()*imgCell.getHeight()*imgCell.getNSlices()*imgCell.getCalibration().pixelDepth;  
        results.write(imgName+"\t"+roiName+"\t"+roiVol+"\t"+cellPop.getNbObjects()+"\t"+cellIntMean+"\t"+cellIntSD+"\t"+cellVolumeMean+
                cellVolumeSD+"\t"+cellVolumeSum+"\t"+sdiF+"\n");
        
        results.flush();
    }
}