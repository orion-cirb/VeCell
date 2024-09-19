package VeCell_Tools;

/* Copyright 2006, 2007 Mark Longair */

/*
This file is part of the ImageJ plugin "Quantile Based Normalization".

The ImageJ plugin "Quantile Based Normalization" is free software;
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option)
any later version.

The ImageJ plugin "Quantile Based Normalization" is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

In addition, as a special exception, the copyright holders give
you permission to combine this program with free software programs or
libraries that are released under the Apache Public License.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import ij.process.ShortProcessor;
import java.awt.image.ColorModel;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import org.apache.commons.io.FilenameUtils;
import util.BatchOpener;

import vib.app.FileGroup;

public class QuantileBasedNormalization {
// implements PlugIn, ActionListener, ItemListener

    /* The idea of this normalization is to rank all of the values
       in an several images, divide each of lists of ranked values
       up into a number of quantiles and replace each value in
       each quantile with the mean of that rank across all the
       images.  If you replace with the rank instead of the mean
       then you get histogram equalization of all the images. */

    /* The only subtlety with the former method is that the mean
       is not going to be integral, so we randomly replace the
       values in a quatile with a mix of the two bytes around the
       mean in such a proportion that the mean will be a close to
       correct as we can get it. */


    class Replacements {

            HashMap<Integer, Long> replacements;
            long totalReplacements;
            int minReplacement = Integer.MAX_VALUE;
            int maxReplacement = Integer.MIN_VALUE;
            Random rng;
            int quantile;

            public Replacements(int possibleValues) {
                    replacements = new HashMap<Integer, Long>();
                    rng=new Random();
            }

            public void addSomeReplacements( long howManyToReplace, int replacement ) {
                    if( replacement < minReplacement )
                            minReplacement = replacement;
                    if( replacement > maxReplacement )
                            maxReplacement = replacement;
                    if( ! replacements.containsKey(replacement) ) {
                            replacements.put(replacement, 0L);
                    }
                    long previousValue = replacements.get(replacement);
                    replacements.put(replacement, previousValue + howManyToReplace);
                    totalReplacements += howManyToReplace;
            }

            public int getRandomReplacement() {
                    if( totalReplacements == 0 ) {
                            return -1;
                    }

                    long index=Math.abs(rng.nextLong()) % totalReplacements;

                    long replacementsSkipped = 0;

                    for( int r = minReplacement; r <= maxReplacement; ++r ) {

                            if ( ! replacements.containsKey(r) )
                                    continue;

                            long indexInThisSlot = index - replacementsSkipped;

                            long numberOfReplacements = replacements.get(r);
                            if( indexInThisSlot < numberOfReplacements ) {
                                    // Then we remove one of these and return
                                    // the value of r.
                                    replacements.put(r, numberOfReplacements - 1);
                                    -- totalReplacements;
                                    return r;
                            } else {
                                    replacementsSkipped += numberOfReplacements;
                            }
                    }
                    return -1;
            }

            @Override
            public String toString() {
                    if( totalReplacements == 0 )
                            return "No replacements left.";
                    else {
                            String result = "" + totalReplacements + " replacements left (in";
                            for( int i = minReplacement; i <= maxReplacement; ++i ) {
                                    long numberOfReplacements = replacements.get(i);
                                    if( numberOfReplacements > 0 )
                                            result += " " + i + " (" + numberOfReplacements + ")";
                            }
                            return result;
                    }

            }

    }

    public void divideIntoQuantiles(int numberOfQuantiles,
                                    long frequencies[],
                                    long pointsInImage,
                                    long [] resultSumValuesInQuantile,
                                    long [] resultNumberOfValuesInQuantile) {

            if (numberOfQuantiles != resultNumberOfValuesInQuantile.length)
                    throw new RuntimeException("BUG: numberOfQuantiles didn't match resultNumberOfValuesInQuantile.length");
            if (numberOfQuantiles != resultSumValuesInQuantile.length)
                    throw new RuntimeException("BUG: numberOfQuantiles didn't match resultSumValuesInQuantile.length");

            for (int q = 0; q < numberOfQuantiles; ++q) {

                    long indexStartThisQuantile = (int) (q * pointsInImage / numberOfQuantiles);
                    long indexStartNextQuantile = (int) (((q + 1) * pointsInImage) / numberOfQuantiles);

                    long pointsInQuantile = indexStartNextQuantile - indexStartThisQuantile;

                    // If this is the last quantile, make sure we actually
                    // include everything...
                    if (q == numberOfQuantiles - 1) {
                            indexStartNextQuantile = pointsInImage;
                    }

                    // Keep track of the sum of the values
                    long cumulativeIncluding = 0;
                    long cumulativeBefore = 0;

                    resultSumValuesInQuantile[q] = 0;
                    resultNumberOfValuesInQuantile[q] = 0;

                    for (int value = 0; value < frequencies.length; ++value) {

                            cumulativeIncluding += frequencies[value];

                            if ((cumulativeIncluding < indexStartThisQuantile) || (cumulativeBefore >= indexStartNextQuantile)) {

                                    // Then there's no overlap...

                            } else {

                                    long startInValues = 0;

                                    if (indexStartThisQuantile > cumulativeBefore) {
                                            startInValues = indexStartThisQuantile - cumulativeBefore;
                                    }

                                    // This is the end inclusive...
                                    long endInValues = frequencies[value] - 1;

                                    if (indexStartNextQuantile < cumulativeIncluding) {
                                            endInValues = (indexStartNextQuantile - cumulativeBefore) - 1;
                                    }
                                    long pointsInOverlap = (endInValues - startInValues) + 1;
                                    // System.out.println("points in overlap: "+pointsInOverlap);
                                    resultNumberOfValuesInQuantile[q] += pointsInOverlap;
                                    resultSumValuesInQuantile[q] += value * pointsInOverlap;
                            }

                            cumulativeBefore += frequencies[value];
                    }
            }
    }

    public void generateReplacements(
            ImagePlus imagePlus,
            int numberOfQuantiles,
            long pointsInImage,
            long [] frequencies,
            long [] sumValuesInQuantile,
            long [] numberOfValuesInQuantile,
            double [] quantileMeans,
            Replacements [] resultRankReplacements,
            Replacements [] resultMeanReplacements) {

            int possibleImageValues = frequencies.length;

            for (int q = 0; q < numberOfQuantiles; ++q) {

                    long [] replacementsInThisQuantile=new long[possibleImageValues];

                    long indexStartThisQuantile = (int) (q * pointsInImage / numberOfQuantiles);
                    long indexStartNextQuantile = (int) (((q + 1) * pointsInImage) / numberOfQuantiles);

                    long pointsInQuantile = indexStartNextQuantile - indexStartThisQuantile;

                    // If this is the last quantile, make sure we actually
                    // include everything...
                    if (q == numberOfQuantiles - 1) {
                            indexStartNextQuantile = pointsInImage;
                    }

                    // Keep track of the sum of the values
                    long cumulativeIncluding = 0;
                    long cumulativeBefore = 0;

                    for (int value = 0; value < frequencies.length; ++value) {

                            cumulativeIncluding += frequencies[value];

                            if ((cumulativeIncluding < indexStartThisQuantile) || (cumulativeBefore >= indexStartNextQuantile)) {

                                    // Then there's no overlap...

                            } else {

                                    long startInValues = 0;

                                    if (indexStartThisQuantile > cumulativeBefore) {
                                            startInValues = indexStartThisQuantile - cumulativeBefore;
                                    }

                                    // This is the end inclusive...
                                    long endInValues = frequencies[value] - 1;

                                    if (indexStartNextQuantile < cumulativeIncluding) {
                                            endInValues = (indexStartNextQuantile - cumulativeBefore) - 1;
                                    }
                                    long pointsInOverlap = (endInValues - startInValues) + 1;
                                    numberOfValuesInQuantile[q] += pointsInOverlap;
                                    sumValuesInQuantile[q] += value * pointsInOverlap;
                                    replacementsInThisQuantile[value] = pointsInOverlap;
                            }

                            cumulativeBefore += frequencies[value];
                    }

                    double mean = quantileMeans[q];

                    int byteLowerThanMean = (int) Math.floor(mean);
                    int byteHigherThanMean = (int) Math.ceil(mean);

                    double proportionLower = Math.ceil(mean) - mean;
                    int lowerBytes = (int) Math.round(proportionLower*(indexStartNextQuantile-indexStartThisQuantile));
                    int higherBytes = (int) (numberOfValuesInQuantile[q] - lowerBytes);

                    long replacementsAddedAlready = 0;

                    for( int i = 0; i < possibleImageValues; ++i ) {

                            long r = replacementsInThisQuantile[i];

                            if( r == 0 )
                                    continue;

                            long howManyLowerToAdd = 0;
                            long howManyHigherToAdd = 0;

                            if( replacementsAddedAlready >= lowerBytes ) {
                                    howManyHigherToAdd = r;
                            } else if( replacementsAddedAlready + r >= lowerBytes ) {
                                    howManyLowerToAdd = lowerBytes - replacementsAddedAlready;
                                    howManyHigherToAdd = r - howManyLowerToAdd;
                            } else {
                                    howManyLowerToAdd = r;
                            }

                            resultMeanReplacements[i].addSomeReplacements(howManyLowerToAdd, byteLowerThanMean);
                            resultMeanReplacements[i].addSomeReplacements(howManyHigherToAdd, byteHigherThanMean);

                            resultRankReplacements[i].addSomeReplacements(r, q);

                            replacementsAddedAlready += r;
                    }
            }
    }

    ImagePlus remapImage(ImagePlus imagePlus,
                         int numberOfQuantiles,
                         boolean replaceWithRankInstead,
                         boolean rescaleRanks,
                         Replacements [] rankReplacements,
                         Replacements [] meanReplacements) {

            int originalImageType = imagePlus.getType();
            int width = imagePlus.getWidth();
            int height = imagePlus.getHeight();
            ImageStack stack = imagePlus.getStack();
            int depth = stack.getSize();
            ImageStack newStack = new ImageStack(width,height);
            for( int z = 0; z < depth; ++z ) {
                    byte [] oldPixelsByte = null, newPixelsByte = null;
                    short [] oldPixelsShort = null, newPixelsShort = null;
                    if (originalImageType == ImagePlus.GRAY16) {
                            oldPixelsShort = (short[])stack.getPixels(z+1);
                            newPixelsShort = new short[width*height];
                    } else {
                            oldPixelsByte = (byte[])stack.getPixels(z+1);
                            newPixelsByte = new byte[width*height];
                    }
                    for( int y = 0; y < height; ++y )
                            for( int x = 0; x < width; ++x ) {
                                    int oldValue;
                                    if (originalImageType == ImagePlus.GRAY16)
                                            oldValue = oldPixelsShort[y*width+x]&0xFFFF;
                                    else
                                            oldValue = oldPixelsByte[y*width+x]&0xFF;
                                    int replacement;
                                    if( replaceWithRankInstead ) {
                                            replacement = rankReplacements[oldValue].getRandomReplacement();
                                            if(rescaleRanks)
                                                    replacement = (255*replacement) / (numberOfQuantiles - 1);
                                    } else {
                                            replacement = meanReplacements[oldValue].getRandomReplacement();
                                    }
                                    if( replacement < 0 ) {
                                            System.out.println("BUG: ran out of replacements for "+oldValue);
                                            replacement = oldValue;
                                    }
                                    if (originalImageType == ImagePlus.GRAY16)
                                            newPixelsShort[y*width+x] = (short)replacement;
                                    else
                                            newPixelsByte[y*width+x] = (byte)replacement;
                            }
                    if (originalImageType == ImagePlus.GRAY16) {
                            ShortProcessor sp=new ShortProcessor(width,height);
                            sp.setPixels(newPixelsShort);
                            newStack.addSlice("",sp);
                    } else {
                            ByteProcessor bp=new ByteProcessor(width,height);
                            bp.setPixels(newPixelsByte);
                            newStack.addSlice("",bp);
                    }

                    IJ.showProgress( z / (double)depth );
            }

            IJ.showProgress(1.0);

            if( ImagePlus.COLOR_256 == imagePlus.getType() ) {
                    ColorModel cm = null;
                    cm = stack.getColorModel();
                    if( cm != null ) {
                            newStack.setColorModel( cm );
                    }
            }

            ImagePlus newImage = new ImagePlus( "normalized "+imagePlus.getTitle(), newStack );
            newImage.setCalibration(imagePlus.getCalibration());

            return newImage;
    }

    static final int POSSIBLE_8_BIT_VALUES = 256;
    static final int POSSIBLE_16_BIT_VALUES = 65536;

    public void processToDirectory( FileGroup fg,
                                    String outputDirectory,
                                    int channelToUse,
                                    int numberOfQuantiles,
                                    boolean replaceWithRankInstead,
                                    boolean rescaleRanks ) {

            File o=new File(outputDirectory);
            if( ! o.exists() ) {
                    IJ.error("The output directory ('"+outputDirectory+"') doesn't exist.");
                    return;
            }
            if( ! o.isDirectory() ) {
                    IJ.error("'"+outputDirectory+"' is not a directory");
                    return;
            }

            int n = fg.size();
            if (n < 1) {
                    IJ.error("No image files selected");
                    return;
            }

            /* First go through each image building totalling the
               frequencies of each value. */

            long frequencies[][] = null;
            long pointsInImage[] = new long[n];

            long [][] sumValuesInQuantile = new long[n][numberOfQuantiles];
            long [][] numberOfValuesInQuantile = new long[n][numberOfQuantiles];

            int possibleImageValues = -1;
            System.out.println("Nb file "+n);
            for (int b = 0; b < n; ++b) {

                    File f = fg.get(b);
                    String path = f.getAbsolutePath();
                    System.out.println("Opening file "+path);
                    ImagePlus [] channels=BatchOpener.open(path);

                    if( channelToUse >= channels.length ) {
                            IJ.error("There is no channel "+channelToUse+" in "+path);
                            return;
                    }

                    ImagePlus imagePlus=channels[channelToUse];

                    int type=imagePlus.getType();
                    if( (type == ImagePlus.GRAY8) || (type == ImagePlus.COLOR_256) ) {
                            if (possibleImageValues > 0) {
                                    if (possibleImageValues != POSSIBLE_8_BIT_VALUES) {
                                            IJ.error("The image '"+path+"' was 8 bit, but the previous images were all 16 bit");
                                    }
                            } else {
                                    possibleImageValues = POSSIBLE_8_BIT_VALUES;
                                    frequencies = new long[n][possibleImageValues];
                            }
                    } else if( type == ImagePlus.GRAY16 ) {
                            if (possibleImageValues > 0) {
                                    if (possibleImageValues != POSSIBLE_16_BIT_VALUES) {
                                            IJ.error("The image '"+path+"' was 16 bit, but the previous images were all 8 bit");
                                    }
                            } else {
                                    possibleImageValues = POSSIBLE_16_BIT_VALUES;
                                    frequencies = new long[n][possibleImageValues];
                            }
                    } else {
                            IJ.error("Error processing '"+path+"': This plugin only works on 8bit (GRAY8 or COLOR_256) images.");
                            return;
                    }

                    String freeMemory = IJ.freeMemory();
                    //System.out.println("free memory is: "+freeMemory);

                    int width=imagePlus.getWidth();
                    int height=imagePlus.getHeight();
                    int depth=imagePlus.getStackSize();

                    ImageStack stack=imagePlus.getStack();

                    IJ.showStatus("Calculating frequencies and quantiles for "+imagePlus.getShortTitle()+" ...");

                    for( int z=0; z<depth; ++z ) {
                            if (possibleImageValues == POSSIBLE_8_BIT_VALUES) {
                                    byte [] pixels=(byte[])stack.getPixels(z+1);
                                    for( int y=0; y<height; ++y )
                                            for( int x=0; x<width; ++x ) {
                                                            int value=pixels[y*width+x]&0xFF;
                                                            ++frequencies[b][value];
                                            }
                            } else {
                                    short [] pixels=(short[])stack.getPixels(z+1);
                                    for( int y=0; y<height; ++y )
                                            for( int x=0; x<width; ++x ) {

                                                            int value=pixels[y*width+x]&0xFFFF;
                                                            ++frequencies[b][value];

                                            }
                            }
                    }

                    pointsInImage[b]= width*height*depth;
                    //System.out.println("Proportion of points to consider: "+((double)pointsInImage[b]/(width*height*depth)));

                    divideIntoQuantiles(numberOfQuantiles,
                                        frequencies[b],
                                        pointsInImage[b],
                                        sumValuesInQuantile[b],
                                        numberOfValuesInQuantile[b]);

                    imagePlus.close();
            }

            System.out.println("Now going on to calculate the mean in each quantile.");

            // Calculate the mean in each quantile (even if we're
            // not going to use it)...

            double [] quantileMeans = new double[numberOfQuantiles];

            for( int q = 0; q < numberOfQuantiles; ++q ) {
                    long sum = 0;
                    long values = 0;
                    for( int b = 0; b < n; ++ b ) {
                            sum += sumValuesInQuantile[b][q];
                            values += numberOfValuesInQuantile[b][q];
                    }
                    quantileMeans[q] = sum / (double)values;
            }

            // Now we go through each image again, remap the
            // values according to the options chosen and write
            // the new image out to the output directory....

            for (int b = 0; b < n; ++b) {

                    File f = fg.get(b);
                    String path = f.getAbsolutePath();
                    //System.out.println("Opening file "+path);
                    ImagePlus [] channels=BatchOpener.open(path);
                    ImagePlus imagePlus=channels[channelToUse];

                    String newLeafName;
                    String leafName = f.getName();
                    int dotIndex=leafName.lastIndexOf(".");
                    if(dotIndex >= 0) {
                            newLeafName = leafName.substring(0,dotIndex) + "-normalized.tif";
                    } else {
                            newLeafName = leafName + "-normalized";
                    }

                    File outputFile=new File(outputDirectory,newLeafName);

                    /* meanReplacements or rankReplacements are
                       the arrays that are ultimately used to get
                       replacement values for the image, as in:

                         rankReplacements[oldValue].getRandomReplacement()

                         The Replacment class stores a particular
                         number of replacment values, and tracks
                         how many are left after removing a random
                         one.
                    */

                    Replacements [] meanReplacements = new Replacements[possibleImageValues];
                    for( int value = 0; value < possibleImageValues; ++value )
                            meanReplacements[value] = new Replacements(possibleImageValues);

                    Replacements [] rankReplacements = new Replacements[possibleImageValues];
                    for( int value = 0; value < possibleImageValues; ++value )
                            rankReplacements[value] = new Replacements(numberOfQuantiles);

                    int width=imagePlus.getWidth();
                    int height=imagePlus.getHeight();
                    int depth=imagePlus.getStackSize();

                    IJ.showStatus("Replacing values in: "+imagePlus.getShortTitle()+" ...");

                    generateReplacements(imagePlus,
                                         numberOfQuantiles,
                                         pointsInImage[b],
                                         frequencies[b],
                                         sumValuesInQuantile[b],
                                         numberOfValuesInQuantile[b],
                                         quantileMeans,
                                         rankReplacements,
                                         meanReplacements);

                    IJ.showProgress(0);

                    ImagePlus newImage = remapImage(imagePlus,
                                                    numberOfQuantiles,
                                                    replaceWithRankInstead,
                                                    rescaleRanks,
                                                    rankReplacements,
                                                    meanReplacements);

                    // newImage.show();

                    boolean saved;
                    if (newImage.getStackSize() == 1)
                            saved = new FileSaver(newImage).saveAsTiff(outputFile.getAbsolutePath());
                    else
                            saved = new FileSaver(newImage).saveAsTiffStack(outputFile.getAbsolutePath());
                    if( ! saved )
                            return;

                    newImage.close();
                    imagePlus.close();

            }

            IJ.showStatus("Normalization complete: files written to: "+outputDirectory);

    }


    public void run(String dir, ArrayList<String> listFiles, String channel) {
            FileGroup fg = new FileGroup("foo");
            for( String f : listFiles) {
                String fileName = dir + FilenameUtils.getBaseName(f)+channel+".tif";
                System.out.println(fileName);
                fg.add(fileName);
            }
            int channelToUse = 0;	
            int numberOfQuantiles = 256;
            boolean replaceWithRankInstead=false;
            boolean rescaleRanks=true;

            processToDirectory( fg,
                                dir,
                                channelToUse,
                                numberOfQuantiles,
                                replaceWithRankInstead,
                                rescaleRanks );


    }
	
}
