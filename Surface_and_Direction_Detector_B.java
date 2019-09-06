import ij.*;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.measure.*;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.*;
import ij.plugin.filter.RankFilters;
import ij.plugin.frame.PlugInFrame;
import ij.plugin.frame.RoiManager;



/** The Surface_Observer plugin generates outlines ignoring enclosed structures and considering the ROI positions as a function of x to evaluate deviation from linear regression. Works for grey level images. */
public class Surface_and_Direction_Detector_B implements PlugIn {

    static double newdirect = 30;
    static double newfilos = 1;
    static double newcover = 30;
    static double sensitivity = 150;
    static double sumx, sumy;

    public void run(String arg) {

        if (IJ.getImage().getBitDepth()!=8){IJ.showMessage("Grey scale measure", "Grey scale Image required");}
        String directory =   IJ.getDir("plugins");
        directory = directory + "Analyze/sketch.jpg";
        directory = directory.replace('\\', '/');
        int lang = 0;
        int extralang = 0;
        int stackcount = 0;

        ImagePlus imp1 = IJ.getImage();
        ImagePlus imp2 = imp1.duplicate();

        imp2.setCalibration(imp1.getCalibration());
        ImageStack stack = imp1.getStack();
        ImageStack stack2 = imp2.getStack();
        ImagePlus Sketch = IJ.openImage (directory);
        int size = stack.getSize();

        String title = imp1.getTitle();
        Calibration cal = imp1.getCalibration();
        ImageProcessor ip = (ImageProcessor)imp1.getProcessor();
        ImageProcessor ip2 = (ImageProcessor)imp2.getProcessor();
        double pixres = cal.pixelWidth;

//Dialog
        GenericDialog gd = new GenericDialog("Membrane cover", IJ.getInstance());
        gd.addMessage("Membrane cover describes the expected size proportion for a structure in stretch of membrane.");
        gd.addImage(Sketch);
        gd.addMessage("Percentage should be larger (>= 30 percent and larger) for large protrusion structures.");
        gd.addMessage("Percentage should be smaller (~ 2 percent) for small filopodia structures.");
        gd.addNumericField("Major protrusions percentage: ", newdirect, 3);
        gd.addNumericField("Filopodia percentage: ", newfilos, 3);
        gd.addNumericField("Image coverage of cell: (minimum percentage of whole image)", newcover, 3);
        gd.addNumericField("Sensitivity of measurement", sensitivity, 3);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        if (gd.invalidNumber()) {
            IJ.showMessage("Error", "Threshold not in range");
            return;
        }
        newdirect = gd.getNextNumber();
        newfilos = gd.getNextNumber();
        newcover = gd.getNextNumber();
        sensitivity = gd.getNextNumber();

//output images setup

        int w = ip.getWidth();
        int h = ip.getHeight();
        ImageStack outstack = stack.duplicate();
        ImageProcessor outstackip = (ImageProcessor)imp1.getProcessor();
        ImagePlus first = IJ.createImage(title+"Outline", "8-bit black", w, h, size);// from automatic selection
        ImageProcessor firstip = (ImageProcessor)first.getProcessor();
        ImageProcessor stackip = (ImageProcessor)imp1.getProcessor();
        ResultsTable rt1 = new ResultsTable();
        ResultsTable rt2 = new ResultsTable();
        first.show();

        for (int i=1;i<size+1;i++){

//prepare image

            stackip = stack2.getProcessor(i);
            stackip.invert();
            imp2.setSlice(i);
            IJ.setAutoThreshold(imp2, "Huang");
            RankFilters desp = new RankFilters();
            desp.rank(stackip, 1.0, desp.MEDIAN,desp.BRIGHT_OUTLIERS,0);

//get outline

            RoiManager manager = RoiManager.getInstance();
            if (manager == null)
            {
                manager = new RoiManager();
            }
            ResultsTable tempResults = new ResultsTable();
            ParticleAnalyzer party1 = new ParticleAnalyzer( ParticleAnalyzer.ADD_TO_MANAGER + ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES + ParticleAnalyzer.INCLUDE_HOLES,Measurements.CENTER_OF_MASS,tempResults,((w*h/100)*newcover),((w-50)*(h-50)));

//use particle analyser

            party1.analyze(imp2);
            int howmany = manager.getCount();
            double covercorrect = newcover;
            if (howmany < 1)
            {
                while (howmany<1 && covercorrect>0)
                {
                    ParticleAnalyzer party2 = new ParticleAnalyzer( ParticleAnalyzer.ADD_TO_MANAGER + ParticleAnalyzer.INCLUDE_HOLES,Measurements.CENTER_OF_MASS,tempResults,((w*h/100)*covercorrect),((w-50)*(h-50)));
                    party2.analyze(imp2);
                    howmany = manager.getCount();
                    covercorrect =  covercorrect-0.5;
                }
            }
            party1.analyze(imp2);
            if (howmany<1 && covercorrect<0.1){continue;}

//extract ROI

            Roi roix = manager.getRoi(0);
            first.setSlice(i);
            firstip.setColor(255);
            roix.drawPixels(firstip);

//get centre of mass

            double centerx = tempResults.getValue("XM",0);
            double centery = tempResults.getValue("YM",0);
            firstip.setColor(60);
            firstip.drawOval((int) centerx, (int) centery,10,10);
            centerx = centerx / pixres;
            centery = centery / pixres;

//transform ROI and calculate windows

            FloatPolygon dings = new FloatPolygon();
            dings = roix.getFloatPolygon();
            lang = dings.npoints;
            extralang = 2*lang;
            double large = (lang/100)*newdirect;
            double small = (lang/100)*newfilos;


//get radius (distance to centre)

            double radiusnew [] = new double[lang];
            double radiusextra [] = new double[extralang];
            double radtrans = 0;
            int x1[] = new int[lang];
            int y1[] = new int[lang];
            tempResults.reset();

            for (int test=0; test < lang; test++){
                x1[test] = (int) Math.round(dings.xpoints[test]);
                y1[test] = (int) Math.round(dings.ypoints[test]);
                radtrans = Math.pow(Math.abs(x1[test]-centerx),2)+Math.pow(Math.abs(y1[test]-centery),2);
                radiusnew[test] = Math.sqrt(radtrans);
            }
            int testzahl = 0;
            for (int test=0; test < extralang; test++){
                    radiusextra[test] = radiusnew[testzahl];
                    testzahl++;
                    if (testzahl > lang-1){testzahl = testzahl - lang;}
            }

//get rounded radius

            double radiusround [] = new double[lang];
            double sumdist = 0;
            int correctpos = 0;
            for (int test=0; test < lang; test++){
                    for (int i2=0;i2 < (lang/3);i2++){
                        correctpos = test+i2;
                        if (correctpos>= lang) {correctpos = correctpos-lang;}
                        sumdist = sumdist + radiusnew[correctpos];
                    }
                    radiusround[test] = sumdist/(lang/3);
                    sumdist = 0;
            }

//get rounded positions

            double xround [] = new double[lang];
            int xroundint [] = new int[lang];
            double yround [] = new double[lang];
            int yroundint [] = new int[lang];
            double angle;
            int roundsumx = 0; int roundsumy = 0;
            int x1lang [] = new int[extralang];
            int y1lang [] = new int[extralang];
            int intcounter = 0;
            for (int fill = 0; fill < extralang; fill++){
                x1lang [fill] = x1[intcounter];
                y1lang [fill] = y1[intcounter];
                intcounter++;
                if (intcounter > lang-1){intcounter = intcounter - lang;}
            }

            int Broundsumx = 0; int Broundsumy = 0;
            firstip.setColor(180);
            for (int rund = 0; rund < lang; rund++){
                for (int limit = 0; limit < (lang/100*50); limit ++){
                    roundsumx = roundsumx+x1lang[limit+rund];
                    roundsumy = roundsumy+y1lang[limit+rund];
                }
                xround[rund]=roundsumx/(lang/100*50);
                yround[rund]=roundsumy/(lang/100*50);
                firstip.drawPixel((int) xround[rund],(int) yround[rund]);
                roundsumx=0;roundsumy=0;
                Broundsumx = Broundsumx + (int) xround[rund];Broundsumy = Broundsumy + (int) yround[rund];
            }

            roundsumx=0;roundsumy=0;

            for (int test=0; test < lang; test++){
                    angle = Math.acos((Math.abs(centery-dings.ypoints[test]))/radiusnew[test]);
                    if ((centery-dings.ypoints[test])>0){
                        yround[test]=centery+(radiusround[test]*Math.cos(Math.toRadians(angle)));
                    }
                    if ((centerx-dings.xpoints[test])>0){
                        xround[test]=centerx+(radiusround[test]*Math.sin(Math.toRadians(angle)));
                    }
                    if ((centery-dings.ypoints[test])<0){
                        yround[test]=centery-(radiusround[test]*Math.cos(Math.toRadians(angle)));
                    }
                    if ((centerx-dings.xpoints[test])<0){
                        xround[test]=centerx+(radiusround[test]*Math.sin(Math.toRadians(angle)));
                    }
                    yroundint[test] = (int) Math.round(yround[test]);
                    xroundint[test] = (int) Math.round(xround[test]);
                    roundsumx = roundsumx + xroundint[test];
                    roundsumy = roundsumy + yroundint[test];
            }
            manager.reset();

            double roundcenterx = Broundsumx/lang;
            double roundcentery = Broundsumy/lang;

            IJ.log("Slice: "+i+" Horizontal distance of outer/inner center of mass: "+((centerx-roundcenterx)*pixres)+" Vertical distance of outer/inner center of mass: "+((centery-roundcentery)*pixres));

//calculate major number of protrusions

            int slopecount = 0;
            int slopestart = 0;
            int slopelength = 0;
            double slope [] = new double[lang+1];
            slope = slopecalculation(large,extralang,radiusextra);
            double revradius [] = new double[lang+1];
            revradius = reverser(extralang, radiusextra);
            double revslopeA [] = new double[lang+1];
            revslopeA = slopecalculation(large,extralang,revradius);
            double revslope [] = new double[lang+1];
            revslope =  reverser(extralang, revslopeA);

            stackcount = 0;
            stackcount = (i-1)*4;
            int helper = 1;
            int protcount = 0;
            int schalter = 0;

            for (int major = 1; major < lang; major++){
                    if (Math.abs(slope[helper]/revslope[helper])>sensitivity){
                        if (slope[helper]<0){
                            helper++;
                            if (helper > lang){helper = helper - lang;}
                            protcount++;
                            while (slope[helper]<0){
                                helper++;
                            }
                        }
                        if (helper > lang){helper = helper - lang;}
                        while (slope[helper+1]*revslope[helper+1]>slope[helper]*revslope[helper]) {helper++;if (helper > lang){helper = helper - lang;}}
                        while (slope[helper+1]*revslope[helper+1]>slope[helper]*revslope[helper]) {helper++;}
                    }
                    helper++;
            }

//generate major protrusion table

            stackip.setColor(0);
            firstip.setColor(130);
            firstip.setLineWidth(2);
            for (int major = 1; major < lang; major++){
                if (Math.abs(slope[helper]/revslope[helper])>sensitivity){
                    slopecount++;
                    slopestart = helper;
                    if (slope[helper]<0){
                        helper++;
                        if (helper > lang){helper = helper - lang;}
                        while (slope[helper]<0){
                            helper++;
                            }
                        }
                        if (helper > lang){helper = helper - lang;}
                        while (slope[helper+1]*revslope[helper+1]>slope[helper]*revslope[helper]) {helper++;if (helper > lang){helper = helper - lang;}}
                        slopelength = helper - slopestart;
                        if (schalter < (protcount)){
                            rt1.incrementCounter();
                            rt1.addValue("X position slice "+i, x1[(int) Math.round(slopestart+(slopelength))]);
                            rt1.addValue("Y position slice "+i, y1[(int) Math.round(slopestart+(slopelength))]);
                            rt1.addValue("Broad "+i, slopelength*pixres);
                            rt1.addValue("Length "+i, radiusnew[(int) Math.round(slopestart+(slopelength))]*pixres);
                        }
                        schalter++;
                        firstip.drawLine(((int) Math.round(roundcenterx)),((int) Math.round(roundcentery)),x1[(int) Math.round(slopestart+(slopelength))],y1[(int) Math.round(slopestart+(slopelength))]);//(int) Math.round(slopestart+(slopelength/2))
                        while (slope[helper+1]*revslope[helper+1]>slope[helper]*revslope[helper]) {helper++;}
                }
                helper++;
            }
            int slopecountlarge = slopecount;
            slopecount = 0;
            for (int clear = 0; clear < lang; clear++){slope[clear] = 0;revslopeA[clear] = 0;}

//calculate minor protrusions and generate table

            slope = slopecalculation(small,lang,radiusnew);
            revslopeA = slopecalculation(small,lang,revradius);
            revslope =  reverser(lang, revslopeA);

            int minorcount = 0;
            for (int minor = 0; minor < lang; minor++){
                if (revslope[minor]+slope[minor]<0){
                    slopecount++;
                    slopestart = minor;
                    while (slope[minor]+revslope[minor]<0 && minor<(lang-1)) {
                         minor++;
                    }
                    slopelength = minor - slopestart;
                    rt2.incrementCounter();
                    rt2.addValue("X position slice "+i, x1[(int) Math.round(slopestart+(slopelength/2))]);
                    rt2.addValue("Y position slice "+i, y1[(int) Math.round(slopestart+(slopelength/2))]);
                    rt2.addValue("Length "+i, slopelength*pixres);
                    rt2.addValue("Broad "+i, radiusnew[(int) Math.round(slopestart+(slopelength/2))]*pixres);
                }
            }

//Output results

            IJ.log("Slice: "+i+" Major Protrusions: "+protcount+" Filopodia: "+slopecount);
            for (int clear = 0; clear < lang; clear++){slope[clear] = 0;revslopeA[clear] = 0;}
            manager.close();
        }

//show results tables

        rt1.show("Major protrusions");
        rt2.show("Filopodia");
        stackip.invert();
        }

//procedures for slope calculation

    public double [] reverser(int lang, double [] array) {
        int negcorrect = lang;
        double revarray [] = new double[lang];
        for (int revers = 0; revers < lang; revers++){
            negcorrect--;
            revarray[revers] = array[negcorrect];
        }
        return revarray;
    }

    public double [] slopecalculation (double direct, int lang, double [] radius){
        double slopeint [] = new double[lang+1];
        double sumxy = 0;
        double sumxsquare = 0;
        double sumxint = 0;
        double sumyint = 0;
        int poscorrect = 0;
        for (int major = 0; major < lang; major++){
            for (int direction = 0; direction < direct; direction ++){
                poscorrect = major + direction;
                if (poscorrect > lang-1) {poscorrect = poscorrect - lang;}
                sumxy = (direction*radius[poscorrect])+sumxy;
                sumxint = direction + sumxint;
                sumxsquare = direction * direction;
                sumyint = radius[poscorrect]+sumyint;
            }
            slopeint[major] = (((direct)*sumxy)-(sumxint*sumyint))/(((direct)*sumxsquare)-(sumxint*sumxint));
            sumxy = 0; sumxsquare = 0; sumxint = 0; sumyint = 0;
        }
        return slopeint;
    }
}
