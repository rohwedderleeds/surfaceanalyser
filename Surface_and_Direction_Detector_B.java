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
        //copy image and prepare threshold image
        stackip = stack2.getProcessor(i);

        stackip.invert();

        //-----------------
        imp2.setSlice(i);

        //-----------------
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
        //---------------------------

        //---------------------------
        party1.analyze(imp2);

        int howmany = manager.getCount();
        //IJ.log("slice "+i+" step 1 "+ howmany);

        //tempResults = manager.multiMeasure(imp2);
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
        Roi roix = manager.getRoi(0);

        first.setSlice(i);
        firstip.setColor(255);
        //---------------------------

        roix.drawPixels(firstip);
        //firstip.fill(roix);

        double centerx = tempResults.getValue("XM",0);
        double centery = tempResults.getValue("YM",0);
        firstip.setColor(60);
        firstip.drawOval((int) centerx, (int) centery,10,10);
        centerx = centerx / pixres;
        centery = centery / pixres;

        //---------------------------
        FloatPolygon dings = new FloatPolygon();
        dings = roix.getFloatPolygon();
        lang = dings.npoints;
        extralang = 2*lang;
        double large = (lang/100)*newdirect;
        double small = (lang/100)*newfilos;
        //IJ.log("All "+lang+" large "+large+" small "+small);


        // get radius (distance to centre)
        double radiusnew [] = new double[lang];
        double radiusextra [] = new double[extralang];
        double radtrans = 0;
        int x1[] = new int[lang];
        int y1[] = new int[lang];

        //IJ.log(""+centerx+" "+centery);//Correct

        tempResults.reset();
        //IJ.log("slice "+i+" step2 ");

        for (int test=0; test < lang; test++){
            x1[test] = (int) Math.round(dings.xpoints[test]);
            y1[test] = (int) Math.round(dings.ypoints[test]);
            //IJ.log("slice "+i+" x "+x1[test]+ " y "+ y1[test]); Correct
            radtrans = Math.pow(Math.abs(x1[test]-centerx),2)+Math.pow(Math.abs(y1[test]-centery),2);
            radiusnew[test] = Math.sqrt(radtrans);
            //IJ.log(""+radiusnew[test]);//Correct
        }
        int testzahl = 0;
        for (int test=0; test < extralang; test++){
                radiusextra[test] = radiusnew[testzahl];
                testzahl++;
                if (testzahl > lang-1){testzahl = testzahl - lang;}

        }
        //get rounded radius
        double radiusround [] = new double[lang];
        double sumdist = 0; //sumcount,
        int correctpos = 0;
        for (int test=0; test < lang; test++){
                for (int i2=0;i2 < (lang/3);i2++){
                    correctpos = test+i2;
                    if (correctpos>= lang) {correctpos = correctpos-lang;}
                    sumdist = sumdist + radiusnew[correctpos];
                }
                radiusround[test] = sumdist/(lang/3);
                //IJ.log("slice "+i+" radius "+radiusround[test]);correct
                sumdist = 0;
        }

        //get rounded positions
        double xround [] = new double[lang];
        int xroundint [] = new int[lang];
        double yround [] = new double[lang];
        int yroundint [] = new int[lang];
        double angle;
        int roundsumx = 0; int roundsumy = 0;

        //----------------------------------------------rundrund
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



        //-------------------------------------------rundrundende

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
                //firstip.setColor(180);
                //firstip.drawPixel((int) xround[test],(int) yround[test]);
                roundsumx = roundsumx + xroundint[test];
                roundsumy = roundsumy + yroundint[test];
        }
        manager.reset();

        //double roundcenterx = roundsumx/lang;
        double roundcenterx = Broundsumx/lang;
        //double roundcentery = roundsumy/lang;
        double roundcentery = Broundsumy/lang;

        //---------------------------- wieder aktivieren
        IJ.log("Slice: "+i+" Horizontal distance of outer/inner center of mass: "+((centerx-roundcenterx)*pixres)+" Vertical distance of outer/inner center of mass: "+((centery-roundcentery)*pixres));
        //---------------------------- wieder aktivieren
        //calculate major protrusions

        //------------------------------------------------------------------------------------------------------

    int slopecount = 0;
    int slopestart = 0;
    int slopelength = 0;
    double slope [] = new double[lang+1];
    slope = slopecalculation(large,extralang,radiusextra);
    //IJ.log("slice "+i+" step4 ");
    double revradius [] = new double[lang+1];
    revradius = reverser(extralang, radiusextra);
    //IJ.log("slice "+i+" step5 ");
    double revslopeA [] = new double[lang+1];
    revslopeA = slopecalculation(large,extralang,revradius);
    //IJ.log("slice "+i+" step6 ");
    double revslope [] = new double[lang+1];
    revslope =  reverser(extralang, revslopeA);
    //IJ.log("slope:"+slope[test]+" revslop:"+revslope[test]+" radius:"+radiusnew[test]);
    //IJ.log(""+radiusnew[test])
    //IJ.log("slice "+i+" step6 ");
    //for (int clear = 0; clear < lang; clear++){IJ.log(""+clear+" "+radiusnew[clear]+" "+slope[clear]+" "+revslopeA[clear]);}
    //for (int test=0; test < extralang; test++){IJ.log(""+radiusextra[test]+" "+slope[test]+" "+revslope[test]+" "+i);}

//calculate major protrusions

    stackcount = 0;
    stackcount = (i-1)*4;
    int helper = 1;
//show major protrusions

        for (int major = 1; major < lang; major++){
                if (Math.abs(slope[helper]/revslope[helper])>sensitivity){
                    slopecount++;
                    slopestart = helper;
                    //while (slope[major]/revslope[major]>100) {major++;}
                    //if (slope[major]<0){while (slope[major]<0){major++;}
                    if (slope[helper]<0){
                            helper++;
                            if (helper > lang){helper = helper - lang;}
                        //IJ.log("major "+major+" "+i);
                        while (slope[helper]<0){

                            helper++;
                            }
                        }
                        if (helper > lang){helper = helper - lang;}
                    while (slope[helper+1]*revslope[helper+1]>slope[helper]*revslope[helper]) {helper++;if (helper > lang){helper = helper - lang;}}//IJ.log("major "+helper+" "+i);
                    rt1.incrementCounter();
                    slopelength = helper - slopestart;
                    rt1.addValue("X position slice "+i, x1[(int) Math.round(slopestart+(slopelength))]);//(int) Math.round(slopestart+(slopelength/2))
                    rt1.addValue("Y position slice "+i, y1[(int) Math.round(slopestart+(slopelength))]);//(int) Math.round(slopestart+(slopelength/2))
                    rt1.addValue("Broad "+i, slopelength*pixres);
                    rt1.addValue("Length "+i, radiusnew[(int) Math.round(slopestart+(slopelength))]*pixres);//(int) Math.round(slopestart+(slopelength/2))
                    stackip.setColor(0);
                    firstip.setColor(130);
                    firstip.setLineWidth(2);
                    //firstip.drawLine(((int) Math.round(centerx)),((int) Math.round(centery)),x1[(int) Math.round(slopestart+(slopelength/2))],y1[(int) Math.round(slopestart+(slopelength/2))]);//(int) Math.round(slopestart+(slopelength/2))
                    //firstip.drawLine(((int) Math.round(centerx)),((int) Math.round(centery)),x1[(int) Math.round(slopestart+(slopelength/2))],y1[(int) Math.round(slopestart+(slopelength/2))]);//(int) Math.round(slopestart+(slopelength/2))
                    firstip.drawLine(((int) Math.round(roundcenterx)),((int) Math.round(roundcentery)),x1[(int) Math.round(slopestart+(slopelength))],y1[(int) Math.round(slopestart+(slopelength))]);//(int) Math.round(slopestart+(slopelength/2))
                    //while (slope[major]/revslope[major]>1000) {major++;}
                    while (slope[helper+1]*revslope[helper+1]>slope[helper]*revslope[helper]) {helper++;}

                }
                helper++;
        }


        //---------------------------- wieder aktivieren
        //IJ.log(" Major Protrusions: "+slopecount);
        //---------------------------- wieder aktivieren
        int slopecountlarge = slopecount;


        slopecount = 0;
        for (int clear = 0; clear < lang; clear++){slope[clear] = 0;revslopeA[clear] = 0;}
        //calculate minor protrusions

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
        //---------------------------- wieder aktivieren
        IJ.log("Major Protrusions: "+slopecountlarge+" Filopodia: "+slopecount);
        //---------------------------- wieder aktivieren
        for (int clear = 0; clear < lang; clear++){slope[clear] = 0;revslopeA[clear] = 0;}

    manager.close();

    }
    //write results tables
    rt1.show("Major protrusions");
    rt2.show("Filopodia");
    stackip.invert();
    //ImagePlus outstackbild = new ImagePlus(title,stackip);
    //outstackbild.setCalibration(imp1.getCalibration());
    //outstackbild.show();
    }


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
        //IJ.log("major "+major+" step3A ");
        for (int direction = 0; direction < direct; direction ++){

            poscorrect = major + direction;
            if (poscorrect > lang-1) {poscorrect = poscorrect - lang;}
            //IJ.log("poscorrect "+poscorrect);
            sumxy = (direction*radius[poscorrect])+sumxy;
            sumxint = direction + sumxint;
            sumxsquare = direction * direction;
            sumyint = radius[poscorrect]+sumyint;
        }
    slopeint[major] = (((direct)*sumxy)-(sumxint*sumyint))/(((direct)*sumxsquare)-(sumxint*sumxint));
    sumxy = 0; sumxsquare = 0; sumxint = 0; sumyint = 0;
    //IJ.log("slope "+slopeint[major]+" step3A ");
    }
    return slopeint;
}

}
