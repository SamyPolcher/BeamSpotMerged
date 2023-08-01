package analysis;

import java.io.ByteArrayOutputStream;

import java.io.PrintStream;
import objects.Event;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.IDataSet;

import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.graphics.IDataSetPlotter;
import org.jlab.jnp.utils.benchmark.ProgressPrintout;
import org.jlab.utils.options.OptionParser;


import objects.Track;
import objects.Event;
import modules.DCModule;
import modules.CDModule;

public class BeamSpot {
    
    private final boolean debug = false;
    private boolean fastmode = false;
    ByteArrayOutputStream pipeOut = new ByteArrayOutputStream();
    private static PrintStream outStream = System.out;
    private static PrintStream errStream = System.err;
     
    ArrayList<Module>    modules = new ArrayList<>();

    private static String OPTSTAT = "";
    
    public BeamSpot(Boolean runCD, Boolean runDC, String opts, double[] thetaBins, float Ztarget, float zmin, float zmax, int NphiBins) {
        this.init(runCD, runDC, opts, thetaBins, Ztarget, zmin, zmax, NphiBins);
    }
    

    private void init(Boolean runCD, Boolean runDC, String opts, double[] thetaBins, float Ztarget, float zmin, float zmax, int NphiBins) {
        OPTSTAT = opts;
        GStyle.getH1FAttributes().setOptStat(opts);
        GStyle.getAxisAttributesX().setTitleFontSize(24);
        GStyle.getAxisAttributesX().setLabelFontSize(18);
        GStyle.getAxisAttributesY().setTitleFontSize(24);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        GStyle.getAxisAttributesZ().setLabelFontSize(14);
//        GStyle.getAxisAttributesX().setLabelFontName("Arial");
//        GStyle.getAxisAttributesY().setLabelFontName("Arial");
//        GStyle.getAxisAttributesZ().setLabelFontName("Arial");
//        GStyle.getAxisAttributesX().setTitleFontName("Arial");
//        GStyle.getAxisAttributesY().setTitleFontName("Arial");
//        GStyle.getAxisAttributesZ().setTitleFontName("Arial");
        GStyle.setGraphicsFrameLineWidth(2);
        GStyle.getH1FAttributes().setLineWidth(1);
        
        if(runCD){
            CDModule cd = new CDModule();
            cd.init();
            this.modules.add(cd);
        }
        if(runDC){
            DCModule dc = new DCModule();
            dc.setThetaBins(thetaBins);
            dc.setTargetZ(Ztarget, zmin, zmax);
            dc.setBinsPerSector(NphiBins);
            dc.init();
            this.modules.add(dc);
        }
    }
    
    private void processEvent(DataEvent de) {
        Event event = new Event(de);
        for(Module m : modules) m.processEvent(event);
    }

    private void analyzeHistos() {
        for(Module m : modules) m.analyzeHistos();
    }

    public JTabbedPane plotHistos() {
        JTabbedPane panel = new JTabbedPane();
        for(Module m : modules) {
            EmbeddedCanvasTabbed canvas = m.plotHistos();
            for(String name : m.getCanvasNames()) {
                for(EmbeddedPad p : canvas.getCanvas(name).getCanvasPads()) {
                    for(IDataSetPlotter dsp: p.getDatasetPlotters()) {
                        IDataSet ds = dsp.getDataSet();
                        if(ds instanceof H1F) {
                            H1F h1 = (H1F) ds;
                            h1.setOptStat(OPTSTAT);
                        }
                    }
                }            
            }
            panel.add(m.getName(), canvas);
        }
        return panel;
    }
    
    // public void plotDC() {
    //     modules.get(1).plot(true);
    // }
    
    public void readHistos(String fileName) {
        System.out.println("Opening file: " + fileName);
        PrintStream pipeStream = new PrintStream(pipeOut);
        System.setOut(pipeStream);
        System.setErr(pipeStream);
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        System.out.println(dir.getDirectoryList());
        dir.cd();
        dir.pwd();
        for(Module m : modules) {
            m.readDataGroup(dir);
        }
        System.setOut(outStream);
        System.setErr(errStream);
    }

    public void saveHistos(String fileName) {
        System.out.println("\n>>>>> Saving histograms to file " + fileName);
        PrintStream pipeStream = new PrintStream(pipeOut);
        System.setOut(pipeStream);
        System.setErr(pipeStream);
        TDirectory dir = new TDirectory();
        for(Module m : modules) {
            m.writeDataGroup(dir);
        }
        dir.writeFile(fileName);
        System.setOut(outStream);
        System.setErr(errStream);
	    System.out.println("\n>>>>> "+fileName+" written");	
    }

    private void printHistos() {
        System.out.println("\n>>>>> Printing canvases to directory plots");
        for(Module m : modules) {
            m.printHistos("plots");
        }
    }
    
    public void writeCCDB(String outputPrefix) {
        for(Module m : modules) {
            m.writeCCDB(outputPrefix);
        }
    }
    
    private void testHistos() {
        for(Module m : modules) {
            // Not implemented
            m.testHistos();
        }
    }
    
    public static void main(String[] args) {
        
        OptionParser parser = new OptionParser("Beam Spot [options] file1 file2 ... fileN");
        parser.setRequiresInputList(false);

        // general options
        parser.addOption("-o"          ,"",     "histogram file name prefix");
        parser.addOption("-n"          ,"-1",   "maximum number of events to process");
        parser.addOption("-X"          ,"0",    "do NOT save histograms in a hipo file");
        parser.addOption("-CD"         ,"1",    "set to 0 to deactivate CD beamspot analysis");
        parser.addOption("-DC"         ,"1",    "set to 0 to deactivate DC beamspot analysis");
        parser.addOption("-x0"         ,"0",    "x position of the beam used in the reconstruction (average raster position if raster is used)");
        parser.addOption("-y0"         ,"0",    "y position of the beam used in the reconstruction (average raster position if raster is used)");

        // histogram settings
        parser.addOption("-histo"      ,"0",    "read histogram from hipo file (0/1)");
        parser.addOption("-plot"       ,"1",    "display histograms (0/1)");
        parser.addOption("-print"      ,"0",    "print canvases (0/1)");
        parser.addOption("-stats"      ,"",     "histogram stat option (e.g. \"10\" will display entries)");
        
        // DC analysis settigs
        parser.addOption("-Zvertex", "25.4", "nominal Z position of the Foil");
        parser.addOption("-zmin", "0.", "lower bound of the foil fit window in z");
        parser.addOption("-zmax", "0.", "upper bound of the foil fit window in z");
        parser.addOption("-Nphi", "10", "phi bins per sector");

        parser.parse(args);
        
        String namePrefix  = parser.getOption("-o").stringValue();        
        String histoName   = "histo.hipo";
        if(!namePrefix.isEmpty()) {
            histoName  = namePrefix + "_" + histoName;
        }

        int     maxEvents     = parser.getOption("-n").intValue();
        boolean saveHistos    = (parser.getOption("-X").intValue()==0);
        boolean runCD         = (parser.getOption("-CD").intValue()==1);
        boolean runDC         = (parser.getOption("-DC").intValue()==1);
        boolean readHistos    = (parser.getOption("-histo").intValue()!=0);            
        boolean openWindow    = (parser.getOption("-plot").intValue()!=0);
        boolean printHistos   = (parser.getOption("-print").intValue()!=0);
        String  optStats      = parser.getOption("-stats").stringValue(); 
        float  zTarget        = (float)parser.getOption("-Zvertex").doubleValue();
        float  zmin           = (float)parser.getOption("-zmin").doubleValue();
        float  zmax           = (float)parser.getOption("-zmax").doubleValue();
        int NphiBins          = parser.getOption("-Nphi").intValue();
        
        // beam position in cooking
        float  x0             = (float)parser.getOption("-x0").doubleValue();
        float  y0             = (float)parser.getOption("-y0").doubleValue();
        System.out.printf("nominal beam position (%2.3f, %2.3f) cm\n", x0, y0);
        Event.setAvgBeamPos(x0, y0);

        if(!openWindow) System.setProperty("java.awt.headless", "true");
        
        // double[] thetaBins = new double[]{10,20,30};
        // double[] thetaBins = new double[]{10,11,12,13,14,16,18,22,30};
        double[] thetaBins = new double[]{10, 15, 20, 25, 30};
        BeamSpot bs = new BeamSpot(runCD, runDC, optStats, thetaBins, zTarget, zmin, zmax, NphiBins);
        
        List<String> inputList = parser.getInputList();
        if(inputList.isEmpty()==true){
            parser.printUsage();
            System.out.println("\n >>>> error: no input file is specified....\n");
            System.exit(0);
        }

        if(readHistos) {
            for(String inputFile : inputList){
                bs.readHistos(inputFile);
            }
            bs.analyzeHistos();
        }
        else{

            ProgressPrintout progress = new ProgressPrintout();

            int counter = -1;
            for(String inputFile : inputList){
                HipoDataSource reader = new HipoDataSource();
                reader.open(inputFile);
                
                while (reader.hasEvent()) {

                    counter++;

                    DataEvent event = reader.getNextEvent();
                    bs.processEvent(event);
                    
                    progress.updateStatus();
                    if(maxEvents>0){
                        if(counter>=maxEvents) break;
                    }
                }
                progress.showStatus();
                reader.close();
                if(maxEvents>0){
                    if(counter>=maxEvents) break;
                }
            }   
            bs.analyzeHistos();
            if(saveHistos) bs.saveHistos(histoName);
        }
        
        bs.writeCCDB(namePrefix);

        if(openWindow) {
	    System.out.println("Starting plots");
            JFrame frame = new JFrame("Beam Spot");
            frame.setSize(1400, 900);
            frame.add(bs.plotHistos());
            frame.setLocationRelativeTo(null);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setVisible(true);
            if(printHistos) bs.printHistos();
	        // needed to print all DC plots in a readable way, needed until I can better figure out how groot works
            // if(runDC) bs.plotDC();
        }
    }

}
