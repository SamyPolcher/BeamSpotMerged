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
    
    public BeamSpot(String opts, double[] thetaBins, float fitScale, float Ztarget, int NphiBins) {
        this.init(opts, thetaBins, fitScale, Ztarget, NphiBins);
    }
    

    private void init(String opts, double[] thetaBins, float fitScale, float Ztarget, int NphiBins) {
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
        
        CDModule cd = new CDModule();
        cd.init();
        
        DCModule dc = new DCModule();
        dc.setThetaBins(thetaBins);
        dc.setTargetZ(Ztarget);
        dc.setBinsPerSector(NphiBins);
        dc.setFitRangeScale(fitScale);
        dc.init();
        
        this.modules.add(cd);
        this.modules.add(dc);
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
    }

    private void printHistos() {
        System.out.println("\n>>>>> Printing canvases to directory plots");
        for(Module m : modules) {
            m.printHistos("plots");
        }
    }
    
    private void testHistos() {
        for(Module m : modules) {
            m.testHistos();
        }
    }
    
    public static void main(String[] args) {
        
        OptionParser parser = new OptionParser("Beam Spot [options] file1 file2 ... fileN");
        parser.setRequiresInputList(false);
        // valid options for event-base analysis
        parser.addOption("-o"          ,"",     "histogram file name prefix");
        parser.addOption("-n"          ,"-1",   "maximum number of events to process");
        parser.addOption("-x"          ,"0",   "do NOT save histograms in a hipo file");
        // histogram settings
        parser.addOption("-histo"      ,"0",    "read histogram from hipo file (0/1)");
        parser.addOption("-plot"       ,"1",    "display histograms (0/1)");
        parser.addOption("-print"      ,"0",    "print canvases (0/1)");
        parser.addOption("-stats"      ,"",     "histogram stat option (e.g. \"10\" will display entries)");
        // DC analysis settigs
        parser.addOption("-scale", "1.0", "Fit range scale factor");
        parser.addOption("-Zvertex", "25.4", "Nominal Z of Target/Foil");
        parser.addOption("-Nphi", "10", "Phi bins per sector");

        parser.parse(args);
        
        String namePrefix  = parser.getOption("-o").stringValue();        
        String histoName   = "histo.hipo";
        if(!namePrefix.isEmpty()) {
            histoName  = namePrefix + "_" + histoName; 
        }
        int     maxEvents     = parser.getOption("-n").intValue();
        boolean saveHistos    = (parser.getOption("-x").intValue()!=0);
        boolean readHistos    = (parser.getOption("-histo").intValue()!=0);            
        boolean openWindow    = (parser.getOption("-plot").intValue()!=0);
        boolean printHistos   = (parser.getOption("-print").intValue()!=0);
        String  optStats      = parser.getOption("-stats").stringValue(); 
        float  fitScale       = (float)parser.getOption("-scale").doubleValue();
        float  zTarget        = (float)parser.getOption("-Zvertex").doubleValue();
        int NphiBins          = parser.getOption("-Nphi").intValue();
        
        
        if(!openWindow) System.setProperty("java.awt.headless", "true");
        
        double[] thetaBins = new double[]{10,11,12,13,14,16,18,22,30};
        
        BeamSpot bs = new BeamSpot(optStats, thetaBins, fitScale, zTarget, NphiBins);
        
        List<String> inputList = parser.getInputList();
        if(inputList.isEmpty()==true){
            parser.printUsage();
            System.out.println("\n >>>> error: no input file is specified....\n");
            System.exit(0);
        }

        if(readHistos) {
            bs.readHistos(inputList.get(0));
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
            }    
            bs.analyzeHistos();
            if(saveHistos) bs.saveHistos(histoName);
        }

        if(openWindow) {
            JFrame frame = new JFrame("Beam Spot");
            frame.setSize(1400, 900);
//            frame.add(bs.plotHistos());
            frame.add(bs.modules[1].plot());
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);
            if(printHistos) bs.printHistos();
        }
    }

}