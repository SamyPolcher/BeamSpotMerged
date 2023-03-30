package analysis;

import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import org.jlab.utils.benchmark.BenchmarkTimer;
import org.jlab.utils.options.OptionParser;
import org.jlab.io.hipo.HipoDataSource;

import objects.Track;
import objects.Event;
import modules.DCModule;

public class BeamSpot {

    public static void main(String[] args) {
        OptionParser cli = new OptionParser("BeamSpot");
        cli.addOption("-H", "0", "Interpret inputs as HIPO histogram files (instead of DSTs) and add them together");
        cli.addOption("-T", "0", "Interpret input as a TXT histogram file (instead of HIPO DSTs)");
        cli.addOption("-O", "BeamSpot", "String prefix for output file names");
        cli.addOption("-X", "0", "Run with no output files");
        cli.addOption("-B", "0", "Batch mode, no graphics");
        cli.addOption("-R", "1.0", "Fit range scale factor");
        cli.addOption("-Z", "25.4", "Nominal Z of Target/Foil");
        cli.addOption("-N", "10", "Phi bins per sector");
        cli.parse(args);
        
        if (cli.getInputList().size()==0) {
            System.err.println(cli.getUsageString());
            System.err.println("BeamSpot:  ERROR: No input files specified.");
            System.exit(1);
        }
        
        DCModule bs = new DCModule();
        
        // set the theta bin edges
        bs.setThetaBins( new double[]{10,11,12,13,14,16,18,22,30} );
        bs.setCheckSlices(true);
        bs.setFitRangeScale((float)cli.getOption("-R").doubleValue());
        
        bs.setTargetZ((float)cli.getOption("-Z").doubleValue());
        bs.setBinsPerSector(cli.getOption("-N").intValue());
        
        // call the method to properly setup all the parameters
        bs.init();
    
        // loop over input files
        int n = 0;
        BenchmarkTimer bt = new BenchmarkTimer();
        
        for( String s : cli.getInputList() ) {
            HipoDataSource reader = new HipoDataSource();
            bt.resume();
            reader.open( s );
            while(reader.hasEvent() ) {
                Event event = new Event(reader.getNextEvent());
                bs.processEvent( event );
                n += 1;
            }
        bt.pause();
        System.out.println(String.format("### EVENT RATE:  %.4f kHz",n/bt.getSeconds()/1000));
        reader.close();
        }
        
        // run the analysis
        bs.analyzeHistos();
        bs.plot(true);
    }
}
