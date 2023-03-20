package analysis;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import org.jlab.groot.math.*;
import org.jlab.groot.ui.TCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.utils.benchmark.BenchmarkTimer;
import org.jlab.utils.options.OptionParser;
import org.jlab.groot.base.AxisAttributes;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.graphics.GraphicsAxis;
import javax.swing.JFrame;
import java.awt.Dimension;

import objects.Particle;
import objects.Track;
import objects.FitFunc;

import modules.HIPOFile;
import modules.TXTFile;

public class DCbeamSpot {

  String outputPrefix = "DCbeamSpot";

  // histograms
  // ----------------------------------------- 

  // general 1D histograms
  H1F h1_z;   // z vertex distribution
  H1F h1_phi; // phi distribution at vertex

  // containers for theta bins
  ArrayList<H2F> a_h2_z_phi;           // phi vs z at vertex
  ArrayList<GraphErrors> a_g_results;  // mean position of the target window versus phi
  ArrayList<Func1D> a_fits;            // fit functions of the target window position modulation in phi

  // containers for z slice fits
  ArrayList< ArrayList<H1F> > a_hz;

  // graphs for plotting the results as a function of theta
  GraphErrors gZ;   // Z 
  GraphErrors gR;   // R
  GraphErrors gP;   // phi
  GraphErrors gX;   // x 
  GraphErrors gY;   // y

  float fitRangeScale = 1.0f;
  float targetZ = 25.4f;
  int binsPerSector = 10;
  boolean isInit;

  // settings
  // -----------------------------------------
  boolean check_slices;
  double[] theta_bins;

  // 
  // ----------------------------------------- 
  public DCbeamSpot(String outputPrefix) {
    this.outputPrefix = outputPrefix;

    a_h2_z_phi = new ArrayList<H2F>();
    a_g_results = new ArrayList<GraphErrors>();
    a_fits = new ArrayList<Func1D>();

    h1_z   = new H1F( "h1_z",   "z vertex"  , 200, -20, 50 );
    h1_phi = new H1F( "h1_phi", "phi vertex", 180, -30, 330 );

    a_hz = new ArrayList<ArrayList<H1F>>();

    check_slices = true;
    isInit = false;
  }

  // initialize histograms
  // -----------------------------------------
  public void init() {
    final float zmin = (int)(targetZ - 4.4);
    final float zmax = (int)(targetZ + 15.6);
    final int bins = (int)(6*binsPerSector);
    
    for( int i = 0; i<theta_bins.length-1; i++ ){
      H2F h = new H2F("h2_z_phi_"+i, "#theta = "+(theta_bins[i]+theta_bins[i+1])/2,100,zmin,zmax,bins,-30,330);
      h.setTitleX("Z vertex (cm)");
      h.setTitleY("#phi (degrees)");
      GraphErrors g = new GraphErrors();
      g.setTitle("#theta = "+(theta_bins[i]+theta_bins[i+1])/2);
      g.setTitleX("#phi (degrees)");
      g.setTitleY("Z vertex (cm)");
      a_h2_z_phi.add( h );
      a_g_results.add( g );
    }
    isInit = true;
  }  

  // setters
  // -----------------------------------------
  public void setCheckSlices( boolean t ) { check_slices = t; }

  public void setThetaBins( double[] bins ) { theta_bins = bins; }

  public void setFitRangeScale( float s ) { fitRangeScale = s; }

  public void setTargetZ( float z ) { targetZ = z; }

  public void setBinsPerSector( int n ) { binsPerSector = n; }
  
  public void setA_h2_z_phi( ArrayList<H2F> a_h) { 
	  if(isInit == true) { a_h2_z_phi = a_h; }
  }
  
  
  // getters
  // -----------------------------------------
  public double[] getThetaBins() { return theta_bins; }
  
  public int getBinsPerSector() { return binsPerSector;}

  public String getOutputPrefix() { return outputPrefix; }

  public ArrayList<H2F> getA_h2_z_phi() { return a_h2_z_phi; }
  
  public GraphErrors getGZ() { return gZ;}
  public GraphErrors getGR() { return gR;}
  public GraphErrors getGP() { return gP;}
  public GraphErrors getGX() { return gX;}
  public GraphErrors getGY() { return gY;}

  
  // processEvent   
  // ----------------------------------------- 
  public boolean processEvent( DataEvent event ){

    if( event.hasBank( "REC::Particle" ) == false ||
        event.hasBank( "REC::Track" )    == false    ) return false;

    DataBank bpart = event.getBank( "REC::Particle" );
    DataBank btrk  = event.getBank( "REC::Track" );

    // interfaces to the banks
    Track trk = new Track();
    Particle part = new Particle();

    // loop over tracks
    for( int i=0; i < btrk.rows(); i++){

      trk.load( btrk, i );

      // check the quality of the track
      if( trk.checkTrack() == false ) continue;

      // get the corresponding particle
      part.load( bpart, trk.getPindex() );

      // check the quality of the particle
      if( part.checkElectron() == false ) continue;

      // find theta bin
      int bin = Arrays.binarySearch( theta_bins, part.getTheta() );
      bin = -bin -2;
      if( bin < 0 || bin >= theta_bins.length - 1 ) continue;

      // fill the histograms
      h1_z.fill( part.getVz() );
      h1_phi.fill( part.getPhi() );
      //System.out.println( "Vz: "+part.getVz() + " Phi: "+ part.getPhi());
      a_h2_z_phi.get(bin).fill( part.getVz(), part.getPhi() );
    } // end loop over tracks
    return true;
  }


  // analysis
  // ------------------------------------
  public void analyze() {

    // loop over theta bins
    for( int i=0; i<theta_bins.length-1; i++ ){
      analyze( i );
    } // end loop

    // extract the results and organize them as a function of theta
    double[] Z = new double[ theta_bins.length -1 ];
    double[] R = new double[ theta_bins.length -1 ];
    double[] P = new double[ theta_bins.length -1 ];
    double[] X = new double[ theta_bins.length -1 ];
    double[] Y = new double[ theta_bins.length -1 ];
    double[] T = new double[ theta_bins.length -1 ];
    double[] ET = new double[ theta_bins.length -1 ];
    double[] EZ = new double[ theta_bins.length -1 ];
    double[] ER = new double[ theta_bins.length -1 ];
    double[] EP = new double[ theta_bins.length -1 ];
    double[] EX = new double[ theta_bins.length -1 ];
    double[] EY = new double[ theta_bins.length -1 ];

    for( int i=0; i<theta_bins.length-1; i++ ){
      Func1D f = a_fits.get(i);

      // average theta
      T[i] = (theta_bins[i] + theta_bins[i+1])/2.;
      ET[i] = 0.0;

      Z[i] = f.getParameter( 0 );
      EZ[i] = f.parameter( 0 ).error();

      R[i] = f.getParameter(1) * Math.tan( Math.toRadians( T[i] ) );
      ER[i] = f.parameter(1).error() * Math.tan( Math.toRadians( T[i] ) );

      P[i] = Math.IEEEremainder( Math.toDegrees(f.getParameter(2))+180, 360) + 180;
      EP[i] = Math.toDegrees( f.parameter(2).error() );

      X[i] = R[i] * Math.cos( f.getParameter(2) ); 
      Y[i] = R[i] * Math.sin( f.getParameter(2) );

      EX[i] = Math.sqrt( Math.pow(Math.cos(f.getParameter(2))*ER[i],2) +
          Math.pow(R[i]*Math.sin(f.getParameter(2))*f.parameter(2).error(),2) );

      EY[i] = Math.sqrt( Math.pow(Math.sin(f.getParameter(2))*ER[i],2) +
          Math.pow(R[i]*Math.cos(f.getParameter(2))*f.parameter(2).error(),2) );

      // munge the signs for more human-friendly plots:
      if (R[i] < 0) P[i] += 180;
      R[i] = Math.abs(R[i]);
    }

    gZ = new GraphErrors("Z",   T, Z, ET, EZ );
    gR = new GraphErrors("R",   T, R, ET, ER );
    gP = new GraphErrors("Phi", T, P, ET, EP );
    gX = new GraphErrors("X",   T, X, ET, EX );
    gY = new GraphErrors("Y",   T, Y, ET, EY );

    fitPol0( gZ );
    fitPol0( gR );
    fitPol0( gP );
    fitPol0( gX );
    fitPol0( gY );
  }

  // analysis of one theta bin
  // ------------------------------------
  public void analyze( int i_theta_bin ) {

    GraphErrors g_results = a_g_results.get(i_theta_bin);
    H2F h2_z_phi = a_h2_z_phi.get(i_theta_bin);

    // loop over the phi bins of the 2D histogram phi vs z
    // and fit with a gaussian around the target window position

    // peak validity window:
    final double xmin = 20.;
    final double xmax = 31.;

    // for debug 
    ArrayList<H1F> z_slices = new ArrayList<H1F>();
    int ic = 0;

    // loop  over the phi bins
    for( int i=0;i<h2_z_phi.getYAxis().getNBins(); i++ ){

      // get the phi slice
      H1F h = h2_z_phi.sliceY( i );
      h.setTitle("");

      if( h.integral() < 10 ) continue;  // to skip empty bins

      // check if the maximum is in the  expected range for the target window
      final double hmax = h.getAxis().getBinCenter( h.getMaximumBin() ) ;
      if( hmax < xmin || hmax > xmax ) continue;

      // check the entries around the peak
      final double rms = getRMSInInterval( h, hmax - 5. , hmax + 5. );
      double rmin = h.getAxis().getBinCenter( h.getMaximumBin() ) - 2.0*rms*fitRangeScale;
      double rmax = h.getAxis().getBinCenter( h.getMaximumBin() ) + 1.5*rms*fitRangeScale;

      // truncate fit range if out of bounds:
      if (rmin < h.getAxis().getBinCenter(1)) rmin = h.getAxis().getBinCenter(1);
      if (rmax > h.getAxis().getBinCenter(h.getAxis().getNBins()-1))
        rmax = h.getAxis().getBinCenter(h.getAxis().getNBins()-1);

      // skip if there are not enough entries
      if( h.integral( h.getAxis().getBin(rmin) , h.getAxis().getBin(rmax) ) < 50 ) continue;

      // the fit function of the target window peak, a gaussian for simplicity
      // the fit range is +- RMS around the peak
      F1D func = new F1D( "func"+i, "[amp]*gaus(x,[mean],[sigma]) + [c] + [d]*x", rmin, rmax ); 
      func.setParameter(0, h.getBinContent( h.getMaximumBin() ) );
      func.setParameter(1, h.getAxis().getBinCenter( h.getMaximumBin() )  ); 
      func.setParameter(2, rms/2. );
      func.setParameter(3, 1. );
      func.setParameter(4, .01 );
      func.setOptStat(110);
      DataFitter.fit( func, h, "Q" );

      // skip if Gaussian amplitude too small:
      if (func.getParameter(0) < 8) continue;

      // skip if Gaussian sigma too small:
      if (Math.abs(func.getParameter(2)) < 0.1) continue;

      // skip if Gaussian sigma too big:
      if (Math.abs(func.getParameter(2)) > 2) continue;

      // skip if chi-square bad:
      if (func.getChiSquare()/func.getNDF() < 0.05) continue;
      if (func.getChiSquare()/func.getNDF() > 10) continue;

      // store the fir result in the corresponding graph
      g_results.addPoint( 
          h2_z_phi.getYAxis().getBinCenter( i ),
          func.getParameter(1),
          0,
          func.parameter(1).error() );

      z_slices.add( h );

    } // end loop over bins

    // debug
    a_hz.add( z_slices );

    // extract the modulation of the target z position versus phi by fitting the graph
    // the function is defined below
    FitFunc func = new FitFunc( "f1", -30., 330. );
    func.setParameter(0,28.0);
    func.setParameter(1,2.0);
    func.setParameter(2, 0.);
    func.setLineWidth(3);
    DataFitter.fit( func, g_results,"Q");
    func.setLineColor(2);
    func.setOptStat(11110);
    func.show();

    // store the fit function
    a_fits.add( func );
  }

  // useful functions
  // ---------------- 
  private void fitPol0( GraphErrors g ){
    double y = 0.;
    double ey = 0.;
    for( int i=0; i<g.getDataSize(0); i++ ) y += g.getDataY(i);
    y /= g.getDataSize(0);

    for( int i=0; i<g.getDataSize(0); i++ ) ey += (g.getDataY(i)-y)*(g.getDataY(i)-y);
    ey /= g.getDataSize(0);
    ey = Math.sqrt( ey );

    F1D f = new F1D( "fff"+g.getName(), "[mean]", g.getDataX(0), g.getDataX( g.getDataSize(0)-1 ) );
    //System.out.println( " ++++++++++++ " + f.getName() + " " + y + " " + ey);
    f.setParameter(0,y);
    f.parameter(0).setError( 2*ey );
    DataFitter.fit( f, g, "Q" );
    f.setOptStat(10);
    f.setLineColor(2);
    f.setLineWidth(2);
    f.show();
  }


  private double getMeanInInterval( H1F h, double min, double max ){

    double s = 0.;
    double n = 0.;
    int bmin = h.getAxis().getBin( min );
    int bmax = h.getAxis().getBin( max );

    for ( int i=bmin; i <= bmax; i++ ){
      double X = h.getAxis().getBinCenter(i);
      double Y = h.getBinContent(i);
      s += X * Y;
      n += Y;
    }
    return s/n;
  }

  private double getRMSInInterval( H1F h, double min, double max ){
    double m = getMeanInInterval( h, min, max );

    double s = 0.;
    double n = 0.;
    int bmin = h.getAxis().getBin( min );
    int bmax = h.getAxis().getBin( max );

    for ( int i=bmin; i <= bmax; i++ ){
      double X = h.getAxis().getBinCenter(i);
      double Y = h.getBinContent(i);
      s += (X-m)*(X-m) * Y;
      n += Y;
    }
    return Math.sqrt( s/n );
  }

  public void zoom(GraphErrors g, GraphicsAxis a) {
	final double min = g.getMin(); 
	final double max = g.getMax();
	a.setRange(min - 0.3*(max-min), max + 0.3*(max-min));
  }


  // plots
  // ------------------------------------
  public void plot(boolean write) {

    EmbeddedCanvasTabbed czfits = new EmbeddedCanvasTabbed( false );
    for( int i=0; i<theta_bins.length-1; i++ ){
      String cname = String.format("%.1f",(theta_bins[i]+theta_bins[i+1])/2);
      czfits.addCanvas( cname );
      EmbeddedCanvas ci = czfits.getCanvas( cname );
      ci.divide(7,8);
      int j=0;
      for( H1F h : a_hz.get(i) ){
        ci.cd(j).setAxisTitleSize(18);
        Func1D func = h.getFunction();
        func.setLineColor( 2 );
        func.setLineWidth( 2 );
        func.setOptStat(1110);
        ci.setAxisLabelSize(8);
        ci.setAxisLabelSize(8);
        ci.setAxisTitleSize(8);
        ci.setAxisTitleSize(8);
        ci.draw( h );
        //F1D fg = new F1D( "fg"+h.getName(), "[amp]*gaus(x,[mean],[sigma])", func.getMin(), func.getMax() );
        //fg.setParameter(0, func.getParameter(0) );
        //fg.setParameter(1, func.getParameter(1) );
        //fg.setParameter(2, func.getParameter(2) );
        //fg.setLineColor(4);
        //ci.draw(fg,"same");
        F1D fb = new F1D( "fb"+h.getName(), "[c]+[d]*x", func.getMin(), func.getMax() );
        fb.setParameter(0, func.getParameter(3) );
        fb.setParameter(1, func.getParameter(4) );
        fb.setLineColor(5);
        fb.setLineWidth(2);
        ci.draw(fb,"same");
        // can't do this because it adds to the legend ...
        //ci.draw(func,"same");
        j++;
      }
    }

    JFrame czframe = new JFrame("Beam Spot - Gaussian Fits");
    czframe.add(czfits);
    czframe.pack();
    czframe.setMinimumSize( new Dimension( 1400,904 ) );
    czframe.setVisible(true);
    
    EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed( "Parameters" );
    for( int i=0; i<theta_bins.length-1; i++ ){
      String cname = String.format("%.1f",(theta_bins[i]+theta_bins[i+1])/2);
      canvas.addCanvas( cname );
      EmbeddedCanvas ci = canvas.getCanvas( cname );
      ci.divide(2,1);
      ci.cd(0).setAxisTitleSize(18);
      ci.draw( a_h2_z_phi.get(i) );
      ci.cd(1).setAxisTitleSize(18);
      ci.draw( a_g_results.get(i) );
    }

    // plot the results as a function of theta
    gZ.setTitleX( "#theta (degrees)" );
    gZ.setTitleY( "z_0 (cm)" );
    gR.setTitleX( "#theta (degrees)" );
    gR.setTitleY( "r_0 (cm)" );
    gP.setTitleX( "#theta (degrees)" );
    gP.setTitleY( "#phi_0 (degrees)" );
    gX.setTitleX( "#theta (degrees)" );
    gX.setTitleY( "x_0 (cm)" );
    gY.setTitleX( "#theta (degrees)" );
    gY.setTitleY( "y_0 (cm)" );

    EmbeddedCanvas cp = canvas.getCanvas( "Parameters" );
    cp.divide(2,3);
    cp.cd(0).setAxisTitleSize(18);
    cp.draw( gX );
    this.zoom(gX, cp.getPad(0).getAxisY());
    cp.cd(1).setAxisTitleSize(18);
    cp.draw( gY );
    this.zoom(gY, cp.getPad(1).getAxisY());
    cp.cd(2).setAxisTitleSize(18);
    cp.draw( gZ );
    this.zoom(gZ, cp.getPad(2).getAxisY());
    cp.cd(3).setAxisTitleSize(18);
    cp.draw( gP );
    this.zoom(gP, cp.getPad(3).getAxisY());
    cp.cd(4).setAxisTitleSize(18);
    cp.draw( gR );
    this.zoom(gR, cp.getPad(4).getAxisY());

    canvas.setActiveCanvas( "Parameters" );

    JFrame frame = new JFrame("BeamSpot - Modulation Fits");
    frame.add(canvas);
    frame.pack();
    frame.setMinimumSize( new Dimension( 800, 700 ) );
    frame.setVisible(true);
    
    System.out.println("cc");

    // save plots as png files
    if (write){
      for( int i=0; i<theta_bins.length-1; i++ ){
        String cname = String.format("%.1f",(theta_bins[i]+theta_bins[i+1])/2);
        EmbeddedCanvas ci = canvas.getCanvas( cname );
        ci.save( outputPrefix+"_bin"+i+".png");
      }
      cp.save(outputPrefix+"_results.png");
    }
    System.out.println("cc");
  }


  // ------------------------------------
  // ##### MAIN ######
  // ------------------------------------
  public static void main( String[] args ){

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

    DCbeamSpot bs = new DCbeamSpot(cli.getOption("-O").stringValue());
    HIPOFile hipo = new HIPOFile(bs);
	TXTFile txt = new TXTFile(bs);

    // set the theta bin edges
    bs.setThetaBins( new double[]{10,11,12,13,14,16,18,22,30} );

    bs.setCheckSlices(true);

    bs.setFitRangeScale((float)cli.getOption("-R").doubleValue());

    if( cli.getOption("-T").stringValue().equals("0") ) { // options Z and N do not work with txt histograms
      bs.setTargetZ((float)cli.getOption("-Z").doubleValue());

      bs.setBinsPerSector(cli.getOption("-N").intValue());
    }
    else {
      System.out.println( " WARNING: -N and -Z options are incompatible with -T ");
    }

    // call the init method to properly setup all the parameters
    bs.init();

    // fill the a_h2_z_phi with the necessary data for the analysis
    if( !cli.getOption("-H").stringValue().equals("0") ) {
      hipo.readHistograms(cli.getInputList());
    }
    else if( !cli.getOption("-T").stringValue().equals("0") ) {
      txt.readHistogramsFromTXT(cli.getInputList());
    }
    else {
      // loop over input files
      int n = 0;
      BenchmarkTimer bt = new BenchmarkTimer();
      for( String s : cli.getInputList() ) {
        HipoDataSource reader = new HipoDataSource();
        bt.resume();
        reader.open( s );
        while(reader.hasEvent() ) {
          DataEvent event = reader.getNextEvent();
          bs.processEvent( event );
          n += 1;
        }// end loop on events
        bt.pause();
        System.out.println(String.format("### EVENT RATE:  %.4f kHz",n/bt.getSeconds()/1000));
        reader.close();
      	}// end loop on input files
      
        // if the Run "with no output files" has not been selected the data is stored
        if (cli.getOption("-X").stringValue().equals("0")) {
	        txt.saveHistogramsToTXT();
	        hipo.saveHistograms();
        }
    }
      
    // run the analysis
    bs.analyze();
    txt.writeResults();
      
      
    if (cli.getOption("-B").stringValue().equals("0")) {
      bs.plot(cli.getOption("-X").stringValue().equals("0"));
    }
    
  }

}