package modules;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.swing.JFrame;
import java.awt.Dimension;

import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataBank;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.graphics.GraphicsAxis;

import analysis.Module;
import objects.Track;
import objects.Event;

public class DCModule  extends Module {
    
    String outputPrefix = "DCbeamSpot";

    float fitRangeScale = 1.0f;
    float targetZ = 25.4f;
    int binsPerSector = 10;
    int NphiBin = binsPerSector * 6;
    boolean isInit;

    // settings
    // -----------------------------------------
    boolean check_slices;
    double[] theta_bins;

    // 
    // ----------------------------------------- 
    public DCModule() {
      super("DCVertex");
      this.outputPrefix = "DCVertex";
      check_slices = true;
      isInit = false;
    }

    // initialize histograms, graphs and fit functions
    // -----------------------------------------
    @Override
    public void createHistos() {
        
      final float zmin = (int)(targetZ - 4.4);
      final float zmax = (int)(targetZ + 15.6);
      NphiBin = (int)(6*binsPerSector);
      
      // containers for theta bins
      DataGroup dg_z_phi = new DataGroup(1, theta_bins.length-1);    // phi vs z at vertex
      DataGroup dg_peak = new DataGroup(1, theta_bins.length-1);     // mean position of the target window versus phi
      DataGroup dg_fits = new DataGroup(1, theta_bins.length-1);     // fit functions of the target window position modulation in phi
      
      // containers for z slice fits
      DataGroup dg_z_slice = new DataGroup(NphiBin, theta_bins.length-1);
      
      // containers for general 1D histograms, z and phi distributions
      DataGroup dg_distrib = new DataGroup(1, 2);
      
      // containers for for plotting the fits results as a function of theta
      DataGroup dg_fit_results = new DataGroup(1, 5);
      
      H1F h1_z   = new H1F( "vz",   "z vertex"  , 200, -20, 50 ); dg_distrib.addDataSet(h1_z, 0);
      H1F h1_phi = new H1F( "phi", "phi distribution", 180, -30, 330 ); dg_distrib.addDataSet(h1_phi, 1);
        
      // graphs for plotting the results as a function of theta
      GraphErrors gZ = new GraphErrors("gZ");  // Z 
      GraphErrors gR = new GraphErrors("gR");  // R
      GraphErrors gP = new GraphErrors("gP");  // phi
      GraphErrors gX = new GraphErrors("gX");  // x 
      GraphErrors gY = new GraphErrors("gY");  // y

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
      
      dg_fit_results.addDataSet(gZ, 0);
      dg_fit_results.addDataSet(gR, 1);
      dg_fit_results.addDataSet(gP, 2);
      dg_fit_results.addDataSet(gX, 3);
      dg_fit_results.addDataSet(gY, 4);
      
      
      for( int i = 0; i<theta_bins.length-1; i++ ){
          System.out.println("createHistos " + i);
          H2F h2 = histo2D("z_phi_"+i, "Z vertex (cm)", "#phi (degrees)", 100, zmin, zmax, NphiBin, -30, 330); 
          h2.setTitle("#theta = "+(theta_bins[i]+theta_bins[i+1])/2);
          dg_z_phi.addDataSet(h2, i);
          
          GraphErrors g = new GraphErrors("g_"+i);
          g.setTitle("#theta = "+(theta_bins[i]+theta_bins[i+1])/2);
          g.setTitleX("#phi (degrees)");
          g.setTitleY("Z vertex (cm)");
          dg_peak.addDataSet(g, i);
          
          F1D f = new F1D( "fit_"+i, "[z0] - [A] * cos( x * 3.1415 / 180.0 - [phi0] )", -30, 330 );
          f.setParameter(0,28.0);
          f.setParameter(1,2.0);
          f.setParameter(2, 0.);
          f.setLineWidth(3);
          f.setLineColor(2);
          f.setOptStat(11110);
          dg_fits.addDataSet(f, i);
          
          for(int j=0; j<NphiBin; j++) {
              H1F h1 = new H1F("slice_"+i+"_"+j, "", NphiBin, -30, 330);
              dg_z_slice.addDataSet(h1, NphiBin*i+j);
          }   
      }
      
      this.getHistos().put("distribution",  dg_distrib);
      this.getHistos().put("z_phi", dg_z_phi);
      this.getHistos().put("peak_position", dg_peak);
      this.getHistos().put("fit",  dg_fits);
      this.getHistos().put("fit_result",  dg_fit_results);
      this.getHistos().put("z_slice",  dg_z_slice);
}

    // setters
    // -----------------------------------------
    public void setCheckSlices( boolean t ) { check_slices = t; }

    public void setThetaBins( double[] bins ) { theta_bins = bins; }

    public void setFitRangeScale( float s ) { fitRangeScale = s; }

    public void setTargetZ( float z ) { targetZ = z; }

    public void setBinsPerSector( int n ) { binsPerSector = n; }
        
    
    // getters
    // -----------------------------------------
    public double[] getThetaBins() { return theta_bins; }
    
    public int getBinsPerSector() { return binsPerSector;}

    public String getOutputPrefix() { return outputPrefix; }
    
    /*
    public ArrayList<H2F> getA_h2_z_phi() { return a_h2_z_phi; }
    
    public GraphErrors getGZ() { return gZ;}
    public GraphErrors getGR() { return gR;}
    public GraphErrors getGP() { return gP;}
    public GraphErrors getGX() { return gX;}
    public GraphErrors getGY() { return gY;}
    */

    
    @Override
    public boolean checkTrack(Track trk) {
        if(trk.getDetector()!=2 || trk.charge()>=0) return false;
        // if(trk.getNDF()<1 || trk.getChi2()/trk.getNDF()>30 || trk.pt()<0.2) return false;
        if(trk.getId()!=11 || trk.p()<1.5) return false;
        return true;
    }
    
    @Override
    public void fillHistos(Event event) {

//      DataBank bpart = event.getBank( "REC::Particle" );
//      DataBank btrk  = event.getBank( "REC::Track" );
//
//      if( bpart == null || btrk == null) return;
      
      for(Track track : event.getTracks()) {
          if(checkTrack(track)) {
              
              // compute phi and theta
              float phi = (float) Math.toDegrees( Math.atan2( track.py(), track.px()) );
              if( phi < 0 ) phi += 360.0;  // transform the phi interval from [-180,180) to [0,360)
              if( phi > 330) phi -= 360.0; // pop the split sector back together

              float theta = (float) Math.toDegrees( Math.atan2( Math.sqrt( track.px()*track.px() + track.py()*track.py()), track.pz() ) );
              
              // find theta/phi bin
              int bin = Arrays.binarySearch( theta_bins, theta );
              bin = -bin -2;
              if( bin < 0 || bin >= theta_bins.length - 1 ) continue;
              
              // fill histograms
              this.getHistos().get("distribution").getH1F("vz").fill(track.vz());
              this.getHistos().get("distribution").getH1F("phi").fill(phi);
              this.getHistos().get("z_phi").getH2F("z_phi_"+bin).fill(track.vz(), phi);
          }
      }
    }


    // analysis
    // ------------------------------------
    public void analyze() {
        
      GraphErrors gZ = this.getHistos().get("fit_result").getGraph("gZ");
      GraphErrors gR = this.getHistos().get("fit_result").getGraph("gR");
      GraphErrors gP = this.getHistos().get("fit_result").getGraph("gP");
      GraphErrors gX = this.getHistos().get("fit_result").getGraph("gX");
      GraphErrors gY = this.getHistos().get("fit_result").getGraph("gY");

      // loop over theta bins
      for( int i=0; i<theta_bins.length-1; i++ ){
        
          F1D f = this.getHistos().get("fit").getF1D("fit_"+i);
          analyze( i, f );
          
          double theta = (theta_bins[i] + theta_bins[i+1])/2.;
          double Etheta = 0.0;
          
          double Z = f.getParameter( 0 );
          double EZ = f.parameter( 0 ).error();

          double R = f.getParameter(1) * Math.tan( Math.toRadians( theta ) );
          double ER = f.parameter(1).error() * Math.tan( Math.toRadians( theta ) );

          double P = Math.IEEEremainder( Math.toDegrees(f.getParameter(2))+180, 360) + 180;
          double EP = Math.toDegrees( f.parameter(2).error() );

          double X = R * Math.cos( f.getParameter(2) ); 
          double Y = R * Math.sin( f.getParameter(2) );

          double EX = Math.sqrt( Math.pow(Math.cos(f.getParameter(2))*ER,2) + Math.pow(R*Math.sin(f.getParameter(2))*f.parameter(2).error(),2) );

          double EY = Math.sqrt( Math.pow(Math.sin(f.getParameter(2))*ER,2) + Math.pow(R*Math.cos(f.getParameter(2))*f.parameter(2).error(),2) );

          // munge the signs for more human-friendly plots:
          if (R < 0)  P = Math.IEEEremainder( P + 180, 360 );
          R = Math.abs(R);
          
          gZ.addPoint(theta, Z, Etheta, EZ);
          gR.addPoint(theta, R, Etheta, ER);
          gP.addPoint(theta, P, Etheta, EP);
          gX.addPoint(theta, X, Etheta, EX);
          gY.addPoint(theta, Y, Etheta, EY);
          
      }

      fitPol0( gZ );
      fitPol0( gR );
      fitPol0( gP );
      fitPol0( gX );
      fitPol0( gY );
    }


    // analysis of one theta bin
    // ------------------------------------

    public void analyze( int i_theta_bin, F1D fitFunc ) {

      GraphErrors g_peak = this.getHistos().get("peak_position").getGraph("g_"+i_theta_bin);
      H2F h2_z_phi = this.getHistos().get("z_phi").getH2F("z_phi_"+i_theta_bin);

      // loop over the phi bins of the 2D histogram phi vs z
      // and fit with a gaussian around the target window position

      // peak validity window:
      final double xmin = targetZ - 6.;
      final double xmax = targetZ + 6.;

      // loop  over the phi bins
      for( int i=0; i<h2_z_phi.getYAxis().getNBins(); i++ ){

        // fill the z slice for a given phi and theta bin
        H1F h = this.getHistos().get("z_slice").getH1F("slice_"+i_theta_bin+"_"+i);
        for (int y = 0; y < h2_z_phi.getXAxis().getNBins(); y++) {
            h.setBinContent(y, h2_z_phi.getBinContent(y, i));
        }

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

        // store the fit result in the corresponding graph
        g_peak.addPoint( 
            h2_z_phi.getYAxis().getBinCenter( i ),
            func.getParameter(1),
            0,
            func.parameter(1).error() );

      }

      // extract the modulation of the target z position versus phi by fitting the graph, the function is defined in createHistos()
      DataFitter.fit( fitFunc, g_peak, "Q");
      fitFunc.show();
      
      if(g_peak.getFunction() == null) System.out.println("g_peak its empty " + i);
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

      // check tthat the min and max are inside the axis range
      if( max > h.getAxis().max() ) max = h.getAxis().max() - 0.00001;
      if( min < h.getAxis().min() ) min = h.getAxis().min() + 0.00001;

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

      // check tthat the min and max are inside the axis range
      if( max > h.getAxis().max() ) max = h.getAxis().max() - 0.00001;
      if( min < h.getAxis().min() ) min = h.getAxis().min() + 0.00001;
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
        
        for( int j=0; j<NphiBin; j++ ){
          H1F h = this.getHistos().get("z_slice").getH1F("slice_"+i+"_"+j);
          ci.cd(j).setAxisTitleSize(18);
          Func1D func = h.getFunction();
          if(func == null) System.out.println("its empty " + i);
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
        }
      }

      JFrame czframe = new JFrame("Beam Spot - Gaussian Fits");
      czframe.add(czfits);
      czframe.pack();
      czframe.setMinimumSize( new Dimension( 1400,904 ) );
      czframe.setVisible(true);

      EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed( "Parameters" );
      
      for( int i=0; i<theta_bins.length-1; i++ ){
          
        GraphErrors g_peak = this.getHistos().get("peak_position").getGraph("g_"+i);
        H2F h2_z_phi = this.getHistos().get("z_phi").getH2F("z_phi_"+i);
        String cname = String.format("%.1f",(theta_bins[i]+theta_bins[i+1])/2);
        canvas.addCanvas( cname );
        EmbeddedCanvas ci = canvas.getCanvas( cname );
        ci.divide(2,1);
        ci.cd(0).setAxisTitleSize(18);
        ci.draw( h2_z_phi );
        ci.cd(1).setAxisTitleSize(18);
        ci.draw( g_peak );
      }
      
      GraphErrors gZ = this.getHistos().get("fit_result").getGraph("gZ");
      GraphErrors gR = this.getHistos().get("fit_result").getGraph("gR");
      GraphErrors gP = this.getHistos().get("fit_result").getGraph("gP");
      GraphErrors gX = this.getHistos().get("fit_result").getGraph("gX");
      GraphErrors gY = this.getHistos().get("fit_result").getGraph("gY");
      
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

      // save plots as png files
      if (write){
        for( int i=0; i<theta_bins.length-1; i++ ){
          String cname = String.format("%.1f",(theta_bins[i]+theta_bins[i+1])/2);
          EmbeddedCanvas ci = canvas.getCanvas( cname );
          ci.save( outputPrefix+"_bin"+i+".png");
        }
        cp.save(outputPrefix+"_results.png");
      }
    }
}

