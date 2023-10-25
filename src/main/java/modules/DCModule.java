package modules;

import java.lang.Math;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.swing.JFrame;
import java.awt.Dimension;
import java.util.Vector;

import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataBank;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.DataVector;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.Axis;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.graphics.GraphicsAxis;
import org.jlab.groot.data.TDirectory;

import analysis.Module;
import objects.Track;
// import src.main.java.analysis.String;
// import trash.IOException;
// import trash.PrintWriter;
import objects.Event;

public class DCModule  extends Module {
   
    // Attributes
    // -----------------------------------------
    int binsPerSector = 10;
    float targetZ = 25.4f;
    boolean relative = false;

    // peak validity for the fit window:
    float fitMin = 0.f;
    float fitMax = 0.f;
    
    // bins used for the beamspot analysis
    double[] theta_bins;
    double[] phi_bins;
    double[] z_bins;

    // vectors to keep the track line equation + raster position
    // Vector<Float> ox = new Vector();
    // Vector<Float> oy = new Vector();
    // Vector<Float> oz = new Vector();
    // Vector<Float> cx = new Vector();
    // Vector<Float> cy = new Vector();
    // Vector<Float> cz = new Vector();
    // Vector<Float> xraster = new Vector();
    // Vector<Float> yraster = new Vector();


    // 
    // ----------------------------------------- 
    public DCModule() {
      super("DCVertex");
    }
    

    // setters
    // -----------------------------------------

    public void setThetaBins( double[] bins ) { theta_bins = bins; }

    public void setTargetZ( float z, float zmin, float zmax ) {
      targetZ = z;
      fitMin = zmin;
      fitMax = zmax;
      if(zmin==0) fitMin = targetZ - 6.f;
      if(zmax==0) fitMax = targetZ + 6.f;
    }

    public void setBinsPerSector( int n ) { binsPerSector = n; }
    public void setRelative(boolean b) {relative = b;}    
    
    // getters
    // -----------------------------------------
    public double[] getThetaBins() { return theta_bins; }
    
    public int getBinsPerSector() { return binsPerSector;}
    
    
    // initialize histograms, graphs and fit functions
    // -----------------------------------------
    @Override
    public void createHistos() {
        
      final float zmin = (int)(targetZ - 8.);
      final float zmax = (int)(targetZ + 16);
      int NphiBin = (int)(6*binsPerSector);
      
      // containers for theta bins
      DataGroup dg_z_phi = new DataGroup(1, theta_bins.length-1);    // phi vs z at vertex
      DataGroup dg_peak = new DataGroup(1, theta_bins.length-1);     // mean position of the target window versus phi
      
      // containers for z slice fits
      DataGroup dg_z_slice = new DataGroup(NphiBin, theta_bins.length-1);
      
      // containers for general 1D histograms, z and phi distributions, xb, yb distrib
      DataGroup dg_distrib = new DataGroup(1, 7);
      
      // containers for for plotting the fits results as a function of theta
      DataGroup dg_fit_results = new DataGroup(1, 5);
      
      H1F h1_z   = new H1F( "vz",   "z vertex"  , 200, -20, 50 );         dg_distrib.addDataSet(h1_z, 0);
      H1F h1_phi = new H1F( "phi", "phi distribution", 180, -30, 330 );   dg_distrib.addDataSet(h1_phi, 1);
      H1F h1_xb  = histo1D("xb", "xb (cm)", "Counts", 1000, -1, 1, 43);  dg_distrib.addDataSet(h1_xb, 2);
      H1F h1_yb  = histo1D("yb", "yb (cm)", "Counts", 1000, -1, 1, 43);  dg_distrib.addDataSet(h1_yb, 3);
      H1F h1_vx  = histo1D("vx", "vx (cm)", "Counts", 1000, -1, 1, 43);  dg_distrib.addDataSet(h1_vx, 4);
      H1F h1_vy  = histo1D("vy", "vy (cm)", "Counts", 1000, -1, 1, 43);  dg_distrib.addDataSet(h1_vy, 5);
      H1F h1_theta  = histo1D("theta", "theta (deg)", "Counts", 1000, -100, 100, 43);  dg_distrib.addDataSet(h1_theta, 6);
        
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
      
      Axis phiAxis = new Axis(NphiBin, -30, 330);
      Axis zAxis = new Axis(100, zmin, zmax );
      
      phi_bins = phiAxis.getLimits();
      z_bins = zAxis.getLimits();
      
      for( int i = 0; i<theta_bins.length-1; i++ ){
          
          H2F h2 = new H2F("z_phi_"+i, "#theta = "+(theta_bins[i]+theta_bins[i+1])/2, z_bins, phi_bins); 
          h2.setTitleX("Z vertex (cm)");
          h2.setTitleY("#phi (degrees)");
          dg_z_phi.addDataSet(h2, i);
          
          GraphErrors g = new GraphErrors("g_"+i);
          g.setTitle("#theta = "+(theta_bins[i]+theta_bins[i+1])/2);
          g.setTitleX("#phi (degrees)");
          g.setTitleY("Z vertex (cm)");
          dg_peak.addDataSet(g, i);
          
          for(int j=0; j<NphiBin; j++) {
              H1F h1 = new H1F("slice_"+i+"_"+j, "", z_bins);
              dg_z_slice.addDataSet(h1, NphiBin*i+j);
          }   
      }

      // containers for track equation
      // DataGroup dg_tracks = new DataGroup(1, 6);
      // DataVector ox = new DataVector(); DataVector oy = new DataVector(); DataVector oz = new DataVector();
      // DataVector cx = new DataVector(); DataVector cy = new DataVector(); DataVector cz = new DataVector();
      // dg_tracks.addDataSet(ox, 0); dg_tracks.addDataSet(oy, 1); dg_tracks.addDataSet(oz, 2);
      // dg_tracks.addDataSet(cx, 3); dg_tracks.addDataSet(cy, 4); dg_tracks.addDataSet(cz, 5);

      
      this.getHistos().put("distribution",  dg_distrib);
      this.getHistos().put("z_phi", dg_z_phi);
      this.getHistos().put("peak_position", dg_peak);
      this.getHistos().put("fit_result",  dg_fit_results);
      this.getHistos().put("z_slice",  dg_z_slice);
      // this.getHistos().put("track_eq",  dg_tracks);
    }

    
    @Override
    public boolean checkTrack(Track trk) {
        // if(trk.getNDF()<1 || trk.getChi2()/trk.getNDF()>30 || trk.pt()<0.2) return false;
        // if(trk.charge()>=0 || trk.p()<1.) return false;
        if(trk.pid()!=11 || trk.p()<1.) return false;
        return true;
    }
    
    
    @Override
    public void fillHistos(Event event) {
      
      for(Track track : event.getFDTracks()) {

          if(checkTrack(track)) {
              
              // compute phi and theta
              float phi = (float) Math.toDegrees( Math.atan2( track.py(), track.px()) );
              if( phi < 0 ) phi += 360.0;  // transform the phi interval from [-180,180) to [0,360)
              if( phi > 330) phi -= 360.0; // pop the split sector back together

              float theta = (float) Math.toDegrees( Math.atan2( Math.sqrt( track.px()*track.px() + track.py()*track.py()), track.pz() ) );
              this.getHistos().get("distribution").getH1F("theta").fill(theta);

              // find theta/phi bin
              int thetaBin = Arrays.binarySearch( theta_bins, theta );
              thetaBin = -thetaBin -2;
              if( thetaBin < 0 || thetaBin >= theta_bins.length - 1 ) continue;
              
              int phiBin = Arrays.binarySearch( phi_bins, phi );
              phiBin = -phiBin -2;
              if( phiBin < 0 || phiBin >= phi_bins.length - 1 ) continue;

              Point3D vertex;
              if(relative){
                vertex = new Point3D(track.vx(), track.vy(), track.vz());
              }else{
                // find vz of closest approach to a 0,0 beam position
                Line3D t = new Line3D(new Point3D(track.vx(), track.vy(), track.vz()),
                                                   new Vector3D(track.px(), track.py(), track.pz()));
                // Line3D b = new Line3D(track.xb(), track.yb(), 0, 0, 0, 1);
                Line3D b = new Line3D(-0.2, -0.0, 0, 0, 0, 1);
                // Line3D b = new Line3D(-0.13, -0.08, 0, 0, 0, 1);
                vertex = t.distance(b).lerpPoint(0);
              }
              
              // fill histograms
              this.getHistos().get("distribution").getH1F("vz").fill(vertex.z());
              this.getHistos().get("distribution").getH1F("phi").fill(phi);
              this.getHistos().get("distribution").getH1F("xb").fill(track.xb());
              this.getHistos().get("distribution").getH1F("yb").fill(track.yb());
              this.getHistos().get("distribution").getH1F("vx").fill(vertex.x());
              this.getHistos().get("distribution").getH1F("vy").fill(vertex.y());

              this.getHistos().get("z_phi").getH2F("z_phi_"+thetaBin).fill(vertex.z(), phi);
              this.getHistos().get("z_slice").getH1F("slice_"+ thetaBin+"_"+phiBin).fill(vertex.z());

              // this.getHistos().get("track_eq").getData(0).add(track.ox());
              // this.getHistos().get("track_eq").getData(1).add(track.oy());
              // this.getHistos().get("track_eq").getData(2).add(track.oz());
              // this.getHistos().get("track_eq").getData(3).add(track.cx());
              // this.getHistos().get("track_eq").getData(4).add(track.cy());
              // this.getHistos().get("track_eq").getData(5).add(track.cz());

              // ox.add(track.ox()); oy.add(track.oy()); oz.add(track.oz());
              // cx.add(track.cx()); cy.add(track.cy()); cz.add(track.cz());
              // xraster.add(track.xb()); yraster.add(track.yb());

          }
      }
    }


    @Override
    public void addDataGroup(DataGroup dg, String key){
        if(key=="distribution"){
          this.fillFromDir1D(dg, "distribution", "vz");
          this.fillFromDir1D(dg, "distribution", "phi");
          this.fillFromDir1D(dg, "distribution", "xb");
          this.fillFromDir1D(dg, "distribution", "yb");
          this.fillFromDir1D(dg, "distribution", "vx");
          this.fillFromDir1D(dg, "distribution", "vy");
        }
        else if(key=="z_phi"){
          for(int i=0; i<theta_bins.length-1; i++){
            this.fillFromDir2D(dg, "z_phi", "z_phi_"+i);
          }
        }
        else if(key=="z_slice"){
          for(int i=0; i<theta_bins.length-1; i++){
            for(int j=0; j<phi_bins.length-1; j++){
              this.fillFromDir1D(dg, "z_slice", "slice_"+i+"_"+j);
            }
          }
        }
        else System.out.println(key+" DataGroup not added to "+this.getName()+" existing groups");
    }

    // public void readDataGroup(TDirectory dir) {

    //   System.out.println("in readdatagroup");
    //   String folder = "/" + this.getName() + "/";
    //   dir.cd(folder);
    //   System.out.println("Reading from: " + folder);

    //   this.fillFromDir1D(dir, "distribution", "vz");
    //   this.fillFromDir1D(dir, "distribution", "phi");

    //   for(int i=0; i<theta_bins.length-1; i++){
    //     this.fillFromDir2D(dir, "z_phi", "z_phi_"+i);
    //     for(int j=0; j<phi_bins.length-1; j++){
    //       this.fillFromDir1D(dir, "z_slice", "slice_"+i+"_"+j);
    //     }

    //   }
    // }




    // analysis
    // ------------------------------------
    public void analyzeHistos() {
        
      GraphErrors gZ = this.getHistos().get("fit_result").getGraph("gZ");
      GraphErrors gR = this.getHistos().get("fit_result").getGraph("gR");
      GraphErrors gP = this.getHistos().get("fit_result").getGraph("gP");
      GraphErrors gX = this.getHistos().get("fit_result").getGraph("gX");
      GraphErrors gY = this.getHistos().get("fit_result").getGraph("gY");

      // loop over theta bins
      for( int i=0; i<theta_bins.length-1; i++ ){

          GraphErrors g_peak = this.getHistos().get("peak_position").getGraph("g_"+i);
          analyze( i );
          
          Func1D f = g_peak.getFunction();
          if(f == null) System.out.println("g_peak fit is empty " + i);
          
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

    public void analyze( int thetaBin) {

      GraphErrors g_peak = this.getHistos().get("peak_position").getGraph("g_"+thetaBin);

      // loop over the phi bins of the 2D histogram phi vs z
      // and fit with a gaussian around the target window position

      // loop  over the phi bins
      for( int i=0; i<phi_bins.length - 1; i++ ){

        H1F h = this.getHistos().get("z_slice").getH1F("slice_"+thetaBin+"_"+i);
        H2F h2_z_phi = this.getHistos().get("z_phi").getH2F("z_phi_"+thetaBin);

        if( h.integral() < 10 ) continue;  // to skip empty bins

        // compute rms over fitting interval
        final double rms = getRMSInInterval( h, fitMin , fitMax );

        // skip if there are not enough entries
        if( h.integral( h.getAxis().getBin(fitMin) , h.getAxis().getBin(fitMax) ) < 30 ) continue;

        // the fit function of the target window peak, a gaussian for simplicity
        F1D func = new F1D( "func_"+thetaBin+"_"+i, "[amp]*gaus(x,[mean],[sigma]) + [c] + [d]*x", fitMin, fitMax );
        func.setParameter(0, h.getBinContent( h.getMaximumBin() ) );

        func.setParameter(1, this.getMeanInInterval(h, fitMin, fitMax) );
        func.setParLimits(1, fitMin, fitMax);

        func.setParameter(2, rms );
        func.setParLimits(2, 0.3, 1);

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
        // g_peak.addPoint( 
            // h2_z_phi.getYAxis().getBinCenter( i ), func.getParameter(1),
            // 0, func.parameter(1).error() );
        // }

        g_peak.addPoint( 
            h2_z_phi.getYAxis().getBinCenter( i ), func.getParameter(1),
            0, func.getParameter(2));
      }

      // extract the modulation of the target z position versus phi by fitting the graph, the function is defined in createHistos()
      F1D f = new F1D( "fit_"+thetaBin, "[z0] - [A] * cos( x * 3.1415 / 180.0 - [phi0] )", -30, 330 );
      f.setParameter(0,28.0);
      f.setParameter(1,2.0);
      f.setParameter(2, 0.);
      f.setLineWidth(3);
      f.setLineColor(2);
      f.setOptStat(11110);

      DataFitter.fit( f, g_peak, "Q");
      f.show();
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
    
    @Override
    public boolean plotGroup(String name) {
        if(name.equals("z_slice")) return false;
        else return true;
    }
    
    @Override
    public void drawHistos() {

      // EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed( "distributions" );
      // EmbeddedCanvasTabbed canvas = this.getCanvas();

      // distributions
      H1F hvz = this.getHistos().get("distribution").getH1F("vz");
      H1F hphi = this.getHistos().get("distribution").getH1F("phi");
      H1F hxb = this.getHistos().get("distribution").getH1F("xb");
      H1F hyb = this.getHistos().get("distribution").getH1F("yb");
      H1F hvx = this.getHistos().get("distribution").getH1F("vx");
      H1F hvy = this.getHistos().get("distribution").getH1F("vy");

      this.addCanvas( "distributions" );
      EmbeddedCanvas cdis = this.getCanvas().getCanvas( "distributions" );

      cdis.divide(3,2);
      cdis.cd(0).setAxisTitleSize(18);
      cdis.draw( hvz );
      cdis.cd(1).setAxisTitleSize(18);
      cdis.draw( hphi );
      cdis.cd(2).setAxisTitleSize(18);
      cdis.draw( hxb );
      cdis.cd(3).setAxisTitleSize(18);
      cdis.draw( hyb );
      cdis.cd(4).setAxisTitleSize(18);
      cdis.draw( hvx );
      cdis.cd(5).setAxisTitleSize(18);
      cdis.draw( hvy );
      
      // final results on all bins
      GraphErrors gZ = this.getHistos().get("fit_result").getGraph("gZ");
      GraphErrors gR = this.getHistos().get("fit_result").getGraph("gR");
      GraphErrors gP = this.getHistos().get("fit_result").getGraph("gP");
      GraphErrors gX = this.getHistos().get("fit_result").getGraph("gX");
      GraphErrors gY = this.getHistos().get("fit_result").getGraph("gY");

      this.addCanvas( "parameters" );
      EmbeddedCanvas cp = this.getCanvas().getCanvas( "parameters" );

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

      // results in each theta bin
      for( int i=0; i<theta_bins.length-1; i++ ){
          
        GraphErrors g_peak = this.getHistos().get("peak_position").getGraph("g_"+i);
        H2F h2_z_phi = this.getHistos().get("z_phi").getH2F("z_phi_"+i);
        String cname = String.format("%.1f",(theta_bins[i]+theta_bins[i+1])/2);
        this.addCanvas( cname );
        EmbeddedCanvas ci = this.getCanvas().getCanvas( cname );
        ci.divide(2,1);
        ci.cd(0).setAxisTitleSize(18);
        ci.getPad(0).getAxisZ().setLog(true);
        ci.draw( h2_z_phi );
        ci.cd(1).setAxisTitleSize(18);
        ci.draw( g_peak );
      }

      // zvertex distribution fits in each bin
      for( int i=0; i<theta_bins.length-1; i++ ){
        
        String cname = String.format("fits_%.1f",(theta_bins[i]+theta_bins[i+1])/2);
        this.addCanvas( cname );
        EmbeddedCanvas ci = this.getCanvas().getCanvas( cname );
        ci.divide(4,5);
        
        for( int j=0; j<phi_bins.length-1; j++ ){
          H1F h = this.getHistos().get("z_slice").getH1F("slice_"+i+"_"+j);
          ci.cd(j).setAxisTitleSize(18);
          Func1D func = h.getFunction();
          if(func != null) {
              func.setLineColor( 2 );
              func.setLineWidth( 2 );
              func.setOptStat(1110);
          } // else System.out.println("fit for z slice " + i + ":" + j + " is empty");
          ci.setAxisLabelSize(8);
          ci.setAxisLabelSize(8);
          ci.setAxisTitleSize(8);
          ci.setAxisTitleSize(8);
          ci.getPad(j).getAxisY().setLog(true);

          ci.draw( h );
          
          if(func != null) {
              F1D fb = new F1D( "fb"+h.getName(), "[c]+[d]*x", func.getMin(), func.getMax() );
              fb.setParameter(0, func.getParameter(3) );
              fb.setParameter(1, func.getParameter(4) );
              fb.setLineColor(5);
              fb.setLineWidth(2);
              ci.draw(fb,"same");
          }
        }
      }
      this.getCanvas().setActiveCanvas( "parameters" );
      // this.setCanvas(canvas);
    }
    
    @Override
    public void writeCCDB(String outputPrefix) {
        try {
            GraphErrors gX = this.getHistos().get("fit_result").getGraph("gX");
            GraphErrors gY = this.getHistos().get("fit_result").getGraph("gY");

            //relative offset
            Double vx = gX.getFunction().getParameter(0);
            Double vy = gY.getFunction().getParameter(0);
            Double e_vx = gX.getFunction().parameter(0).error();
            Double e_vy = gY.getFunction().parameter(0).error();

            Double xb = this.getHistos().get("distribution").getH1F("xb").getMean();
            Double yb = this.getHistos().get("distribution").getH1F("yb").getMean();
            
            //absolute position
            System.out.println("Writing to: "+outputPrefix+"_DC_ccdb_table.txt ...");
            PrintWriter wr = new PrintWriter( outputPrefix+"_DC_ccdb_table.txt" );
            wr.printf( "# x y ex ey in cm\n" );
            wr.printf( "0 0 0 " );
            wr.printf(  "%.3f %.3f %.3f %.3f\n", xb+vx, yb+vy, e_vx, e_vy);
            wr.printf( "# absolute position with respect to an average beam position: %3f %3f \n", xb, yb);
            wr.close();
            System.out.println(wr);

            System.out.printf("\nDrift Chambers\n");
            System.out.printf("Absolute beamspot x: (%2.3f +/- %2.3f) cm, y: (%2.3f +/- %2.3f) cm\n", xb+vx, e_vx, yb+vy, e_vy);
            System.out.printf("  with respect to average beam position: (%2.3f, %2.3f) cm\n", xb, yb);
        } catch ( IOException e ) {}
    }
}

