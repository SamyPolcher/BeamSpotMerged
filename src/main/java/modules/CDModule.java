package modules;

import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.List;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

import objects.Track;
import objects.Event;
import analysis.Module;

/**
 *
 * @author devita
 */
public class CDModule extends Module {
    
    private final double PMIN = 0.0;
    private final double PMAX = 2.0;
    private final double PHIMIN = -180.0;
    private final double PHIMAX = 180.0;
    private final double THETAMIN = 40.0;
    private final double THETAMAX = 100.0;
    private final double VXYMIN = -1.0;
    private final double VXYMAX =  1.0;
    private final double VZMIN = -12; //-26;
    private final double VZMAX =  15;//26;
    
    private final double CHI2PIDCUT = 10;
    
    // To store the position of the beam with respect to the default cooking value of (0, 0) and their error 
    List<Double> vx = new ArrayList<Double> ();
    List<Double> vy = new ArrayList<Double> ();
    List<Double> e_vx = new ArrayList<Double> ();
    List<Double> e_vy = new ArrayList<Double> ();
    
    public CDModule() {
        super("CDVertex_");
    }
    
    public DataGroup createVertexGroup(int col) {
        H1F hi_d0       = histo1D("hi_d0", "d0 (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H2F hi_d0phi    = histo2D("hi_d0phi", "#phi (deg)", "d0 (cm)", 30, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        GraphErrors gr  = new GraphErrors("gr_d0phi");
        gr.setTitleX("#phi (deg)");
        gr.setTitleY("d0 (cm)");
        gr.setMarkerColor(2);
        H1F hi_vz       = histo1D("hi_vz", "vz (cm)", "Counts", 100, VZMIN, VZMAX, col);
        H2F hi_vxy      = histo2D("hi_vxy", "vx (cm)", "vy (cm)", 100, VXYMIN, VXYMAX, 100, VXYMIN, VXYMAX);
        H1F hi_xb       = histo1D("hi_xb", "xb (cm)", "Counts", 10000, VXYMIN, VXYMAX, col);
        H1F hi_yb       = histo1D("hi_yb", "yb (cm)", "Counts", 10000, VXYMIN, VXYMAX, col);
        H2F hi_vxphi    = histo2D("hi_vxphi", "#phi (deg)", "vx (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vyphi    = histo2D("hi_vyphi", "#phi (deg)", "vy (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vzphi    = histo2D("hi_vzphi", "#phi (deg)", "vz (cm)", 100, PHIMIN, PHIMAX, 100, VZMIN, VZMAX);

        DataGroup dgVertex = new DataGroup(4,2);
        dgVertex.addDataSet(hi_d0,       0);
        dgVertex.addDataSet(hi_d0phi,    1);
        dgVertex.addDataSet(gr,          2);
        dgVertex.addDataSet(hi_vz,       3);
        dgVertex.addDataSet(hi_xb,       4);
        dgVertex.addDataSet(hi_yb,       4);
        dgVertex.addDataSet(hi_vxy,      4);
        dgVertex.addDataSet(hi_vxphi,    5);
        dgVertex.addDataSet(hi_vyphi,    6);
        dgVertex.addDataSet(hi_vzphi,    7);
        return dgVertex;
    }


    @Override
    public boolean checkTrack(Track trk) {
        if(trk.getNDF()<1 || trk.getChi2()/trk.getNDF()>30 || trk.pt()<0.2) return false;
        return true;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("UNegatives", this.createVertexGroup(43));
        this.getHistos().put("UPositives", this.createVertexGroup(47));
        this.getHistos().put("Negatives",  this.createVertexGroup(44));
        this.getHistos().put("Positives",  this.createVertexGroup(42));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Track> trackPos = new ArrayList<>();
        List<Track> trackNeg = new ArrayList<>();
        List<Track> utrackPos = new ArrayList<>();
        List<Track> utrackNeg = new ArrayList<>();
        for(Track track : event.getCDTracks()) {
            if(checkTrack(track)) {
                if(track.charge()>0) trackPos.add(track);
                else                 trackNeg.add(track);
            }
        }
        for(Track track : event.getCDUTracks()) {
            if(checkTrack(track)) {
                if(track.charge()>0) utrackPos.add(track);
                else                 utrackNeg.add(track);
            }
        }
        this.fillGroup(this.getHistos().get("Positives"), trackPos);
        this.fillGroup(this.getHistos().get("Negatives"), trackNeg);
        this.fillGroup(this.getHistos().get("UPositives"), utrackPos);
        this.fillGroup(this.getHistos().get("UNegatives"), utrackNeg);
    }
    
    public void fillGroup(DataGroup group, List<Track> tracks) {
        for(Track track : tracks) {
            group.getH1F("hi_d0").fill(track.d0());
            group.getH2F("hi_d0phi").fill(Math.toDegrees(track.phi()),track.d0());
            group.getH1F("hi_vz").fill(track.vz());
            group.getH2F("hi_vxy").fill(track.vx(),track.vy());
            group.getH2F("hi_vxphi").fill(Math.toDegrees(track.phi()),track.vx());
            group.getH2F("hi_vyphi").fill(Math.toDegrees(track.phi()),track.vy());
            group.getH2F("hi_vzphi").fill(Math.toDegrees(track.phi()),track.vz());
        }
    }
   
    @Override
    public void analyzeHistos() {
        this.analyzeGroup("Positives", false);
        this.analyzeGroup("Negatives", false);
        this.analyzeGroup("UPositives", true);
        this.analyzeGroup("UNegatives", true);
    }
    
    private void analyzeGroup(String name, boolean findVertex) {
        H2F h2         = this.getHistos().get(name).getH2F("hi_d0phi");
        GraphErrors gr = this.getHistos().get(name).getGraph("gr_d0phi");
        F1D f1 = this.fitD0Phi(h2, gr);
        
        if(findVertex) {
            System.out.printf("\nAnalyzing Vertex group: " + name +"\n");
            System.out.printf("d0(phi) = p0 sin(p1 x + p2):\n");
            for(int i=0; i<f1.getNPars(); i++)
                System.out.printf("\t p%d = (%.4f +/- %.4f)\n", i, f1.getParameter(i), f1.parameter(i).error());
           
            vx.add( -10*f1.getParameter(0)*Math.cos(f1.getParameter(2)) );
            vy.add( 10*f1.getParameter(0)*Math.sin(f1.getParameter(2)) );
            e_vx.add( 10*Math.sqrt(Math.pow(f1.parameter(0).error()*Math.cos(f1.getParameter(2)),2)+
                                      Math.pow(f1.getParameter(0)*Math.sin(f1.getParameter(2))*f1.parameter(2).error(),2)) );
            e_vy.add( 10*Math.sqrt(Math.pow(f1.parameter(0).error()*Math.sin(f1.getParameter(2)),2)+
                                      Math.pow(f1.getParameter(0)*Math.cos(f1.getParameter(2))*f1.parameter(2).error(),2)) );
            System.out.printf("x_offset: (%2.3f +/- %2.3f) mm, y_offset: (%2.3f +/- %2.3f) mm\n", vx.get(vx.size() - 1), e_vx.get(e_vx.size() - 1), 
                    vy.get(vy.size() - 1), e_vy.get(e_vy.size() - 1)); // convert to mm
        }
    }
    
    private F1D fitD0Phi(H2F h2, GraphErrors gr) {
        this.fitSlices(h2, gr);
        F1D f1 = new F1D("f1","[p0]*sin([p1]*x+[p2])", PHIMIN, PHIMAX);
        f1.setParameter(0, (gr.getMax()-gr.getMin())/2.0);
        f1.setParameter(1, Math.PI/180);
        f1.setParLimits(1, Math.PI/180*0.99, Math.PI/180*1.01);
        DataFitter.fit(f1, gr, "Q");
        return f1;
    }

    
    public void fitSlices(H2F h2, GraphErrors gr) {
        gr.reset();
        ParallelSliceFitter psf = new ParallelSliceFitter(h2);
        psf.setBackgroundOrder(0);
        this.toDevNull();
        psf.fitSlicesX();
        this.restoreStdOutErr();
        
        GraphErrors amp   = psf.getAmpSlices();
        GraphErrors mean  = psf.getMeanSlices();
        GraphErrors sigma = psf.getSigmaSlices();
        for(int i=0; i<mean.getDataSize(0); i++) {
            if(amp.getDataY(i)>10 && Math.abs(sigma.getDataY(i))<h2.getSlicesX().get(i).getRMS())
                gr.addPoint(mean.getDataX(i), mean.getDataY(i), 0, sigma.getDataY(i));
        }
        if(gr.getDataSize(0)<2) {
            System.out.println(gr.getName());
            gr.addPoint(PHIMIN, VXYMIN, 0, 0);
            gr.addPoint(PHIMAX, VXYMAX, 0, 0);
        }
    }

    
    @Override
    public void writeCCDB(String outputPrefix) {
        try {
            double Vx=0, Vy=0, E_Vx2=0, E_Vy2=0;
            int N = vx.size();
            for(int i=0; i<N; i++) {
                Vx += 1./N * vx.get(i);
                Vy += 1./N * vy.get(i);
                E_Vx2 += 1./(N*N) * (e_vx.get(i)*e_vx.get(i));
                E_Vy2 += 1./(N*N) * (e_vy.get(i)*e_vy.get(i));
            }
            System.out.println("Writing to: "+outputPrefix+"_CD_ccdb_table.txt ...");
            PrintWriter wr = new PrintWriter( outputPrefix+"_CD_ccdb_table.txt" );
            wr.printf( "# x y ex ey\n" );
            wr.printf( "0 0 0 " );
            wr.printf(  "%.2f %.2f %.2f %.2f\n", Vx, Vy, Math.sqrt(E_Vx2), Math.sqrt(E_Vy2));
            wr.close();
        } catch ( IOException e ) {}
    }

    @Override
    public void setPlottingOptions(String name) {
        this.getCanvas(name).setGridX(false);
        this.getCanvas(name).setGridY(false);
        this.setLogZ(name);
        EmbeddedPad pad = this.getCanvas(name).getCanvasPads().get(4);
        pad.getAxisX().setRange(VXYMIN, VXYMAX);
        pad.getAxisY().setRange(VXYMIN, VXYMAX);
    }
    
}
