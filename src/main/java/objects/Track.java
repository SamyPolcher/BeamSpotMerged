package objects;

import analysis.Constants;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Track extends Particle {

    private double beta = -1;
    private double chi2pid = Double.POSITIVE_INFINITY;
    private double chi2;
    private float xb = 0;
    private float yb = 0;
    private int recStatus=0;
    private int status=0;
    private int index = -1;
    private int sector=0;
    private int NDF;
    private double solenoid = -1;

    // For DC tracks, line equation of the line extrapolation of the track to the vertex
    // ox, oy, oz, is the point in REC::Traj of smallest z
    // cx, cy, cy cosine of the track at that point
    private float ox = 0;
    private float oy = 0;
    private float oz = 0;
    private float cx = 0;
    private float cy = 0;
    private float cz = 0;
    
    
    Track(int pid, double px, double py, double pz, double vx, double vy, double vz) {
        super(pid, px, py, pz, vx, vy, vz);
    }

    Track(int pid, double px, double py, double pz, double vx, double vy, double vz, float ixb, float iyb) {
        super(pid, px, py, pz, vx, vy, vz);
        this.xb = ixb;
        this.yb = iyb;
    }
    
    public static Track readParticle(DataBank recPart, DataBank recTrack, int row) {
        int pid    = recPart.getInt("pid", row);
        int charge = recPart.getByte("charge", row);
        if(pid==0) {
            pid = charge==0 ? 22 : charge*211;
        }
        Track t = new Track(pid,
                    recPart.getFloat("px", row),
                    recPart.getFloat("py", row),
                    recPart.getFloat("pz", row),
                    recPart.getFloat("vx", row),
                    recPart.getFloat("vy", row),
                    recPart.getFloat("vz", row));
        t.setBeta(recPart.getFloat("beta", row));
        t.setChi2pid(recPart.getFloat("chi2pid", row));
        t.setRECStatus(recPart.getShort("status", row));
        if(recTrack!=null) {
            for(int j=0; j<recTrack.rows(); j++) {
                if(recTrack.getShort("pindex", j)==row) {
                    t.setIndex(recTrack.getShort("index", j));
                    t.setSector(recTrack.getByte("sector", j));
                    t.setNDF(recTrack.getShort("NDF", j));
                    t.setChi2(recTrack.getFloat("chi2", j));
                    t.setStatus(recTrack.getShort("status", j));
                    break;
                }
            }
        }       
        return t;
    }
    
    public static Track readParticle(DataBank recPart, DataBank recTrack, DataBank recUTrack, int row) {
        Track t = readParticle(recPart, recTrack, row);
        if(recUTrack!=null) {
            for(int j=0; j<recUTrack.rows(); j++) {
                if(recUTrack.getShort("index", j)==t.getIndex()) {
                    t.setVector(t.pid(),
                                recUTrack.getFloat("px", j),
                                recUTrack.getFloat("py", j),
                                recUTrack.getFloat("pz", j),
                                recUTrack.getFloat("vx", j),
                                recUTrack.getFloat("vy", j),
                                recUTrack.getFloat("vz", j));
                                break;
                }
            }
        }       
        return t;
    }

     public static Track readTrack(DataBank bank, int row) {

        int pid = bank.getInt("pid", row);
        int charge = bank.getByte("q", row);
        double pt =  bank.getFloat("pt", row);
        double tandip = bank.getFloat("tandip", row);
        double phi0 = bank.getFloat("phi0", row);
        double d0 = bank.getFloat("d0", row);
        float x0 = bank.getFloat("xb", row);
        float y0 = bank.getFloat("yb", row);

        Track t = new Track(211*charge,
            pt*Math.cos(phi0),
            pt*Math.sin(phi0),
            pt*tandip,
            -d0 * Math.sin(phi0) + x0,
            d0 * Math.cos(phi0) + y0,
            bank.getFloat("z0", row),
            x0,
            y0
        );
        t.setNDF(bank.getShort("ndf", row));
        t.setChi2(bank.getFloat("chi2", row));
        t.setStatus(bank.getShort("status", row));
        return t;
    }
    
    public void addScale(DataBank config) {
        this.solenoid = config.getFloat("solenoid", 0);
    }

    public void setxbyb(float xb, float yb){
        this.xb = xb;
        this.yb = yb;
    }

    public float xb(){
        return this.xb;
    }

    public float yb(){
        return this.yb;
    }

    public void addTraj(DataBank traj, int pindex){
        float minz = 9999;
        for (int i = 0; i < traj.rows(); i++) {
            
            int rowpindex = traj.getInt("pindex", i);
            if(rowpindex != pindex){ continue;}
            float rowz = traj.getFloat("z", i);

            if( rowz < minz ){
                minz = rowz;
                ox = traj.getFloat("x", i);
                oy = traj.getFloat("y", i);
                oz = rowz;

                cx = traj.getFloat("cx", i);
                cy = traj.getFloat("cy", i);
                cz = traj.getFloat("cz", i);
            }
        }
    }

    public float ox(){
        return this.ox;
    }

    public float oy(){
        return this.oy;
    }

    public float oz(){
        return this.oz;
    }

    public float cx(){
        return this.cx;
    }

    public float cy(){
        return this.cy;
    }

    public float cz(){
        return this.cz;
    }

    public int getNDF() {
        return NDF;
    }

    public void setNDF(int NDF) {
        this.NDF = NDF;
    }

    public int getIndex() {
        return index;
    }


    public void setIndex(int index) {
        this.index = index;
    }

    public double getChi2() {
        return chi2;
    }

    public void setChi2(double chi2) {
        this.chi2 = chi2;
    }

    public int getStatus() {
        return status;
    }

    public void setStatus(int status) {
        this.status = status;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    public double getChi2pid() {
        return chi2pid;
    }

    public void setChi2pid(double chi2pid) {
        this.chi2pid = chi2pid;
    }

    public int getDetector() {
        return (int) Math.abs(this.recStatus)/1000;
    }

    public int getRECStatus() {
        return recStatus;
    }

    public void setRECStatus(int status) {
        this.recStatus = status;
    }

    public int getSector() {
        return sector;
    }

    public void setSector(int sector) {
        this.sector = sector;
    }

    public double pt() {
        return this.p()*Math.sin(this.theta());
    }
    
    public double rho() {
        return PhysicsConstants.speedOfLight()*Constants.B/this.pt()/1E5*10;
    }
   
    public double tandip() {
        return 1/Math.tan(this.theta());
    }
    
    public double d0() {
        double kappa = -this.charge()/this.pt();
        double xcen = this.vx() + Math.signum(kappa) * Constants.ALPHA * this.py();
        double ycen = this.vy() - Math.signum(kappa) * Constants.ALPHA * this.px();
        double phi0 = Math.atan2(ycen, xcen);
        if (Math.signum(kappa) < 0) {
            phi0 = Math.atan2(-ycen, -xcen);
        }
        double drh0 = (xcen-xb)*Math.cos(phi0) + (ycen-yb)*Math.sin(phi0) - Constants.ALPHA/ kappa;
        return -drh0;
//        return Math.signum(this.vy()/Math.cos(this.phi()))*Math.sqrt(this.vx()*this.vx()+this.vy()*this.vy());
    }


    public double tx() {
        return this.px()/this.py();
    }
    
    public double tz() {
        return this.pz()/this.py();
    }
    

    public double deltaPhi(Track o) {
        double dphi = this.phi()-o.phi();
        if(Math.abs(this.phi()-o.phi())>2*Math.PI) dphi -= Math.signum(this.phi()-o.phi())*2*Math.PI;
        return dphi;
    }
    
    public boolean match(Particle p) {
        double dp = (this.p()-p.p())/p.p();
        double dth = Math.toDegrees(this.theta()-p.theta());
        double dph = Math.toDegrees(this.phi()-p.phi());
        if(Math.abs(dph)>Math.PI) dph -= Math.signum(dph)*2*Math.PI;
        
        if(this.solenoid!=0 && this.charge()!=p.charge()) return false;
        else if(this.solenoid!=0 && Math.abs(dp)>Constants.NSIGMA*Constants.SIGMA_P) return false;
        else if(Math.abs(dth)>Constants.NSIGMA*Constants.SIGMA_THETA) return false;
        else if(Math.abs(dph)>Constants.NSIGMA*Constants.SIGMA_PHI) return false;
        else return true;
    }
    
    public boolean hasCharge(int charge) {
        if(charge!=0) 
            return this.charge()==charge;
        else
            return true;
    }
}
