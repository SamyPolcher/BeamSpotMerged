package objects;
import org.jlab.io.base.DataBank;

public class Particle {
	
    public int pid;    
    public float px;     
    public float py;     
    public float pz;
    
    public float P;
    public float phi;
    public float theta;
    
    public float vx;
	public float vy;
    public float vz;     
    public float vt;
    
    public byte charge; 
    public float beta;   
    public float chi2pid;
    public short status; 

    public Particle() { }
    
    public void load( DataBank bank, int i){
      pid       = bank.getInt( "pid",i);
      px        = bank.getFloat( "px",i);
      py        = bank.getFloat( "py",i);
      pz        = bank.getFloat( "pz",i);
      vx        = bank.getFloat( "vx",i);
      vy        = bank.getFloat( "vy",i);
      vz        = bank.getFloat( "vz",i);
      vt        = bank.getFloat( "vt",i);
      charge    = bank.getByte( "charge",i);
      beta      = bank.getFloat( "beta",i);
      chi2pid   = bank.getFloat( "chi2pid",i);
      status    = bank.getShort( "status",i);
      
      // compute phi and theta and P
      P = (float) Math.sqrt( px * px + py * py + pz * pz );
      
      phi = (float) Math.toDegrees( Math.atan2(py, px) );
      if( phi < 0 ) phi += 360.0;  // transform the phi interval from [-180,180) to [0,360)
      if( phi > 330) phi -= 360.0; // pop the split sector back together

      theta = (float) Math.toDegrees( Math.atan2( Math.sqrt(px*px + py*py), pz ) );
    }
    
    public boolean checkElectron (){

        // the particle should be an electron
        if ( pid != 11 ) return false;

        // the particle momentum must be bigger than 1.5 GeV/c
        if ( P < 1.5 ) return false;

        // TODO additional cuts
        return true; 
    }
    
    public int getPid() {
		return pid;
	}
	public void setPid(int pid) {
		this.pid = pid;
	}
	public float getPx() {
		return px;
	}
	public void setPx(float px) {
		this.px = px;
	}
	public float getPy() {
		return py;
	}
	public void setPy(float py) {
		this.py = py;
	}
	public float getPz() {
		return pz;
	}
	public void setPz(float pz) {
		this.pz = pz;
	}
	public float getVx() {
		return vx;
	}
	public void setVx(float vx) {
		this.vx = vx;
	}
	public float getVy() {
		return vy;
	}
	public void setVy(float vy) {
		this.vy = vy;
	}
	public float getVz() {
		return vz;
	}
	public void setVz(float vz) {
		this.vz = vz;
	}
	public float getVt() {
		return vt;
	}
	public void setVt(float vt) {
		this.vt = vt;
	}
	public byte getCharge() {
		return charge;
	}
	public void setCharge(byte charge) {
		this.charge = charge;
	}
	public float getBeta() {
		return beta;
	}
	public void setBeta(float beta) {
		this.beta = beta;
	}
	public float getChi2pid() {
		return chi2pid;
	}
	public void setChi2pid(float chi2pid) {
		this.chi2pid = chi2pid;
	}
	public short getStatus() {
		return status;
	}
	public void setStatus(short status) {
		this.status = status;
	}
	public float getP() {
		return P;
	}
	public void setP(float p) {
		P = p;
	}
	public float getPhi() {
		return phi;
	}
	public void setPhi(float phi) {
		this.phi = phi;
	}
	public float getTheta() {
		return theta;
	}
	public void setTheta(float theta) {
		this.theta = theta;
	}
	
  }
