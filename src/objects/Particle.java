package objects;
import org.jlab.io.base.DataBank;

class Particle {
    public int pid;    
    public float px;     
    public float py;     
    public float pz;     
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
    }
  }
