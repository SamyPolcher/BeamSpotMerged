package objects;
import org.jlab.io.base.DataBank;


public class Track {
	
    public short pindex; 
    public byte  detector;
    public byte  sector;
    public short status;
    public byte  q;
    public float chi2;
    public short NDF;

    public Track() { }
    
    public void load( DataBank bank, int i){
      pindex   = bank.getShort( "pindex",i);     
      detector = bank.getByte(  "detector",i);   
      sector   = bank.getByte(  "sector",i);     
      status   = bank.getShort( "status",i);     
      q        = bank.getByte(  "q",i);          
      chi2     = bank.getFloat( "chi2",i);       
      NDF      = bank.getShort( "NDF",i);        
    }
    
    public boolean checkTrack (){
        // only use FD tracks
        if( detector !=  org.jlab.detector.base.DetectorType.DC.getDetectorId() ) return false;

        // only negative tracks
        if( q > 0 ) return false;

        // TODO additional cuts
        return true; 
      }
    
    
// ------------------ getters and setters ------------------
    
	public short getPindex() {
		return pindex;
	}
	public void setPindex(short pindex) {
		this.pindex = pindex;
	}
	public byte getDetector() {
		return detector;
	}
	public void setDetector(byte detector) {
		this.detector = detector;
	}
	public byte getSector() {
		return sector;
	}
	public void setSector(byte sector) {
		this.sector = sector;
	}
	public short getStatus() {
		return status;
	}
	public void setStatus(short status) {
		this.status = status;
	}
	public byte getQ() {
		return q;
	}
	public void setQ(byte q) {
		this.q = q;
	}
	public float getChi2() {
		return chi2;
	}
	public void setChi2(float chi2) {
		this.chi2 = chi2;
	}
	public short getNDF() {
		return NDF;
	}
	public void setNDF(short nDF) {
		NDF = nDF;
	}
	
  }