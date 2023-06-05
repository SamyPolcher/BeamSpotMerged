package analysis;

import org.jlab.clas.pdg.PhysicsConstants;

/**
 *
 * @author devita
 */
public class Constants {
    
    public static final int BMTSECTORS = 3;
    public static final int BMTREGIONS = 3;
    public static final int SVTLAYERS = 6;
    public static final int[] SVTSECTORS = {10, 10, 14, 14, 18, 18};
    public static final int SVTSTRIPS = 200;
    public static final int[] BMTCLAYERS = {1,4,6};
    public static final int[] BMTZLAYERS = {2,3,5};
    public static final int[] BMTCSTRIPS = {896, 1024, 1152};
    public static final int[] BMTZSTRIPS = {640, 640, 768};
    public static final int[] BMTZPITCH  = {487, 536, 529};
    public static final int[] BMTMEANPHI = {-150, 90, -30};
    public static final double[] BMTRADIUS = {14.7746, 16.2646, 17.7646, 19.2646, 20.7646, 22.2646}; // radius at half drift
    
    public static final double B = 5; //field in Tesla
    public static final double ALPHA = 1E6/PhysicsConstants.speedOfLight()/B;

    public static boolean MODE = false;
    public static int PID = 2212;
    public static int CHARGE = 1;
    public static double THMIN = 30;
    
    public static final int NSIGMA = 5;
    public static final double SIGMA_P = 5; // 3;
    public static final double SIGMA_PHI = 0.3;//0.57; //0.3; //0.12;
    public static final double SIGMA_THETA = 0.5;//0.57; //0.3; //0.2;

    public static double TARGETPOS = -3; // cm
    public static double[] BEAMSPOT = new double[2];

    public static boolean getFASTMODE() {
        return MODE;
    }

    public static void setFASTMODE(boolean MODE) {
        Constants.MODE = MODE;
    }
    
    public static void setPID(int PID) {
        Constants.PID = PID;
    }

    public static void setCharge(int PID) {
        Constants.CHARGE = (int) Math.signum(PID);
    }

    public static void setBEAMSPOT(double[] BEAMSPOT) {
        Constants.BEAMSPOT = BEAMSPOT;
    }
    
    public static double getXBeam() {
        return BEAMSPOT[0];
    }
    
    public static double getYBeam() {
        return BEAMSPOT[1];
    }
}
