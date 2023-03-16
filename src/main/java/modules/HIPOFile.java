package modules;

import java.util.ArrayList;
import java.util.List;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.GraphicsAxis;

import analysis.DCbeamSpot;

public class HIPOFile {
	
	DCbeamSpot bs;
	
	public HIPOFile(DCbeamSpot bs) {
		this.bs = bs;
	}

	public void saveHistograms() {
		
		String outputPrefix = bs.getOutputPrefix();
		ArrayList<H2F> a_h2_z_phi =  bs.getA_h2_z_phi();
		
	  	System.out.println("Writing to: "+outputPrefix+"_histos.hipo ...");
	    TDirectory d = new TDirectory();
	    d.mkdir("/slices");
	    d.cd("/slices");
	    for (H2F h : a_h2_z_phi) d.addDataSet(h);
	    d.writeFile(outputPrefix+"_histos.hipo");
	  }

	public static H2F userRebin(H2F h, int nbins) {
	    if (h.getYAxis().getNBins() < nbins) {
	      System.err.println("User Binning Ignored:  Not enough bins to rebin.");
	    }
	    else {
	      if (h.getYAxis().getNBins() % nbins != 0) {
	        final String msgfmt = "User Binning Ignored:  # of existing bins (%d) is not a multiple of requested bins (%d).";
	        System.err.println(String.format(msgfmt,h.getYAxis().getNBins(),nbins));
	      }
	      else {
	        h.rebinX((int)((float)h.getYAxis().getNBins() / nbins));
	      }
	    }
	    return h;
	  }

	 
	public void readHistograms(String filename) {
		
		double[] theta_bins = bs.getThetaBins();
		ArrayList<H2F> a_h2_z_phi =  bs.getA_h2_z_phi();
		int binsPerSector = bs.getBinsPerSector();
		
	    System.out.println("Reading from: "+filename+" ...");
	    TDirectory d = new TDirectory();
	    d.readFile(filename);
	    d.cd();
	    for (int i = 0; i<theta_bins.length-1; i++) {
	      H2F h = userRebin((H2F)d.getObject("/slices/h2_z_phi_"+i), binsPerSector*6);
	      // histograms saved in a HIPO file don't retain full
	      // attributes, so here we refill:
	      for (int ix = 0; ix<h.getXAxis().getNBins(); ix++) {
	        for (int iy = 0; iy<h.getYAxis().getNBins(); iy++) {
	          for (int iz = 0; iz<h.getBinContent(ix,iy); iz++) {
	            double x = h.getXAxis().getBinCenter(ix);
	            double y = h.getYAxis().getBinCenter(iy);
	            a_h2_z_phi.get(i).fill(x,y);
	          }
	        }
	      }
	    }
	    bs.setA_h2_z_phi(a_h2_z_phi);
	  }

	  public void readHistograms(List<String> filenames) {
	    for (String f : filenames) readHistograms(f);
	  }

}
