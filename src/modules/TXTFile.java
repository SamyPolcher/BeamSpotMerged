package modules;

import java.io.BufferedReader;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jlab.groot.data.H2F;
import org.jlab.groot.math.Axis;

import analysis.DCbeamSpot;


public class TXTFile {

	DCbeamSpot bs;
	
	public TXTFile(DCbeamSpot bs) {
		this.bs = bs;
	}
	
	// save 2D histograms to txt
	public void saveHistogramsToTXT() {
		
		String outputPrefix = bs.getOutputPrefix();
		double[] theta_bins = bs.getThetaBins();
		ArrayList<H2F> a_h2_z_phi =  bs.getA_h2_z_phi();
	    
		try {
	      System.out.println("Writing to: "+outputPrefix+"_histos.txt ...");
	      FileWriter wr  = new FileWriter(outputPrefix+"_histos.txt");
	
	      // write a header with the theta binning
	      wr.write( "# theta bin edges\n# "); 
	      for( int i=0; i<theta_bins.length-1; i++){
	        wr.write( theta_bins[i] +",");
	      }
	      wr.write( theta_bins[theta_bins.length-1] +"\n");
	      // then write the x and y axis bins
	      H2F h0 = a_h2_z_phi.get(0);
	
	      wr.write( "# " );
	      Axis x = h0.getXAxis();
	      for( int j=0;j<x.getLimits().length-1;j++){
	        wr.write( x.getLimits()[j] + "," );
	      }
	      wr.write( x.getLimits()[x.getLimits().length-1] + "\n" );
	      wr.write( "# " );
	      Axis y = h0.getYAxis();
	      for( int j=0;j<y.getLimits().length-1;j++){
	        wr.write( y.getLimits()[j] + "," );
	      }
	      wr.write( y.getLimits()[y.getLimits().length-1] + "\n" );
	
	      // for each theta bin, write a header for the theta bin
	      for( int i=0; i<theta_bins.length-1; i++){
	
	        wr.write( "# " + theta_bins[i] + "," + theta_bins[i+1] + "\n" );
	
	        H2F h = a_h2_z_phi.get(i);
	
	        for( int j=0; j < x.getNBins() ; j++ ){
	          for( int k=0; k < y.getNBins()-1 ; k++ ){
	            wr.write( h.getBinContent(j,k) + "," );
	          }
	          wr.write( h.getBinContent(j, y.getNBins()-1) + "\n" );
	        }
	      }
	      wr.close();
	    } catch ( IOException e ) {
	    }
	  }
	
	// read 2D histograms from TXT
	public void readHistogramsFromTXT(String filename){
		
		double[] theta_bins = bs.getThetaBins();
		ArrayList<H2F> a_h2_z_phi =  bs.getA_h2_z_phi();
		
	    try {
	      System.out.println("Reading from: "+filename+" ...");
	      BufferedReader br = new BufferedReader( new FileReader(filename) );

	      H2F h = null;
	      String line = "";

	      int i=0; // x bin counter

	      while ( ( line = br.readLine() ) != null ){
	        String[] ll = line.split(",");
	        if( ll.length == 2 ){
	          float fl = Float.parseFloat( ll[1] );
	          int bin = Arrays.binarySearch( theta_bins, fl );
	          h = a_h2_z_phi.get(bin-1);
	          i = 0; // reset x counter
	        }
	        if( line.startsWith("#") == false ){
	          //if( ll.length != h.getYAxis().getNBins() ) {
	          //  System.out.println( i + " error " + ll  ); 
	          //}

	          for( int j=0; j<h.getYAxis().getNBins();j++ ){
	            float f = Float.parseFloat(ll[j]);
	            h.setBinContent(i,j, f + h.getBinContent(i,j) );
	          }
	          i++;
	        }

	      }      

	    } catch ( IOException e ) {
	      System.err.println( e );
	      System.exit(-1);
	    }
	    
	    bs.setA_h2_z_phi(a_h2_z_phi);
	  }

	  public void readHistogramsFromTXT(List<String> filenames) {
	    for (String f : filenames) readHistogramsFromTXT(f);
	  }
}
