package objects;

import org.jlab.groot.math.Func1D;

public class FitFunc extends Func1D {

    public FitFunc( String name, double min, double max ) {
      super( name, min, max );
      addParameter( "z_0" );
      addParameter( "A" );
      addParameter( "#phi_0" );
    }

    @Override
    public double evaluate( double x ) {
      return this.getParameter(0) - this.getParameter(1) * Math.cos( x * Math.PI/180.0 - this.getParameter(2) );
    }
}
