package modules;

import modules.CDbeamSpot;
import modules.DCbeamSpot;
import objects.Particle;

public class BeamSpot {

	public static void main(String[] args) {
		DCbeamSpot DC = new DCbeamSpot("DC");
		CDbeamSpot CD = new CDbeamSpot("CD");
		//Particle part = new Particle("part");
		
		DC.Print();
		CD.Print();
		//part.Print();

	}

}
