package objects;

import analysis.Constants;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */

public class Event {
    private final boolean debug = false;
    
    private int run;
    private int event;
    private double startTime;

    private final List<Track> particles  = new ArrayList<>();
    private final List<Track> tracks     = new ArrayList<>();
    private final List<Track> utracks     = new ArrayList<>();
  
    private DataEvent hipoEvent;

    public Event(DataEvent event) {
        this.hipoEvent = event;
        this.readEvent(event);
        if(debug) System.out.println("Read event with " + tracks.size() + " particles");
    }
    


    private DataBank getBank(DataEvent de, String bankName) {
        DataBank bank = null;
        if (de.hasBank(bankName)) {
            bank = de.getBank(bankName);
        }
        return bank;
    }

    private DataBank getBank(DataEvent de, String bankName1, String bankName2) {
        DataBank bank = null;
        if (de.hasBank(bankName1)) {
            bank = de.getBank(bankName1);
        }
        else if (de.hasBank(bankName2)) {
            bank = de.getBank(bankName2);
        }
        return bank;
    }

    private void readHeader(DataEvent event) {
        DataBank head = this.getBank(event, "RUN::config");
        if(head!=null) {
            this.run   = head.getInt("run", 0);
            this.event = head.getInt("event", 0);
        }
    }
    
    private void readParticles(DataEvent event) {
        DataBank recPart   = this.getBank(event, "REC::Particle");
        DataBank recTrack  = this.getBank(event, "REC::Track");
        DataBank runConfig = this.getBank(event, "RUN::config");
        if(recPart!=null && recTrack!=null) {
            for (int i = 0; i < recPart.rows(); i++) {    
                Track track = Track.readParticle(recPart, recTrack, i);
                if(runConfig!=null) track.addScale(runConfig);
                particles.add(track);
            }
        }
    }
    
    private void readStartTime(DataEvent event) {
        DataBank recEvent = this.getBank(event, "REC::Event");
        if(recEvent!=null) {
            startTime = recEvent.getFloat("startTime",0);
        }
    }
    
    private void readTracks(DataEvent event) {
        DataBank recPart   = this.getBank(event, "REC::Particle");
        DataBank recTrack  = this.getBank(event, "REC::Track");
        DataBank urecTrack = this.getBank(event, "REC::UTrack");
        DataBank runConfig = this.getBank(event, "RUN::config");
       
          if(recPart!=null && recTrack!=null) {
            for (int i = 0; i < recPart.rows(); i++) {    
                Track track = Track.readParticle(recPart, recTrack, i);
                if(runConfig!=null) track.addScale(runConfig);
                tracks.add(track);
                System.out.println(recPart.getInt("pid", row));
            }            
            if(urecTrack!=null) {
                for (int i = 0; i < recPart.rows(); i++) {    
                    Track track = Track.readParticle(recPart, recTrack, urecTrack, i);
                    if(runConfig!=null) track.addScale(runConfig);
                    utracks.add(track);
                }            
            }
        }
    }
    
    private void readEvent(DataEvent de) {
        this.readHeader(de);
        this.readParticles(de);
        this.readStartTime(de);
        this.readTracks(de);
    }

    public List<Track> getTracks() {
        return tracks;
    }
 
    public List<Track> getUTracks() {
        return utracks;
    }
 
   public List<Track> getParticles() {
        return particles;
    }

    public int getRun() {
        return run;
    }

    public int getEvent() {
        return event;
    }
    
    public double getStartTime() {
        return startTime;
    }

    public DataEvent getHipoEvent() {
        return hipoEvent;
    }
    
    
}
