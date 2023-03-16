#!/bin/sh

SCRIPT_DIR=`dirname $0`

#java -Dsun.java2d.pmoffscreen=false -Xmx2048m -Xms1024m -Xdiag -cp "$SCRIPT_DIR/../target/BeamSpotMerged-0.0.1-jar-with-dependencies.jar" org.clas.beamspot.analysis.DCbeamSpot $*
java -jar "$SCRIPT_DIR/../target/beam-spot-0.0.1-SNAPSHOT-jar-with-dependencies.jar" $*
