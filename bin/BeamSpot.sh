#!/bin/sh

SCRIPT_DIR=`dirname $0`

java -Dsun.java2d.pmoffscreen=false -Xmx2048m -Xms1024m -Xdiag -cp "$SCRIPT_DIR/../target/BeamSpotSimple-0.0.1-SNAPSHOT-jar-with-dependencies.jar" analysis.BeamSpot $*
