#!/bin/sh
#
# jemboss/runJemboss.sh.  Generated from runJemboss.sh.in by configure.
# EMBOSS
# Version: $Revision: 1.2 $
# Modified $Date: 2012/07/23 13:47:14 $ by $Author: rice $

export JEMBOSS_HOME="/home/aap599/software/JAMg/3rd_party/EMBOSS/../share/EMBOSS/jemboss";

# Override an existing CLASSPATH environment variable

export CLASSPATH="";

CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/activation.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/client.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/jakarta-regexp-1.2.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/jalviewApplet.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/jcert.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/jemboss.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/jnet.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/jsse.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/mail.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/axis.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/commons-discovery.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/commons-logging.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/jaxrpc.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/log4j-1.2.4.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/saaj.jar";
CLASSPATH="${CLASSPATH}:${JEMBOSS_HOME}/lib/axis/wsdl4j.jar";

cd "${JEMBOSS_HOME}";

# Add local to run Jemboss in 'standalone' mode:
# java org.emboss.jemboss.Jemboss local &

case "${1}" in
  local)
    java org/emboss/jemboss/Jemboss local &
    ;;
  *)
    java org/emboss/jemboss/Jemboss &
    ;;
esac
