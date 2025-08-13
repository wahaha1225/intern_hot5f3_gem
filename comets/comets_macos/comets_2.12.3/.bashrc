#
# default .bashrc
# 03/31/13
#
# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi
 
umask 022


export JAVA_HOME=/usr/java/default/


export COMETS_HOME=/projectnb/qedksk/comets/comets_2.12.3
export PATH=$PATH:$COMETS_HOME
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$COMETS_HOME/lib/or-tools
export PYTHONPATH=$PYTHONPATH:$COMETS_HOME/lib/cometspy








