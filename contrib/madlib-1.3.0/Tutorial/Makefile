# -------------------------------------------------------------------
# MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
#
# See the Copyright.txt and License.txt files for license information. 
# You should have received a copy of these files along with MAdLib. 
# If not, see <http://www.madlib.be/license/>
#
# Please report all bugs and problems to <contrib@madlib.be>
#
# Authors: Gaetan Compere, Jean-Francois Remacle
# -------------------------------------------------------------------


INC =  ${MAdLib_INCLUDES}\
       ${DASH}I$(MAdROOT)/Tutorial\
       ${DASH}I$(MAdROOT)/Geo\
       ${DASH}I$(MAdROOT)/Mesh\
       ${DASH}I$(MAdROOT)/Common\
       ${DASH}I$(MAdROOT)/Adapt\
       ${DASH}I$(MAdROOT)/Adapt/constraint\
       ${DASH}I$(MAdROOT)/Adapt/operator\
       ${DASH}I$(MAdROOT)/Adapt/output\
       ${DASH}I$(MAdROOT)/Adapt/quality\
       ${DASH}I$(MAdROOT)/Adapt/repositioning\
       ${DASH}I$(MAdROOT)/Adapt/sizeField\
       ${DASH}I$(MAdROOT)/Adapt/utils

CXXFLAGS = ${OPTIM} ${MAdLib_DEFS} ${FLAGS} ${INC} ${SYSINCLUDE}

SRC = MAdLibInterface.cc

OBJ = ${SRC:.cc=${OBJEXT}}

.SUFFIXES: ${OBJEXT} .cc

.cc${OBJEXT}:
	${CXX} ${CXXFLAGS} ${DASH}c $< ${DASH}o $@

build: ${OBJ}

clean:
	${RM} */*.o *.o *.obj

purge:

depend:
	(sed '/^# DO NOT DELETE THIS LINE/q' Makefile && \
         ${CXX} -MM ${CXXFLAGS} ${SRC} | sed 's/.o:/$${OBJEXT}:/g' \
        ) > Makefile.new
	cp Makefile Makefile.bak
	cp Makefile.new Makefile
	rm -f Makefile.new