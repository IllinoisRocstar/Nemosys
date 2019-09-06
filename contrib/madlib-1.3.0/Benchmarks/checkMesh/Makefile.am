## -------------------------------------------------------------------
## MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain

## See the Copyright.txt and License.txt files for license information. 
## You should have received a copy of these files along with MAdLib. 
## If not, see <http://www.madlib.be/license/>

## Please report all bugs and problems to <contrib@madlib.be>

## Author: Gaetan Compere
## -------------------------------------------------------------------

## Process this file with automake to produce Makefile.in

includedir = @includedir@
MAdLib_DEFS = @MAdLib_DEFS@
MAdLib_LIBS = @MAdLib_LIBS@
MAdLib_INCLUDES = @MAdLib_INCLUDES@
RM = @RM@

bin_PROGRAMS = checkMesh

checkMesh_CPPFLAGS = ${MAdLib_DEFS} ${MAdLib_INCLUDES}\
  -I$(includedir)

checkMesh_SOURCES = main.cc

checkMesh_LDADD = -L$(top_builddir) ${MAdLib_LIBS}

purge:
	${RM} *.msh *.geo *.pos
