## -------------------------------------------------------------------
## MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain

## See the Copyright.txt and License.txt files for license information. 
## You should have received a copy of these files along with MAdLib. 
## If not, see <http://www.madlib.be/license/>

## Please report all bugs and problems to <contrib@madlib.be>

## Author: Gaetan Compere
## -------------------------------------------------------------------

## Process this file with automake to produce Makefile.in

srcdir = @srcdir@
MAdLib_DEFS = @MAdLib_DEFS@
MAdLib_INCLUDES = @MAdLib_INCLUDES@

noinst_LTLIBRARIES = libMAdGeo.la

CPPFLAGS += ${MAdLib_DEFS}

AM_CPPFLAGS =  ${MAdLib_INCLUDES} -I.\
  -I$(srcdir)/../Common\
  -I$(srcdir)/../Mesh\
  -I${srcdir}/../Contrib/mathex

libMAdGeo_la_SOURCES = \
  GmshModel.cc\
  ModelInterface.cc\
  NullModel.cc\
  GmshEntities.cc\
  NullEntities.cc\
  MAdModel.cc

purge:
