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

noinst_LTLIBRARIES = libMAdMesh.la

CPPFLAGS += ${MAdLib_DEFS}

AM_CPPFLAGS =  ${MAdLib_INCLUDES} -I.\
  -I$(srcdir)/../Common\
  -I$(srcdir)/../Geo\
  -I${srcdir}/../Contrib/mathex

libMAdMesh_la_SOURCES = \
  Balance.cc\
  CheckMesh.cc\
  CheckOrientation.cc\
  Mark.cc\
  MeshDataBase.cc\
  MeshDataBaseComm.cc\
  MeshDataBaseCommCheck.cc\
  MeshDataBaseGEntity2Physical.cc\
  MeshDataBaseInterface.cc\
  MeshDataBaseIO.cc\
  MeshDataBaseLoadBalancing.cc\
  MeshDataBaseMigration.cc\
  MeshDataBaseMiniMesh.cc\
  MeshDataBaseParallelInterface.cc\
  MeshDataBaseParallelIO.cc\
  metisAdaptiveRepart.cc\
  ParallelUtils.cc\
  PeriodicInterfaceMigration.cc

purge:
