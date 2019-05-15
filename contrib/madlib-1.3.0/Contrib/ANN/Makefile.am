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
includedir = @includedir@

noinst_LTLIBRARIES = libMAdANN.la

include_HEADERS = \
  include/ANN/ANN.h\
  include/ANN/ANNperf.h\
  include/ANN/ANNx.h

noinst_HEADERS = src/bd_tree.h\
src/kd_fix_rad_search.h\
src/kd_pr_search.h\
src/kd_search.h\
src/kd_split.h\
src/kd_tree.h\
src/kd_util.h\
src/pr_queue.h\
src/pr_queue_k.h

AM_CPPFLAGS = -I${srcdir}/include

libMAdANN_la_SOURCES = \
  src/ANN.cpp\
  src/bd_fix_rad_search.cpp\
  src/bd_pr_search.cpp\
  src/bd_search.cpp\
  src/bd_tree.cpp\
  src/brute.cpp\
  src/kd_dump.cpp\
  src/kd_fix_rad_search.cpp\
  src/kd_pr_search.cpp\
  src/kd_search.cpp\
  src/kd_split.cpp\
  src/kd_tree.cpp\
  src/kd_util.cpp\
  src/perf.cpp

purge:
