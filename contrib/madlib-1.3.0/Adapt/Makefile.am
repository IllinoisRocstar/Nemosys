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
VPATH = @srcdir@

noinst_LTLIBRARIES = libMAdAdapt.la

CPPFLAGS += ${MAdLib_DEFS}

AM_CPPFLAGS =  ${MAdLib_INCLUDES}\
  -I$(srcdir)/../Common\
  -I$(srcdir)/../Geo\
  -I$(srcdir)/../Mesh\
  -I${srcdir}/../Contrib/mathex\
  -I$(srcdir)/constraint\
  -I$(srcdir)/operator\
  -I$(srcdir)/output\
  -I$(srcdir)/quality\
  -I$(srcdir)/repositioning\
  -I$(srcdir)/sizeField\
  -I$(srcdir)/utils

libMAdAdapt_la_SOURCES = \
  ../Contrib/mathex/mathex.cpp\
  AdaptInterface.cc\
  MeshParametersManager.cc\
  constraint/Constraint.cc\
  constraint/ModelConstraintManager.cc\
  operator/DESCOp.cc\
  operator/EdgeCollapseOp.cc\
  operator/EdgeSplitOp.cc\
  operator/EdgeSwapConfig.cc\
  operator/EdgeSwapOp.cc\
  operator/FaceCollapseOp.cc\
  operator/FaceSwapOp.cc\
  operator/MAdOperatorBase.cc\
  operator/OperatorTools.cc\
  operator/RegionRemoveOp.cc\
  operator/SliverFaceHandler.cc\
  operator/SliverRegionHandler.cc\
  operator/VertexMoveOp.cc\
  output/MAdOutput.cc\
  quality/MeanRatioEvaluator.cc\
  quality/OrientedMeanRatioEvaluator.cc\
  quality/MeshQualityManager.cc\
  repositioning/GeoMatcher.cc\
  repositioning/LaplaceSmoothingOp.cc\
  repositioning/MAdElasticityOp.cc\
  repositioning/MobileObject.cc\
  sizeField/AnalyticalSField.cc\
  sizeField/AnisoMeshSize.cc\
  sizeField/BackgroundSF.cc\
  sizeField/DiscreteSF.cc\
  sizeField/IsoMeshSize.cc\
  sizeField/LocalSizeField.cc\
  sizeField/MeshSizeBase.cc\
  sizeField/NullSField.cc\
  sizeField/PWLinearSField.cc\
  sizeField/SizeFieldBase.cc\
  sizeField/SizeFieldManager.cc\
  utils/CallbackManager.cc\
  utils/History.cc\
  utils/MAdStatistics.cc\
  utils/MAdTimeManager.cc\
  utils/NodalDataManager.cc\
  utils/DistanceFunction.cc

purge:
