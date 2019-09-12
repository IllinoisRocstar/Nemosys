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
MAdLib_DIRS = @MAdLib_DIRS@
MAdLib_DEFS = @MAdLib_DEFS@
MAdLib_INCLUDES = @MAdLib_INCLUDES@
MAdLib_BENCHDIRS = @MAdLib_BENCHDIRS@
MAdLib_LIBADD = @MAdLib_LIBADD@
VPATH = @srcdir@
DOXYGEN = @DOXYGEN@
docdir = @docdir@
includedir = @includedir@
allincludedir = @allincludedir@

MADLIB_MAJOR_VERSION = 1
MADLIB_MINOR_VERSION = 3
MADLIB_PATCH_VERSION = 0

MADLIB_VERSION = ${MADLIB_MAJOR_VERSION}.${MADLIB_MINOR_VERSION}.${MADLIB_PATCH_VERSION}

SUBDIRS = ${MAdLib_DIRS}
AUTOMAKE_OPTIONS = subdir-objects

DISTCHECK_CONFIGURE_FLAGS = --enable-benchmarks --enable-ann

dist_doc_DATA = README License.txt Copyright.txt Credits.txt Tutorial/MAdLibInterface.h Tutorial/MAdLibInterface.cpp Tutorial/Makefile doxygen.config

EXTRA_DIST = \
Benchmarks/moveIt/examples/tube/tube.geo \
Benchmarks/moveIt/examples/tube/MyParams.h \
Benchmarks/moveIt/examples/tube/tube.mad \
Benchmarks/moveIt/examples/square/square.geo \
Benchmarks/moveIt/examples/square/square.mad \
Benchmarks/moveIt/examples/cube/cube.geo \
Benchmarks/moveIt/examples/cube/cube.mad \
Benchmarks/moveIt/examples/sphere/spheres.geo \
Benchmarks/moveIt/examples/sphere/spheres.mad \
Benchmarks/moveIt/examples/block/block.step \
Benchmarks/moveIt/examples/block/block.mad \
Benchmarks/moveIt/examples/A319/A319.brep \
Benchmarks/moveIt/examples/A319/A319_nobox.geo \
Benchmarks/moveIt/examples/A319/A319_2D.mad

lib_LTLIBRARIES = libMAdLib.la

libMAdLib_la_LIBADD = ${MAdLib_LIBADD}

libMAdLib_la_DEPENDENCIES = ${MAdLib_LIBADD}

CPPFLAGS += ${MAdLib_DEFS}

AM_CPPFLAGS =  ${MAdLib_INCLUDES}\
       -I$(srcdir)/Geo\
       -I$(srcdir)/Mesh\
       -I$(srcdir)/Common\
       -I$(srcdir)/Adapt\
       -I$(srcdir)/Adapt/constraint\
       -I$(srcdir)/Adapt/operator\
       -I$(srcdir)/Adapt/output\
       -I$(srcdir)/Adapt/quality\
       -I$(srcdir)/Adapt/repositioning\
       -I$(srcdir)/Adapt/sizeField\
       -I$(srcdir)/Adapt/utils


libMAdLib_la_SOURCES =

include_HEADERS = \
  Adapt/AdaptInterface.h\
  Adapt/output/MAdOutput.h\
  Adapt/repositioning/MobileObject.h\
  Adapt/operator/MAdOperations.h\
  Adapt/sizeField/LocalSizeField.h\
  Adapt/sizeField/AnalyticalSField.h\
  Adapt/sizeField/PWLinearSField.h\
  Adapt/sizeField/SizeFieldBase.h\
  Adapt/sizeField/DiscreteSF.h\
  Adapt/sizeField/SizeField.h\
  Adapt/sizeField/BackgroundSF.h\
  Adapt/utils/DistanceToPoints.h\
  Adapt/utils/DistanceFunction.h\
  Adapt/utils/CallbackDefinition.h\
  Geo/ModelInterface.h\
  Mesh/CheckMesh.h\
  Mesh/MeshDataBaseInterface.h\
  Mesh/MeshDataBaseParallelInterface.h\
  Mesh/MeshDataBaseComm.h\
  Mesh/MeshDataBaseCommPeriodic.h\
  Common/MAdLib.h

if Extended_API
allinclude_HEADERS = \
  Adapt/AdaptInterface.h\
  Adapt/constraint/ModelConstraintManager.h\
  Adapt/constraint/Constraint.h\
  Adapt/MeshParametersManager.h\
  Adapt/operator/MAdOperations.h\
  Adapt/operator/VertexMoveOp.h\
  Adapt/operator/SliverRegionHandler.h\
  Adapt/operator/FaceSwapOp.h\
  Adapt/operator/DESCOp.h\
  Adapt/operator/EdgeCollapseOp.h\
  Adapt/operator/EdgeSplitOp.h\
  Adapt/operator/EdgeSwapConfig.h\
  Adapt/operator/EdgeSwapOp.h\
  Adapt/operator/FaceCollapseOp.h\
  Adapt/operator/OperatorTools.h\
  Adapt/operator/SliverFaceHandler.h\
  Adapt/operator/RegionRemoveOp.h\
  Adapt/operator/MAdOperatorBase.h\
  Adapt/output/MAdOutput.h\
  Adapt/quality/ElementEvaluatorBase.h\
  Adapt/quality/ElementStatistics.h\
  Adapt/quality/MeshQualityManager.h\
  Adapt/quality/MeanRatioEvaluator.h\
  Adapt/quality/OrientedMeanRatioEvaluator.h\
  Adapt/repositioning/NodesRepositioningOp.h\
  Adapt/repositioning/GeoMatcher.h\
  Adapt/repositioning/MobileObject.h\
  Adapt/repositioning/LaplaceSmoothingOp.h\
  Adapt/repositioning/MAdElasticityOp.h\
  Adapt/sizeField/NullSField.h\
  Adapt/sizeField/MeshSizeBase.h\
  Adapt/sizeField/IsoMeshSize.h\
  Adapt/sizeField/LocalSizeField.h\
  Adapt/sizeField/AnalyticalSField.h\
  Adapt/sizeField/PWLinearSField.h\
  Adapt/sizeField/AnisoMeshSize.h\
  Adapt/sizeField/SizeField.h\
  Adapt/sizeField/BackgroundSF.h\
  Adapt/sizeField/SizeFieldBase.h\
  Adapt/sizeField/SizeFieldManager.h\
  Adapt/sizeField/DiscreteSF.h\
  Adapt/utils/MAdStatistics.h\
  Adapt/utils/MAdLinearSystemSparskit.h\
  Adapt/utils/DistanceToPoints.h\
  Adapt/utils/DistanceFunction.h\
  Adapt/utils/MAdLinearSystemPETSc.h\
  Adapt/utils/MAdTimeManager.h\
  Adapt/utils/MAdLinearSystemGmm.h\
  Adapt/utils/NodalDataManager.h\
  Adapt/utils/CallbackManager.h\
  Adapt/utils/CallbackDefinition.h\
  Adapt/utils/MAdLinearSystem.h\
  Adapt/utils/History.h\
  Geo/GM_Iterators.h\
  Geo/GmshEntities.h\
  Geo/GmshModel.h\
  Geo/MAdModel.h\
  Geo/ModelInterface.h\
  Geo/NullEntities.h\
  Geo/NullModel.h\
  Geo/PGList.h\
  Geo/Physical.h\
  Mesh/CheckMesh.h\
  Mesh/CheckOrientation.h\
  Mesh/Mark.h\
  Mesh/MeshDataBaseAttachable.h\
  Mesh/MeshDataBaseCommCheck.h\
  Mesh/MeshDataBaseComm.h\
  Mesh/MeshDataBaseCommPeriodic.h\
  Mesh/MeshDataBaseGEntity2Physical.h\
  Mesh/MeshDataBase.h\
  Mesh/MeshDataBaseInterface.h\
  Mesh/MeshDataBaseIO.h\
  Mesh/MeshDataBaseIterators.h\
  Mesh/MeshDataBaseLoadBalancing.h\
  Mesh/MeshDataBaseMiniMesh.h\
  Mesh/MeshDataBaseParallelInterface.h\
  Mesh/MeshDataBaseParallelIO.h\
  Mesh/metisAdaptiveRepart.h\
  Mesh/MshTags.h\
  Mesh/ParallelUtils.h\
  Mesh/PList.h\
  Common/MAdDefines.h\
  Common/MAdFieldEvaluator.h\
  Common/MAdMatrix.h\
  Common/MAdMessage.h\
  Common/MAdResourceManager.h\
  Common/MAdSingleton.h\
  Common/MAdStringFieldEvaluator.h\
  Common/MathUtils.h\
  Common/MAdLib.h\
  Common/MAdMetric.h\
  Common/MAdVector3.h\
  Contrib/mathex/mathex.h
endif

ACLOCAL_AMFLAGS = -I m4


##bench: install
##	@for i in ${MAdLib_BENCHDIRS}; do (cd $$i && ${MAKE}); done

##install-bench: bench
##	@for i in ${MAdLib_BENCHDIRS}; do (cd $$i && ${MAKE} install); done

##uninstall-bench:
##	@for i in ${MAdLib_BENCHDIRS}; do (cd $$i && ${MAKE} uninstall); done

doc:
	cd ${srcdir} && ${DOXYGEN} && cd -

install-doc: doc
	cd ${srcdir} && cp -r doc/ ${prefix}/. && cd -

##purge:
##	for i in ${MAdLib_DIRS}; do (cd $$i && ${MAKE} purge); done
##	for i in ${MAdLib_BENCHDIRS}; do (cd $$i && ${MAKE} purge); done

## GC warning: should be replaced by distcheck or sthing like that
source-tree: ##purge
	rm -rf MAdLib-${MADLIB_VERSION}
	cd ${srcdir} && tar zcf MAdLib.tgz --exclude "*.o" --exclude "*.a" --exclude "libMAdLib*"\
          --exclude "variables" --exclude "config.log" --exclude "config.status"\
          --exclude "autom4*" --exclude "Makefile.distrib" --exclude "*.bak"\
          --exclude "HTML" --exclude "*TAGS*" --exclude ".svn" --exclude "*.pos"\
          --exclude "*.msh" --exclude "toPort" --exclude "poubelle"\
          --exclude "Makefile.buildconfig"\
          Geo Mesh Adapt Common Copyright.txt License.txt Credits.txt\
          README *.in Contrib Makefile* configure Tutorial
	mkdir MAdLib-${MADLIB_VERSION}
	cd MAdLib-${MADLIB_VERSION} && tar zxf ${srcdir}/MAdLib.tgz
	rm -f ${srcdir}/MAdLib.tgz
	cd ${srcdir} && tar zcf MAdLib_Benchmarks.tgz\
          --exclude "*.bak"\
          --exclude "*.msh"\
          --exclude "*.pos"\
          --exclude ".svn"\
          --exclude "cpu" --exclude "*tmp*" --exclude "statistics"\
          --exclude "slivers" --exclude "cpuinfo" --exclude "journal" --exclude "meshSize"\
          Benchmarks/checkMesh Benchmarks/meshInfo Benchmarks/optimize\
          Benchmarks/moveIt/*.h Benchmarks/moveIt/*.cc Benchmarks/moveIt/Makefile*\
          Benchmarks/moveIt/README\
          Benchmarks/moveIt/examples/tube
	cd MAdLib-${MADLIB_VERSION} && tar zxf ${srcdir}/MAdLib_Benchmarks.tgz
	rm -f ${srcdir}/MAdLib_Benchmarks.tgz

## GC warning: should be replaced by distcheck or sthing like that
source: source-tree
	cd MAdLib-${MADLIB_VERSION} && rm -rf ${MADLIB_VERSION_FILE}\
          Common/MAdLibConfig.h 
	tar zcf MAdLib-${MADLIB_VERSION}-source.tgz MAdLib-${MADLIB_VERSION}

MADLIB_PCFILE = MAdLib.pc
pc:
	echo "prefix=${prefix}" > ${MADLIB_PCFILE}
	echo "exec_prefix=${prefix}/bin" >> ${MADLIB_PCFILE}
	echo "libdir=${libdir}" >> ${MADLIB_PCFILE}
	echo "includedir=${includedir}" >> ${MADLIB_PCFILE}
	echo "">> ${MADLIB_PCFILE}
	echo "Name: MAdLib" >> ${MADLIB_PCFILE}
	echo "Version: ${MADLIB_VERSION}" >> ${MADLIB_PCFILE}
	echo "Description: Mesh adaptation library" >> ${MADLIB_PCFILE}
	echo "Libs: -L${libdir} -lMAdLib" >> ${MADLIB_PCFILE}
	echo "Cflags: -I${includedir}" >> ${MADLIB_PCFILE}
	mkdir ${libdir}/pkgconfig
	cp ${MADLIB_PCFILE} ${libdir}/pkgconfig/.
