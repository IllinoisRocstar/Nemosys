dnl -------------------------------------------------------------------
dnl MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain

dnl See the Copyright.txt and License.txt files for license information. 
dnl You should have received a copy of these files along with MAdLib. 
dnl If not, see <http://www.madlib.be/license/>

dnl Please report all bugs and problems to <contrib@madlib.be>

dnl Authors: Gaetan Compere, Jean-Francois Remacle
dnl -------------------------------------------------------------------

dnl Process this file with autoconf to produce the configure script.

dnl Check that this is the MAdLib source tree
dnl AC_INIT(Mesh/MeshDataBaseInterface.h)
AC_INIT([MAdLib], [1.3.0], [contrib@madlib.be])

dnl Initialize Automake
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

dnl Initialize Libtool
LT_INIT()

dnl Initialize Doxygen
dnl DX_INIT_DOXYGEN([MAdLib], [doxygen.config], [doc/])

dnl Parse '--with' command-line options
AC_ARG_WITH(blas-lapack-prefix,
            AC_HELP_STRING([--with-blas-lapack-prefix=PFX],
                           [prefix where BLAS and LAPACK are installed]),
            [BLAS_LAPACK_PREFIX=$withval])
dnl AC_ARG_WITH(gsl-prefix,
dnl             AC_HELP_STRING([--with-gsl-prefix=PFX],
dnl                            [prefix where the GSL is installed]),
dnl             [GSL_PREFIX=$withval])
AC_ARG_WITH(mpi-prefix,
            AC_HELP_STRING([--with-mpi-prefix=PFX],
                           [prefix where MPI is installed]),
            [MPI_PREFIX=$withval])
AC_ARG_WITH(autopack-prefix,
            AC_HELP_STRING([--with-autopack-prefix=PFX],
                           [prefix where Autopack is installed]),
            [AUTOPACK_PREFIX=$withval])
AC_ARG_WITH(metis-prefix,
            AC_HELP_STRING([--with-metis-prefix=PFX],
                           [prefix where Metis is installed]),
            [METIS_PREFIX=$withval])
AC_ARG_WITH(parmetis-prefix,
            AC_HELP_STRING([--with-parmetis-prefix=PFX],
                           [prefix where ParMetis is installed]),
            [PARMETIS_PREFIX=$withval])
AC_ARG_WITH(gmsh-prefix,
            AC_HELP_STRING([--with-gmsh-prefix=PFX],
                           [prefix where Gmsh is installed]),
            [GMSH_PREFIX=$withval])
AC_ARG_WITH(occ-prefix,
            AC_HELP_STRING([--with-occ-prefix=PFX],
                           [prefix where OpenCascade is installed]),
            [OCC_PREFIX=$withval])
AC_ARG_WITH(sparskit-prefix,
            AC_HELP_STRING([--with-sparskit-prefix=PFX],
                           [prefix where Sparskit is installed]),
            [SPARSKIT_PREFIX=$withval])
AC_ARG_WITH(petsc-prefix,
            AC_HELP_STRING([--with-petsc-prefix=PFX],
                           [prefix where PETSc is installed]),
            [PETSC_PREFIX=$withval])
AC_ARG_WITH(parser-prefix,
            AC_HELP_STRING([--with-parser-prefix=PFX],
                           [prefix where the parser is installed]),
            [PARSER_PREFIX=$withval])

dnl Parse '--enable' command line options
AC_ARG_ENABLE(blas-lapack,
              AC_HELP_STRING([--enable-blas-lapack],
                             [use Blas/Lapack for linear algebra (default=yes)]))
AC_ARG_ENABLE(ann,
              AC_HELP_STRING([--enable-ann],
                             [use ANN library (default=no)]))
AC_ARG_ENABLE(mathex,
              AC_HELP_STRING([--enable-mathex],
                             [use Mathex library (default=yes)]))
dnl AC_ARG_ENABLE(gsl,
dnl               AC_HELP_STRING([--enable-gsl],
dnl                              [use GSL as numerical toolkit (default=no)]))
AC_ARG_ENABLE(mpi,
              AC_HELP_STRING([--enable-mpi],
                             [enable MPI support (default=no)]))
AC_ARG_ENABLE(metis,
              AC_HELP_STRING([--enable-metis],
                             [use Metis partitioner (default=no)]))
AC_ARG_ENABLE(parmetis,
              AC_HELP_STRING([--enable-parmetis],
                             [use ParMetis partitioner (default=no)]))
AC_ARG_ENABLE(gmm,
              AC_HELP_STRING([--enable-gmm],
                             [compile gmm++ linear solvers (default=no)]))
AC_ARG_ENABLE(gmsh,
              AC_HELP_STRING([--enable-gmsh],
                             [use Gmsh geometric model (default=no)]))
AC_ARG_ENABLE(occ,
              AC_HELP_STRING([--enable-occ],
                             [use OpenCascade geometric model (through Gmsh) (default=no)]))
AC_ARG_ENABLE(petsc,
              AC_HELP_STRING([--enable-petsc],
                             [use PETSc if available (default=no)]))
AC_ARG_ENABLE(parser,
              AC_HELP_STRING([--enable-parser],
                             [use Parser (restricted, CENAERO, Belgium) if available (default=no)]))
AC_ARG_ENABLE(benchmarks,
              AC_HELP_STRING([--enable-benchmarks],
                             [Compile all benchmark executables: checkMesh, optimize, meshInfo, moveIt (default=no)]))
AC_ARG_ENABLE(moveIt,
              AC_HELP_STRING([--enable-moveIt],
                             [Compile moveIt, the main benchmark executable (default=no)]))
AC_ARG_ENABLE(API-extended,
              AC_HELP_STRING([--enable-extendedAPI],
                             [Whether all header files are included in the API or not (default=no)]))

dnl Disable options when disabled by default and not marked as enabled
if test "x$enable_ann" != "xyes"; then
  enable_ann=no;
fi
dnl if test "x$enable_gsl" != "xyes"; then
dnl   enable_gsl=no;
dnl fi
if test "x$enable_mpi" != "xyes"; then
  enable_mpi=no;
fi
if test "x$enable_metis" != "xyes"; then
  enable_metis=no;
fi
if test "x$enable_parmetis" != "xyes"; then
  enable_parmetis=no;
fi
if test "x$enable_gmm" != "xyes"; then
  enable_gmm=no;
fi
if test "x$enable_gmsh" != "xyes"; then
  enable_gmsh=no;
  enable_occ=no;
fi
if test "x$enable_gmsh" != "xno"; then
  enable_gmm=yes;
fi
if test "x$enable_occ" != "xyes"; then
  enable_occ=no;
fi
if test "x$enable_petsc" != "xyes"; then
  enable_petsc=no;
fi
if test "x$enable_parser" != "xyes"; then
  enable_parser=no;
fi
if test "x$enable_benchmarks" != "xyes"; then
  enable_benchmarks=no;
fi
if test "x$enable_moveIt" != "xyes"; then
  enable_moveIt=no;
fi
if test "x$enable_extendedAPI" != "xyes"; then
  enable_extendedAPI=no;
fi

dnl Set the auxiliary compilation flag to empty
AC_DEFINE([HAVE_64BIT_SIZE_T], [], [Whether we run the compilation on a 64 bits architecture or not])
AC_DEFINE([HAVE_BLAS], [], [Whether we use BLAS or not])
AC_DEFINE([HAVE_LAPACK], [], [Whether we use Lapack or not])
AC_DEFINE([HAVE_NO_DLL], [], [Whether we use DLLs or not])
AC_DEFINE([HAVE_NO_SOCKLEN_T], [], [Whether we use socklens or not])
AC_DEFINE([PARALLEL], [], [Whether we compile for parallel utilization or not])
AC_DEFINE([_HAVE_ANN_], [], [Whether we use the ANN library or not])
AC_DEFINE([_HAVE_BLAS_], [], [Whether we use BLAS or not])
AC_DEFINE([_HAVE_GMM_], [], [Whether we use GMM or not])
AC_DEFINE([_HAVE_GMSH_], [], [Whether we use Gmsh or not])
AC_DEFINE([_HAVE_MATHEX_], [], [Whether we use Mathex or not])
AC_DEFINE([_HAVE_METIS_], [], [Whether we use Metis or not])
AC_DEFINE([_HAVE_MPI_], [], [Whether we use MPI or not])
AC_DEFINE([_HAVE_OCC_], [], [Whether we use OpenCascade or not])
AC_DEFINE([_HAVE_PARMETIS_], [], [Whether we use ParMetis or not])
AC_DEFINE([_HAVE_PARSER_], [], [Whether we use the Cenaero parser or not])
AC_DEFINE([_HAVE_PETSC_], [], [Whether we use PETSc or not])

dnl set a conditional for Automake to know if the API is the normal or extended one
AM_CONDITIONAL([Extended_API],[test "x${enable_extendedAPI}" = "xyes"])

dnl nasty definition: because something wrong in autoconf ?
dnl abssrcdir="../$srcdir"

dnl Get the operating system and machine names
UNAME=`uname`
HOSTNAME=`hostname`

dnl Check for default compilers
AC_PROG_CC
AC_PROG_CXX
if test "x${CC}" = "x" -o "x${CXX}" = "x" ; then
  AC_MSG_ERROR([Could not find required compilers, aborting.])
fi

dnl Change the compiler if PETSc is used
if test "x$enable_petsc" != "xno"; then
  AC_CHECK_LIB(petsc,main,PETSC="yes",[],-lpetsc)
  if test "x${PETSC}" = "xyes"; then
    CXX="mpic++.openmpi"
    LINKER="mpic++.openmpi"
  fi
fi

dnl Set preprocessor and linker
AC_PROG_CPP
LINKER="${CXX}"

dnl Set default compiler flags
CXXFLAGS="${CXXFLAGS} -D_FORTIFY_SOURCE=0 -ansi "

dnl Take care of no-cygwin option before doing any other tests
case "$UNAME" in
  CYGWIN*)
    if test "x$enable_cygwin" != "xyes"; then
      UNAME="${UNAME}-no-cygwin"
      CC="${CC} -mno-cygwin"
      CXX="${CXX} -mno-cygwin"
      LINKER="${LINKER} -mno-cygwin"
    fi
    ;;
esac

dnl Use c++ for all compilation tests
AC_LANG(C++)

dnl See if we need a .exe extension on executables
AC_EXEEXT

dnl defs for the objects compilation
MAdLib_DEFS=""

dnl include flags for the objects compilation
MAdLib_INCLUDES="       -I\$(top_srcdir)/Geo\
       -I\$(top_srcdir)/Mesh\
       -I\$(top_srcdir)/Common\
       -I\$(top_srcdir)/Adapt\
       -I\$(top_srcdir)/Adapt/constraint\
       -I\$(top_srcdir)/Adapt/operator\
       -I\$(top_srcdir)/Adapt/output\
       -I\$(top_srcdir)/Adapt/quality\
       -I\$(top_srcdir)/Adapt/repositioning\
       -I\$(top_srcdir)/Adapt/sizeField\
       -I\$(top_srcdir)/Adapt/utils"

dnl libraries included in libMAdLib
MAdLib_LIBADD="Common/libMAdCommon.la Geo/libMAdGeo.la Mesh/libMAdMesh.la Adapt/libMAdAdapt.la"

dnl librairies to be included in the executables compilation
MAdLib_LIBS="-L${libdir} -Wl,-rpath,${libdir} -lMAdLib"

dnl dirs in which the main compilation process will go
MAdLib_DIRS="Common/ Geo/ Mesh/ Adapt/ ."

dnl dirs in which the benchmark compilation process will go
benchs=""
MAdLib_BENCHDIRS=""
if test "x$enable_benchmarks" != "xno"; then
  MAdLib_BENCHDIRS="Benchmarks/checkMesh Benchmarks/meshInfo Benchmarks/optimize Benchmarks/moveIt"
  benchs="checkMesh, meshInfo, optimize, moveIt"
  MAdLib_DIRS="${MAdLib_DIRS} Benchmarks/checkMesh Benchmarks/meshInfo Benchmarks/optimize Benchmarks/moveIt"
else
  if test "x$enable_moveIt" != "xno"; then
    MAdLib_BENCHDIRS="Benchmarks/moveIt"
    benchs="moveIt"
    MAdLib_DIRS="${MAdLib_DIRS} Benchmarks/moveIt"
  else
    benchs=""
  fi
fi

dnl Check for standard math library (no rule given if found, so gets
dnl added to $LIBS, used for further checks)
AC_CHECK_LIB(m,main)

dnl Check for gmm++ linear solver
if test "x$enable_gmm" != "xno"; then
  AC_CHECK_FILE(${srcdir}/Contrib/gmm/gmm.h,GMM="yes")
  if test "x${GMM}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_GMM_"
    MAdLib_INCLUDES="${MAdLib_INCLUDES} -I\$(top_srcdir)/Contrib/gmm "
    AC_DEFINE(_HAVE_GMM_)
    BO="${BO} Gmm"
    echo "********************************************************************"
    echo "  You are building a version of MAdLib that contains the gmm++"
    echo "  linear solvers. Gmm++ is available under the GNU LGPL."
    echo "  To disable gmm++, run configure again with the --disable-gmm"
    echo "  option."
    echo "********************************************************************"
  fi
fi

dnl Check for Gmsh
if test "x$enable_gmsh" != "xno"; then
  if test "x${GMSH_PREFIX}" != "x"; then
    LDFLAGS="-L${GMSH_PREFIX}/lib ${LDFLAGS}"
  fi
  AC_CHECK_FILE("${GMSH_PREFIX}/include/gmsh/Gmsh.h",GMSH="yes") 
dnl  AC_CHECK_LIB(Gmsh,main,GMSH="yes",[],-lGmsh)
  if test "x${GMSH}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_GMSH_"
    AC_DEFINE(_HAVE_GMSH_)
    BO="${BO} Gmsh"
    if test "x${GMSH_PREFIX}" = "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lGmsh -lfltk_images -lfltk -lfltk_gl -lfltk_forms -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -llapack -lgmp"
    else
      MAdLib_LIBS="${MAdLib_LIBS} -L${GMSH_PREFIX}/lib -Wl,-rpath,${GMSH_PREFIX}/lib -lGmsh -lfltk_images -lfltk -lfltk_gl -lfltk_forms -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -llapack -lgmp"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${GMSH_PREFIX} -I${GMSH_PREFIX}/include"
    fi
    echo "********************************************************************"
    echo "  You are building a version of MAdLib that uses Gmsh and. Gmsh is "
    echo "  available under the GNU GPL."
    echo "  To disable Gmsh, run configure again with the --disable-gmsh "
    echo "  option."
dnl    echo "  Note that you should disable Chaco and Metis in Gmsh or add"
dnl    echo "  the corresponding libraries in the current version of MAdLib."
    echo "********************************************************************"
  fi
fi

dnl Check for OpenCascade
if test "x$enable_occ" != "xno"; then
  if test "x${OCC_PREFIX}" != "x"; then
    LDFLAGS="-L${OCC_PREFIX}/lib ${LDFLAGS}"
  fi
dnl  AC_CHECK_FILE("${OCC_PREFIX}/inc/Geom_Curve.hxx",OCC="yes") 
  AC_CHECK_LIB(TKernel,main,OCC="yes",[],-lTKernel)
  if test "x${OCC}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_OCC_"
    AC_DEFINE(_HAVE_OCC_)
    BO="${BO} OpenCascade"
    if test "x${OCC_PREFIX}" != "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -L${OCC_PREFIX}/lib -Wl,-rpath,${OCC_PREFIX}/lib"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${OCC_PREFIX}/include"
    fi
    # DataExchange (subset; see occ/ros/adm/make/Makefile for more info)
    MAdLib_LIBS="${MAdLib_LIBS} -lTKSTEP -lTKSTEP209 -lTKSTEPAttr -lTKSTEPBase -lTKIGES -lTKXSBase"
    # ModelingAlgorithms
    MAdLib_LIBS="${MAdLib_LIBS} -lTKOffset -lTKFeat -lTKFillet -lTKBool -lTKShHealing"
    MAdLib_LIBS="${MAdLib_LIBS} -lTKMesh -lTKHLR -lTKBO -lTKPrim -lTKTopAlgo -lTKGeomAlgo"
    # ModelingData
    MAdLib_LIBS="${MAdLib_LIBS} -lTKBRep -lTKGeomBase -lTKG3d -lTKG2d"
    # FoundationClasses
    MAdLib_LIBS="${MAdLib_LIBS} -lTKAdvTools -lTKMath -lTKernel"
  fi
fi

dnl Check for ANN, the Approximate Nearest Neighbor library
if test "x$enable_ann" != "xno"; then
  AC_CHECK_FILE(${srcdir}/Contrib/ANN/include/ANN/ANN.h,ANN="yes") 
  if test "x${ANN}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_ANN_"
    AC_DEFINE(_HAVE_ANN_)
    BO="${BO} Ann"
    MAdLib_DIRS="Contrib/ANN ${MAdLib_DIRS}"
    MAdLib_INCLUDES="${MAdLib_INCLUDES} -I\$(top_srcdir)/Contrib/ANN/include/ANN -I\$(top_srcdir)/Contrib/ANN/src"
    MAdLib_LIBADD="Contrib/ANN/libMAdANN.la ${MAdLib_LIBADD}"
    echo "********************************************************************"
    echo "  You are building a version of MAdLib that uses ANN, the"
    echo "  Approximate Nearest Neighbor library. ANN is available under the"
    echo "  GNU LGPL. To disable ANN, run configure again with the"
    echo "  --disable-ann option."
    echo "********************************************************************"
  fi
fi

dnl Check for Mathex
if test "x$enable_mathex" != "xno"; then
  AC_CHECK_FILE(${srcdir}/Contrib/mathex/mathex.h,MATHEX="yes")
  if test "x${MATHEX}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_MATHEX_"
    MAdLib_INCLUDES="${MAdLib_INCLUDES} -I\$(top_srcdir)/Contrib/mathex "
    AC_DEFINE(_HAVE_MATHEX_)
    BO="${BO} Mathex"
    echo "********************************************************************"
    echo "  You are building a version of MAdLib that contains Mathex."
    echo "  Mathex is available under the GNU LGPL. To disable Mathex, run "
    echo "  configure again with the --disable-mathex option."
    echo "********************************************************************"
  fi
fi

dnl Check for PETSc solvers
if test "x$enable_petsc" != "xno"; then
  if test "x${PETSC_PREFIX}" != "x"; then
    LDFLAGS="-L${PETSC_PREFIX}/lib ${LDFLAGS}"
  fi
  AC_CHECK_LIB(petsc,main,PETSC="yes",[],-lpetsc)
  if test "x${PETSC}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_PETSC_"
    AC_DEFINE(_HAVE_PETSC_)
    BO="${BO} PETSc"
    if test "x${PETSC_PREFIX}" = "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lpetsc -lpetscksp -lpetscmat -lpetscvec "
    else
      MAdLib_LIBS="${MAdLib_LIBS} -L${PETSC_PREFIX}/lib -Wl,-rpath,${PETSC_PREFIX}/lib -lpetsc -lpetscksp -lpetscmat -lpetscvec "
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${PETSC_PREFIX} -I${PETSC_PREFIX}/include"
    fi
  fi
fi

dnl Check for Metis partitionner
if test "x$enable_metis" != "xno"; then
  if test "x${METIS_PREFIX}" != "x"; then
    LDFLAGS="-L${METIS_PREFIX} ${LDFLAGS}"
  fi
  AC_CHECK_LIB(metis,main,METIS="yes",[],-lmetis)
  if test "x${METIS}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_METIS_"
    AC_DEFINE(_HAVE_METIS_)
    BO="${BO} Metis"
    if test "x${METIS_PREFIX}" = "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lmetis"
    else
      MAdLib_LIBS="${MAdLib_LIBS} -L${METIS_PREFIX} -Wl,-rpath,${METIS_PREFIX} -lmetis"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${METIS_PREFIX}/Lib"
    fi
  fi
fi

dnl Check for ParMetis partitionner
if test "x$enable_parmetis" != "xno"; then
  if test "x${PARMETIS_PREFIX}" != "x"; then
    LDFLAGS="-L${PARMETIS_PREFIX} ${LDFLAGS}"
  fi
  AC_CHECK_LIB(metis,main,METIS="yes",[],-lmetis)
  AC_CHECK_LIB(parmetis,main,PARMETIS="yes",[],-lparmetis)
  if test "x${METIS}" = "xyes"; then
    if test "x${PARMETIS}" = "xyes"; then
      MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_PARMETIS_"
      AC_DEFINE(_HAVE_PARMETIS_)
      BO="${BO} Parmetis"
      if test "x${PARMETIS_PREFIX}" = "x"; then
        MAdLib_LIBS="${MAdLib_LIBS} -lparmetis -lmetis"
      else
        MAdLib_LIBS="${MAdLib_LIBS} -L${PARMETIS_PREFIX} -Wl,-rpath,${PARMETIS_PREFIX} -lparmetis -lmetis"
        MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${PARMETIS_PREFIX}"
      fi
    fi
  fi
fi

dnl Check for GSL
dnl if test "x$enable_gsl" != "xno"; then
dnl   if test "x${GSL_PREFIX}" != "x"; then
dnl     LDFLAGS="-L${GSL_PREFIX} -L${GSL_PREFIX}/lib ${LDFLAGS}"
dnl   fi
dnl   AC_CHECK_LIB(gsl,main,GSL="yes",[],-lgslcblas)
dnl   if test "x${GSL}" = "xyes"; then
dnl     MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_GSL_"
dnl     AC_DEFINE(_HAVE_GSL_)
dnl     BO="${BO} Gsl"
dnl     if test "x${GSL_PREFIX}" = "x"; then
dnl       MAdLib_LIBS="${MAdLib_LIBS} -lgsl"
dnl     else
dnl       MAdLib_LIBS="${MAdLib_LIBS} -L${GSL_PREFIX} -L${GSL_PREFIX}/lib -lgsl"
dnl       MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${GSL_PREFIX} -I${GSL_PREFIX}/include"
dnl     fi
dnl   fi
dnl fi

dnl Check for C version of BLAS
  if test "x${BLAS_LAPACK_PREFIX}" != "x"; then
    LDFLAGS="-L${BLAS_LAPACK_PREFIX} -L${BLAS_LAPACK_PREFIX}/lib ${LDFLAGS}"
  fi
  AC_CHECK_LIB(cblas,cblas_dgemm,CBLAS="yes" BLAS_LIBS="-lcblas")
  if test "x${CBLAS}" != "xyes"; then
    AC_CHECK_LIB(cblas,cblas_dgemm,CBLAS="yes" BLAS_LIBS="-lcblas -latlas",[],-latlas)
  fi
  if test "x${CBLAS}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_BLAS_"
    AC_DEFINE(_HAVE_BLAS_)
    BO="${BO} Cblas"
  else 
    if test "x${GSL}" = "xyes"; then
      dnl use unoptimized gsl version
      BLAS_LIBS="-lgslcblas"
      MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_BLAS_"
      AC_DEFINE(_HAVE_BLAS_)
      BO="${BO} Cblas"
    fi
  fi

dnl Check for Fortran version of blas and lapack (only used when not
dnl using GSL, or of FourierModel is linked in)
if test "x${FM}" = "xyes" -o "x${GSL}" != "xyes"; then
dnl   AC_PROG_F77
dnl   case "${F77}" in
dnl     *gfortran*)
dnl       F77LIB="-lgfortran"
dnl       ;;
dnl     *g77*)
dnl       F77LIB="-lg2c"
dnl       ;;
dnl     *)
dnl       F77LIB=""
dnl       ;;
dnl   esac
dnl   LDFLAGS="${LDFLAGS} ${F77LIB}"
  AC_CHECK_LIB(atlas,ATL_xerbla,
    AC_CHECK_LIB(f77blas,dgemm_,
     [BLAS="yes" BLAS_LIBS="${BLAS_LIBS} -lf77blas -latlas"],[],-latlas))
  if test "x${BLAS}" != "xyes"; then
    AC_CHECK_LIB(blas,dgemm_,[BLAS="yes" BLAS_LIBS="${BLAS_LIBS} -lblas"])
  fi
  if test "x${BLAS}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_BLAS_"
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_BLAS"
    AC_DEFINE(HAVE_BLAS)
    BO="${BO} Blas"
    AC_CHECK_LIB(lapack,dbdsqr_,
      [LAPACK="yes" BLAS_LIBS="-llapack ${BLAS_LIBS}"],[],${BLAS_LIBS})
    if test "x${LAPACK}" = "xyes"; then
      MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_LAPACK_"
      MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_LAPACK"
      AC_DEFINE(HAVE_LAPACK)
      BO="${BO} Lapack"
    fi
  fi
fi

if test "x${BLAS_LIBS}" != "x"; then
  if test "x${BLAS_LAPACK_PREFIX}" != "x"; then
    MAdLib_LIBS="${MAdLib_LIBS} -L${BLAS_LAPACK_PREFIX} -L${BLAS_LAPACK_PREFIX}/lib ${BLAS_LIBS}"
  else
    MAdLib_LIBS="${MAdLib_LIBS} ${BLAS_LIBS}"
  fi
fi

dnl Check for MPI
if test "x$enable_mpi" = "xyes"; then
  if test "x${MPI_PREFIX}" != "x"; then
    LDFLAGS="-L${MPI_PREFIX}/lib ${LDFLAGS}"
  fi
  dnl check presence of mpi.h rather than the lib, the lib is not always present 
  AC_CHECK_HEADERS(mpi.h,MPI="yes")
  dnl AC_CHECK_LIB(mpi,main,MPI="yes")
  if test "x${MPI}" = "xyes"; then
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_MPI_ -DPARALLEL"
    AC_DEFINE(_HAVE_MPI_)
    AC_DEFINE(PARALLEL)
    BO="${BO} Mpi"
    echo "********************************************************************"
    echo "Warning: MAdLib is configured with MPI enabled. It may be necessary"
    echo "to specify the values of the CXX and LINKER variables."
    echo "********************************************************************"
    if test "x${MPI_PREFIX}" = "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lmpi"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I/usr/include/mpi"
    else
      MAdLib_LIBS="${MAdLib_LIBS} -L${MPI_PREFIX}/lib -lmpi"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${MPI_PREFIX}/include"
    fi
  fi
fi

dnl Check for Autopack
if test "x$enable_mpi" = "xyes"; then
  if test "x${AUTOPACK_PREFIX}" != "x"; then
    LDFLAGS="-L${AUTOPACK_PREFIX}/lib ${LDFLAGS}"
  fi
  AC_CHECK_LIB(autopack,main,AUTOPACK="yes")
  if test "x${AUTOPACK}" = "xyes"; then
    BO="${BO} Autopack"
    if test "x${AUTOPACK_PREFIX}" = "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lautopack"
    else
      MAdLib_LIBS="${MAdLib_LIBS} -L${AUTOPACK_PREFIX}/lib -lautopack"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${AUTOPACK_PREFIX}/include"
    fi
  fi
fi

dnl Check for Parser
if test "x$enable_parser" = "xyes"; then
  if test "x${PARSER_PREFIX}" != "x"; then
    LDFLAGS="-L${PARSER_PREFIX}/lib ${LDFLAGS}"
  fi
  AC_CHECK_LIB(Parser,main,PARSER="yes")
  if test "x${PARSER}" = "xyes"; then
    BO="${BO} Parser"
    MAdLib_DEFS="${MAdLib_DEFS} -D_HAVE_PARSER_"
    AC_DEFINE(_HAVE_PARSER_)
    if test "x${PARSER_PREFIX}" = "x"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lParser"
    else
      MAdLib_LIBS="${MAdLib_LIBS} -L${PARSER_PREFIX}/lib -Wl,-rpath,${PARSER_PREFIX}/lib -lParser"
      MAdLib_INCLUDES="${MAdLib_INCLUDES} -I${PARSER_PREFIX}/include"
    fi
  fi
fi

dnl Finish link line
MAdLib_LIBS="${MAdLib_LIBS} -lm"

dnl Modify defaults according to OS
case "$UNAME" in

  CYGWIN* | MINGW*)
    dnl increase stack size to 16Mb to avoid stack overflows in
    dnl recursive tet classification for large 3D Delaunay grids
    LINKER="${LINKER} -mwindows -Wl,--stack,16777216"
    if test "x$enable_cygwin" != "xyes"; then
      MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_NO_DLL"
      AC_DEFINE(HAVE_NO_DLL)
      BO="${BO} NoDll"
    fi
    if test "x${OCC}" = "xyes"; then
      MAdLib_LIBS="${MAdLib_LIBS} -lwinspool -lws2_32"
    fi
    if test "x$enable_gui" != "xno"; then
      MAdLib_LIBS="${MAdLib_LIBS} Fltk/Win32Icon.res"
    fi
    ;;

  Darwin*)
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_NO_DLL"
    AC_DEFINE(HAVE_NO_DLL)
    BO="${BO} NoDll"
    if test "x$enable_universal" = "xyes"; then
      CXXFLAGS="${CXXFLAGS} -arch ppc -arch i386"
    fi
    if test "x$enable_gui" = "xno"; then
      MAdLib_LIBS="${MAdLib_LIBS} -framework ApplicationServices"
    fi
    ;;

  AIX*)
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_NO_DLL"
    AC_DEFINE(HAVE_NO_DLL)
    BO="${BO} NoDll"
    CXXFLAGS="${CXXFLAGS} -D_BSD"
    ;;

  IRIX*)
    dnl options for native SGI compiler
    case "${CXX}" in
      *CC*)
        dnl add "-mips3 -n32" to FLAGS, AR and LINKER for portable 32 bit exe
        dnl "-DOLDCINCLUDE" is for Netgen
        CXXFLAGS="{CXXFLAGS} -LANG:std -OPT:Olimit=0 -DOLDCINCLUDE"
        AR="${CXX} -ar -o"
        LINKER="${CXX}"
        ;;
    esac
    ;;

  OSF1*)
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_NO_SOCKLEN_T"
    AC_DEFINE(HAVE_NO_SOCKLEN_T)
    BO="${BO} NoSocklenT"
    dnl options for native DEC compiler
    case "${CXX}" in
      *cxx*)
        CXXFLAGS="${CXXFLAGS} -D__USE_STD_IOSTREAM"
        ;;
    esac
    ;;

  SunOS*)
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_NO_DLL"
    AC_DEFINE(HAVE_NO_DLL)
    BO="${BO} NoDll"
    MAdLib_LIBS="${MAdLib_LIBS} -lsocket -lnsl -ldl"
    ;;

  HP-UX*)
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_NO_DLL"
    AC_DEFINE(HAVE_NO_DLL)
    BO="${BO} NoDll"
    ;;

esac

dnl Check sizeof size_t (flag as 64 if not 32)
AC_CHECK_SIZEOF([size_t])
if test $ac_cv_sizeof_size_t != 4; then
  if test $ac_cv_sizeof_size_t != 8; then
    AC_MSG_WARN([Unsupported size of size_t - this may affect FNV hashing.])
  else
    MAdLib_DEFS="${MAdLib_DEFS} -DHAVE_64BIT_SIZE_T"
    AC_DEFINE(HAVE_64BIT_SIZE_T)
    BO="${BO} Have64BitSizeT"
    if test "x${OCC}" = "xyes"; then
      CXXFLAGS="${CXXFLAGS} -D_OCC64"
    fi
  fi
fi

dnl Check for Doxygen
DOX=""
AC_CHECK_PROG(DOX,doxygen,"yes","no")
if test "x${DOX}" = "xyes"; then
  DOXYGEN="doxygen doxygen.config"
else
  DOXYGEN=""
fi
AC_SUBST(DOXYGEN)

dnl what is the command for removing files
RM="rm -f"
AC_SUBST(RM)

AC_CONFIG_HEADER(MAdConfig.h:MAdConfig.h.in)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([Common/Makefile])
AC_CONFIG_FILES([Geo/Makefile])
AC_CONFIG_FILES([Mesh/Makefile])
AC_CONFIG_FILES([Adapt/Makefile])
AC_CONFIG_FILES([variables])
AC_CONFIG_FILES([Contrib/ANN/Makefile])
AC_CONFIG_FILES([Benchmarks/meshInfo/Makefile])
AC_CONFIG_FILES([Benchmarks/checkMesh/Makefile])
AC_CONFIG_FILES([Benchmarks/optimize/Makefile])
AC_CONFIG_FILES([Benchmarks/moveIt/Makefile])

AC_CONFIG_MACRO_DIR([m4])

includedir="${prefix}/include/MAdLib"
allincludedir="${prefix}/include_all/MAdLib"

dnl Write output
AC_SUBST(UNAME)
AC_SUBST(HOSTNAME)
AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST(LINKER)
AC_SUBST(AR)
AC_SUBST(MAdLib_DEFS)
AC_SUBST(MAdLib_DIRS)
AC_SUBST(MAdLib_LIBS)
AC_SUBST(MAdLib_INCLUDES)
AC_SUBST(MAdLib_BENCHDIRS)
AC_SUBST(MAdLib_LIBADD)
AC_SUBST(includedir)
AC_SUBST(allincludedir)

AC_OUTPUT

dnl Print some information
echo "********************************************************************"
echo "MAdLib is configured for"
echo "  - OS           : ${UNAME} on ${HOSTNAME}"
echo "  - C++ compiler : ${CXX}"
echo "  - Exe linker   : ${LINKER}"
echo "  - C++ flags    : ${CXXFLAGS}"
echo "  - Build options: ${BO}"
echo "  - Executables  : ${benchs}"
echo "********************************************************************"
echo "Edit 'variables' and 'MAdConfig.h' to fine-tune the config"
echo "********************************************************************"

