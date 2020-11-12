set(Omega_h_SOURCES
  Omega_h_adapt.cpp
  Omega_h_adj.cpp
  Omega_h_align.cpp
  Omega_h_amr.cpp
  Omega_h_amr_topology.cpp
  Omega_h_amr_transfer.cpp
  Omega_h_any.cpp
  Omega_h_approach.cpp
  Omega_h_array.cpp
  Omega_h_array_ops.cpp
  Omega_h_assoc.cpp
  Omega_h_base64.cpp
  Omega_h_bbox.cpp
  Omega_h_bcast.cpp
  Omega_h_bipart.cpp
  Omega_h_box.cpp
  Omega_h_build.cpp
  Omega_h_build_parser.cpp
  Omega_h_chartab.cpp
  Omega_h_check_collapse.cpp
  Omega_h_class.cpp
  Omega_h_cmdline.cpp
  Omega_h_coarsen.cpp
  Omega_h_coarsen_overshoot.cpp
  Omega_h_coarsen_qualities.cpp
  Omega_h_coarsen_topology.cpp
  Omega_h_collapse_rail.cpp
  Omega_h_comm.cpp
  Omega_h_compare.cpp
  Omega_h_confined.cpp
  Omega_h_conserve.cpp
  Omega_h_dist.cpp
  Omega_h_eigen.cpp
  Omega_h_expr.cpp
  Omega_h_fail.cpp
  Omega_h_fence.cpp
  Omega_h_file.cpp
  Omega_h_filesystem.cpp
  Omega_h_finite_automaton.cpp
  Omega_h_fix.cpp
  Omega_h_future.cpp
  Omega_h_ghost.cpp
  Omega_h_globals.cpp
  Omega_h_gmsh.cpp
  Omega_h_grammar.cpp
  Omega_h_graph.cpp
  Omega_h_hilbert.cpp
  Omega_h_histogram.cpp
  Omega_h_indset.cpp
  Omega_h_inertia.cpp
  Omega_h_input.cpp
  Omega_h_int128.cpp
  Omega_h_int_scan.cpp
  Omega_h_language.cpp
  Omega_h_laplace.cpp
  Omega_h_library.cpp
  Omega_h_linpart.cpp
  Omega_h_malloc.cpp
  Omega_h_map.cpp
  Omega_h_mark.cpp
  Omega_h_math_lang.cpp
  Omega_h_matrix.cpp
  Omega_h_mesh.cpp
  Omega_h_metric.cpp
  Omega_h_metric_input.cpp
  Omega_h_migrate.cpp
  Omega_h_modify.cpp
  Omega_h_owners.cpp
  Omega_h_parser.cpp
  Omega_h_parser_graph.cpp
  Omega_h_pool.cpp
  Omega_h_print.cpp
  Omega_h_profile.cpp
  Omega_h_quality.cpp
  Omega_h_reader.cpp
  Omega_h_recover.cpp
  Omega_h_refine.cpp
  Omega_h_refine_qualities.cpp
  Omega_h_refine_topology.cpp
  Omega_h_regex.cpp
  Omega_h_regress.cpp
  Omega_h_remotes.cpp
  Omega_h_reorder.cpp
  Omega_h_scatterplot.cpp
  Omega_h_shape.cpp
  Omega_h_shared_alloc.cpp
  Omega_h_simplify.cpp
  Omega_h_sort.cpp
  Omega_h_surface.cpp
  Omega_h_swap.cpp
  Omega_h_swap2d.cpp
  Omega_h_swap2d_qualities.cpp
  Omega_h_swap2d_topology.cpp
  Omega_h_swap3d.cpp
  Omega_h_swap3d_qualities.cpp
  Omega_h_swap3d_topology.cpp
  Omega_h_tag.cpp
  Omega_h_timer.cpp
  Omega_h_transfer.cpp
  Omega_h_unmap_mesh.cpp
  Omega_h_vector.cpp
  Omega_h_vtk.cpp
  Omega_h_xml.cpp
  Omega_h_xml_lite.cpp
  Omega_h_yaml.cpp
  )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
else()
  set(Omega_h_SOURCES ${Omega_h_SOURCES}
    Omega_h_random.cpp
    Omega_h_coarsen_flip.cpp
  )
endif()

if (Omega_h_USE_CUDA)
  if (NOT Omega_h_USE_Kokkos)
    set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_for.cpp)
  endif()
else()
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_overlay.cpp)
endif()

if(Omega_h_USE_libMeshb)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_meshb.cpp)
endif()

if(Omega_h_USE_EGADS)
  message(STATUS "Omega_h_USE_EGADS: ${Omega_h_USE_EGADS}")
  if(NOT EGADS_PREFIX)
    message(FATAL_ERROR "Omega_h_USE_EGADS=${Omega_h_USE_EGADS} but EGADS_PREFIX not set!")
  endif()
  find_path(EGADS_INCLUDE_DIR NAMES egads.h
      PATHS
      "${EGADS_PREFIX}/EngSketchPad/include"
      "${EGADS_PREFIX}/include"
      NO_DEFAULT_PATH)
  if(NOT EGADS_INCLUDE_DIR)
    message(FATAL_ERROR "could not find path to \"egads.h\"")
  endif()
  find_library(EGADS_LIBRARY NAMES egads
      PATHS
      "${EGADS_PREFIX}/EngSketchPad/lib"
      NO_DEFAULT_PATH)
  if(NOT EGADS_LIBRARY)
    message(FATAL_ERROR "could not find EGADS_LIBRARY")
  endif()
  message(STATUS "EGADS_INCLUDE_DIR: ${EGADS_INCLUDE_DIR}")
  message(STATUS "EGADS_LIBRARY: ${EGADS_LIBRARY}")
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_egads.cpp)
endif()

if(Omega_h_USE_SEACASExodus)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_exodus.cpp)
endif()

if(Omega_h_USE_DOLFIN)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_dolfin.cpp)
endif()

if (Omega_h_USE_CUDA)
  set_source_files_properties(${Omega_h_SOURCES} PROPERTIES LANGUAGE CUDA)
endif()

add_library(omega_h ${Omega_h_SOURCES})

set_property(TARGET omega_h PROPERTY CXX_STANDARD "11")
set_property(TARGET omega_h PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET omega_h PROPERTY CXX_EXTENSIONS OFF)
set_property(TARGET omega_h PROPERTY VERSION ${Omega_h_VERSION})
set_property(TARGET omega_h PROPERTY SOVERSION ${Omega_h_VERSION_MAJOR})

bob_library_includes(omega_h)

get_filename_component(TPL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../tpl" ABSOLUTE)
target_include_directories(omega_h PRIVATE $<BUILD_INTERFACE:${TPL_PATH}>)

if (Omega_h_USE_OpenMP)
  target_compile_options(omega_h PUBLIC -fopenmp)
endif()

if (Omega_h_USE_CUDA)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # -DNDEBUG causes an "unexpected number of captures" error when compiling
    # Omega_h_coarsen_flip.cpp. Perhaps the error is related to why Omega_h
    # does not let user set CMAKE_BUILD_TYPE if using bob_*_cxx_flags.
    target_compile_options(omega_h PRIVATE
        --compiler-options -W,-Wall,-Wextra,-Werror,-Wno-noexcept-type
        --Werror cross-execution-space-call,deprecated-declarations
        --expt-extended-lambda
        -UNDEBUG)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(omega_h PRIVATE
        --compiler-options /wd6294,/wd26451,/wd26439
        --Werror cross-execution-space-call,deprecated-declarations
        --expt-extended-lambda)
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(omega_h PRIVATE
        -W
        -Wall
        -Wextra
        -Werror
        -Wimplicit-fallthrough
        )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(omega_h PRIVATE
        -W
        -Wall
        -Wextra
        -Werror
        )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(omega_h PRIVATE
        "/wd26451" #arithmetic overflow, long(int * int), shows up in much vector access
        "/wd26439" #noexcept advice. We've determined that noexcept is useless for performance, we don't use it
        )
  endif()
endif()

bob_link_dependency(omega_h PUBLIC Kokkos)

bob_link_dependency(omega_h PUBLIC libMeshb)

if(Omega_h_USE_EGADS)
  target_include_directories(omega_h PUBLIC "${EGADS_INCLUDE_DIR}")
  target_link_libraries(omega_h PUBLIC "${EGADS_LIBRARY}")
endif()

bob_link_dependency(omega_h PUBLIC SEACASExodus)

bob_link_dependency(omega_h PUBLIC ZLIB)

if (Omega_h_USE_dwarf)
  target_include_directories(omega_h PRIVATE "${LIBDWARF_INCLUDE_DIRS}")
  target_link_libraries(omega_h PUBLIC "${LIBDWARF_LIBRARIES}")
endif()

if(Omega_h_USE_DOLFIN)
  target_link_libraries(omega_h PUBLIC Boost::timer dolfin)
  target_include_directories(omega_h SYSTEM PUBLIC ${DOLFIN_INCLUDE_DIRS})
  target_include_directories(omega_h SYSTEM PUBLIC ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set_target_properties(omega_h PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

bob_export_target(omega_h)

bob_config_header("${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h")
set(Omega_h_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h"
  Omega_h_adapt.hpp
  Omega_h_adj.hpp
  Omega_h_affine.hpp
  Omega_h_align.hpp
  Omega_h_amr.hpp
  Omega_h_any.hpp
  Omega_h_array.hpp
  Omega_h_array_ops.hpp
  Omega_h_assoc.hpp
  Omega_h_base64.hpp
  Omega_h_bbox.hpp
  Omega_h_box.hpp
  Omega_h_build.hpp
  Omega_h_class.hpp
  Omega_h_cmdline.hpp
  Omega_h_comm.hpp
  Omega_h_compare.hpp
  Omega_h_defines.hpp
  Omega_h_dist.hpp
  Omega_h_eigen.hpp
  Omega_h_element.hpp
  Omega_h_expr.hpp
  Omega_h_fail.hpp
  Omega_h_fence.hpp
  Omega_h_few.hpp
  Omega_h_file.hpp
  Omega_h_filesystem.hpp
  Omega_h_finite_automaton.hpp
  Omega_h_for.hpp
  Omega_h_functors.hpp
  Omega_h_future.hpp
  Omega_h_grammar.hpp
  Omega_h_graph.hpp
  Omega_h_hypercube.hpp
  Omega_h_indset.hpp
  Omega_h_indset_inline.hpp
  Omega_h_input.hpp
  Omega_h_int128.hpp
  Omega_h_int_iterator.hpp
  Omega_h_int_scan.hpp
  Omega_h_kokkos.hpp
  Omega_h_language.hpp
  Omega_h_library.hpp
  Omega_h_lie.hpp
  Omega_h_macros.h
  Omega_h_malloc.hpp
  Omega_h_map.hpp
  Omega_h_mark.hpp
  Omega_h_math_lang.hpp
  Omega_h_matrix.hpp
  Omega_h_mesh.hpp
  Omega_h_metric.hpp
  Omega_h_mpi.h
  Omega_h_owners.hpp
  Omega_h_parser.hpp
  Omega_h_print.hpp
  Omega_h_profile.hpp
  Omega_h_qr.hpp
  Omega_h_quality.hpp
  Omega_h_random.hpp
  Omega_h_rbtree.hpp
  Omega_h_reader.hpp
  Omega_h_reader_tables.hpp
  Omega_h_recover.hpp
  Omega_h_reduce.hpp
  Omega_h_remotes.hpp
  Omega_h_scalar.hpp
  Omega_h_scan.hpp
  Omega_h_scatterplot.hpp
  Omega_h_shape.hpp
  Omega_h_shared_alloc.hpp
  Omega_h_simplex.hpp
  Omega_h_sort.hpp
  Omega_h_std_vector.hpp
  Omega_h_svd.hpp
  Omega_h_table.hpp
  Omega_h_tag.hpp
  Omega_h_template_up.hpp
  Omega_h_timer.hpp
  Omega_h_vector.hpp
  Omega_h_vtk.hpp
  Omega_h_xml_lite.hpp
  r3d.hpp
  )

if (NOT Omega_h_USE_CUDA)
  set(Omega_h_HEADERS ${Omega_h_HEADERS} Omega_h_overlay.hpp)
endif()

if(Omega_h_USE_DOLFIN)
  set(Omega_h_HEADERS ${Omega_h_HEADERS} Omega_h_dolfin.hpp)
endif()

install(FILES ${Omega_h_HEADERS} DESTINATION include/omega_h)

if (Omega_h_USE_pybind11)
  set(PYBIND11_SOURCES
      PyOmega_h.cpp
      PyOmega_h_defines.cpp
      PyOmega_h_array.cpp
      PyOmega_h_comm.cpp
      PyOmega_h_library.cpp
      PyOmega_h_mesh.cpp
      PyOmega_h_build.cpp
      PyOmega_h_adapt.cpp
      PyOmega_h_file.cpp
      PyOmega_h_class.cpp
      )
  if (Omega_h_USE_DOLFIN)
    set(PYBIND11_SOURCES ${PYBIND11_SOURCES}
        PyOmega_h_dolfin.cpp
       )
  endif()
  pybind11_add_module(PyOmega_h
      ${PYBIND11_SOURCES})
  target_link_libraries(PyOmega_h PUBLIC omega_h)
  if (PYTHON_VERSION_MAJOR EQUAL 3)
    set(PyOmega_h_DEST "lib/python3/dist-packages")
  else()
    set(PyOmega_h_DEST "lib/python/dist-packages")
  endif()
  install(TARGETS PyOmega_h
      ARCHIVE DESTINATION "${PyOmega_h_DEST}"
      LIBRARY DESTINATION "${PyOmega_h_DEST}")
endif()

bob_end_subdir()
