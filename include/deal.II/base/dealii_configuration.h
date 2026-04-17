// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_base_dealii_configuration_h
#define dealii_base_dealii_configuration_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace that contains `constexpr` variables and functions that
 * describe the system as determined by CMake. Each of the variables and
 * functions in this namespace correspond to preprocessor variables and
 * macros also described in the file
 * `<deal.II/base/config.h>`.
 */
namespace configuration
{
  /**
   * A variable that equals the preprocessor value of the DEAL_II_PACKAGE_NAME
   * preprocessor symbol as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
  constexpr const char package_name[] = DEAL_II_PACKAGE_NAME;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_PACKAGE_VERSION preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file.
   */
  constexpr const char package_version[] = DEAL_II_PACKAGE_VERSION;

  /**
   * A variable that equals the preprocessor value of the DEAL_II_VERSION_MAJOR
   * preprocessor symbol as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
  constexpr const unsigned int dealii_version_major = DEAL_II_VERSION_MAJOR;

  /**
   * A variable that equals the preprocessor value of the DEAL_II_VERSION_MINOR
   * preprocessor symbol as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
  constexpr const unsigned int dealii_version_minor = DEAL_II_VERSION_MINOR;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_VERSION_SUBMINOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file.
   */
  constexpr const unsigned int dealii_version_subminor =
    DEAL_II_VERSION_SUBMINOR;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_VECTORIZATION_WIDTH_IN_BITS preprocessor symbol as determined by
   * CMake and stored in the include/deal.II/base/config.h file.
   */
  constexpr const unsigned int vectorization_width_in_bits =
    DEAL_II_VECTORIZATION_WIDTH_IN_BITS;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_COMPILER_VECTORIZATION_LEVEL preprocessor symbol as determined by
   * CMake and stored in the include/deal.II/base/config.h file.
   */
  constexpr const unsigned int compiler_vectorization_level =
    DEAL_II_COMPILER_VECTORIZATION_LEVEL;

  // Configured deal.II features
  /**
   * A variable that equals whether the DEAL_II_WITH_64BIT_INDICES preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
  constexpr const bool with_64bit_indices = true;
#else
  constexpr const bool         with_64bit_indices                       = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_ADOLC preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_ADOLC
  constexpr const bool with_adolc = true;
#else
  constexpr const bool         with_adolc                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_ARPACK preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_ARPACK
  constexpr const bool with_arpack = true;
#else
  constexpr const bool         with_arpack                              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_ARBORX preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_ARBORX
  constexpr const bool with_arborx = true;
#else
  constexpr const bool         with_arborx                              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_ASSIMP preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_ASSIMP
  constexpr const bool with_assimp = true;
#else
  constexpr const bool         with_assimp                              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED
  constexpr const bool feature_boost_bundled_configured = true;
#else
  constexpr const bool         feature_boost_bundled_configured         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_CGAL preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_CGAL
  constexpr const bool with_cgal = true;
#else
  constexpr const bool         with_cgal                                = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_COMPLEX_VALUES preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_COMPLEX_VALUES
  constexpr const bool with_complex_values = true;
#else
  constexpr const bool         with_complex_values                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_GINKGO preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_GINKGO
  constexpr const bool with_ginkgo = true;
#else
  constexpr const bool         with_ginkgo                              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_GMSH preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_GMSH
  constexpr const bool with_gmsh = true;
#else
  constexpr const bool         with_gmsh                                = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_GSL preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_GSL
  constexpr const bool with_gsl = true;
#else
  constexpr const bool         with_gsl                                 = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_HDF5 preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_HDF5
  constexpr const bool with_hdf5 = true;
#else
  constexpr const bool         with_hdf5                                = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_KOKKOS preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_KOKKOS
  constexpr const bool with_kokkos = true;
#else
  constexpr const bool         with_kokkos                              = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_FEATURE_KOKKOS_BUNDLED_CONFIGURED preprocessor symbol is defined as
   * determined by CMake and stored in the include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_FEATURE_KOKKOS_BUNDLED_CONFIGURED
  constexpr const bool feature_kokkos_bundled_configured = true;
#else
  constexpr const bool         feature_kokkos_bundled_configured        = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_LAPACK preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_LAPACK
  constexpr const bool with_lapack = true;
#else
  constexpr const bool         with_lapack                              = false;
#endif

  /**
   * A variable that equals whether the LAPACK_WITH_64BIT_BLAS_INDICES
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef LAPACK_WITH_64BIT_BLAS_INDICES
  constexpr const bool lapack_with_64bit_blas_indices = true;
#else
  constexpr const bool         lapack_with_64bit_blas_indices           = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_LAPACK_WITH_MKL preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_LAPACK_WITH_MKL
  constexpr const bool lapack_with_mkl = true;
#else
  constexpr const bool         lapack_with_mkl                          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_MAGIC_ENUM preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_MAGIC_ENUM
  constexpr const bool with_magic_enum = true;
#else
  constexpr const bool         with_magic_enum                          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_METIS preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_METIS
  constexpr const bool with_metis = true;
#else
  constexpr const bool         with_metis                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_MPI preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_MPI
  constexpr const bool with_mpi = true;
#else
  constexpr const bool         with_mpi                                 = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_MUPARSER preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_MUPARSER
  constexpr const bool with_muparser = true;
#else
  constexpr const bool         with_muparser                            = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_MUMPS preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_MUMPS
  constexpr const bool with_mumps = true;
#else
  constexpr const bool         with_mumps                               = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_FEATURE_MUPARSER_BUNDLED_CONFIGURED preprocessor symbol is defined
   * as determined by CMake and stored in the include/deal.II/base/config.h
   * file.
   */
#ifdef DEAL_II_FEATURE_MUPARSER_BUNDLED_CONFIGURED
  constexpr const bool feature_muparser_bundled_configured = true;
#else
  constexpr const bool         feature_muparser_bundled_configured      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_OPENCASCADE preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_OPENCASCADE
  constexpr const bool with_opencascade = true;
#else
  constexpr const bool         with_opencascade                         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_P4EST preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_P4EST
  constexpr const bool with_p4est = true;
#else
  constexpr const bool         with_p4est                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_PETSC preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_PETSC
  constexpr const bool with_petsc = true;
#else
  constexpr const bool         with_petsc                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_PSBLAS preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_PSBLAS
  constexpr const bool with_psblas = true;
#else
  constexpr const bool         with_psblas                              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_SCALAPACK preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_SCALAPACK
  constexpr const bool with_scalapack = true;
#else
  constexpr const bool         with_scalapack                           = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_SLEPC preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_SLEPC
  constexpr const bool with_slepc = true;
#else
  constexpr const bool         with_slepc                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_SUNDIALS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_SUNDIALS
  constexpr const bool with_sundials = true;
#else
  constexpr const bool         with_sundials                            = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_SYMENGINE preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_SYMENGINE
  constexpr const bool with_symengine = true;
#else
  constexpr const bool         with_symengine                           = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_TASKFLOW preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_TASKFLOW
  constexpr const bool with_taskflow = true;
#else
  constexpr const bool         with_taskflow                            = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_FEATURE_TASKFLOW_BUNDLED_CONFIGURED preprocessor symbol is defined
   * as determined by CMake and stored in the include/deal.II/base/config.h
   * file.
   */
#ifdef DEAL_II_FEATURE_TASKFLOW_BUNDLED_CONFIGURED
  constexpr const bool feature_taskflow_bundled_configured = true;
#else
  constexpr const bool         feature_taskflow_bundled_configured      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_TBB preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_TBB
  constexpr const bool with_tbb = true;
#else
  constexpr const bool         with_tbb                                 = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_FEATURE_TBB_BUNDLED_CONFIGURED
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_FEATURE_TBB_BUNDLED_CONFIGURED
  constexpr const bool feature_tbb_bundled_configured = true;
#else
  constexpr const bool         feature_tbb_bundled_configured           = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_TRILINOS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_TRILINOS
  constexpr const bool with_trilinos = true;
#else
  constexpr const bool         with_trilinos                            = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_UMFPACK preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_UMFPACK
  constexpr const bool with_umfpack = true;
#else
  constexpr const bool         with_umfpack                             = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_FEATURE_UMFPACK_BUNDLED_CONFIGURED preprocessor symbol is defined
   * as determined by CMake and stored in the include/deal.II/base/config.h
   * file.
   */
#ifdef DEAL_II_FEATURE_UMFPACK_BUNDLED_CONFIGURED
  constexpr const bool feature_umfpack_bundled_configured = true;
#else
  constexpr const bool         feature_umfpack_bundled_configured       = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_USE_VECTORIZATION_GATHER
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_USE_VECTORIZATION_GATHER
  constexpr const bool use_vectorization_gather = true;
#else
  constexpr const bool         use_vectorization_gather                 = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_VTK preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_VTK
  constexpr const bool with_vtk = true;
#else
  constexpr const bool         with_vtk                                 = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_ZLIB preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_ZLIB
  constexpr const bool with_zlib = true;
#else
  constexpr const bool         with_zlib                                = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_THREADS preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_THREADS
  constexpr const bool with_threads = true;
#else
  constexpr const bool         with_threads                             = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TBB_WITH_ONEAPI preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TBB_WITH_ONEAPI
  constexpr const bool tbb_with_oneapi = true;
#else
  constexpr const bool         tbb_with_oneapi                          = false;
#endif

  // Compiler bugs
  /**
   * A variable that equals whether the DEAL_II_DELETED_MOVE_CONSTRUCTOR_BUG
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_DELETED_MOVE_CONSTRUCTOR_BUG
  constexpr const bool deleted_move_constructor_bug = true;
#else
  constexpr const bool         deleted_move_constructor_bug             = false;
#endif

  // Compiler features
  /**
   * A variable that equals whether the DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
  constexpr const bool compiler_use_vector_arithmetics = true;
#else
  constexpr const bool         compiler_use_vector_arithmetics          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_VECTOR_ITERATOR_IS_POINTER
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_VECTOR_ITERATOR_IS_POINTER
  constexpr const bool vector_iterator_is_pointer = true;
#else
  constexpr const bool         vector_iterator_is_pointer               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_BUILTIN_EXPECT preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_BUILTIN_EXPECT
  constexpr const bool have_builtin_expect = true;
#else
  constexpr const bool         have_builtin_expect                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_GLIBC_STACKTRACE
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  constexpr const bool have_glibc_stacktrace = true;
#else
  constexpr const bool         have_glibc_stacktrace                    = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_LIBSTDCXX_DEMANGLER
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_LIBSTDCXX_DEMANGLER
  constexpr const bool have_libstdcxx_demangler = true;
#else
  constexpr const bool         have_libstdcxx_demangler                 = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_COMPILER_HAS_DIAGNOSTIC_PRAGMA
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_COMPILER_HAS_DIAGNOSTIC_PRAGMA
  constexpr const bool compiler_has_diagnostic_pragma = true;
#else
  constexpr const bool         compiler_has_diagnostic_pragma           = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_KOKKOS_ENABLE_HIP preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_KOKKOS_ENABLE_HIP
  constexpr const bool kokkos_enable_hip = true;
#else
  constexpr const bool         kokkos_enable_hip                        = false;
#endif

  // CPU features
  /**
   * A variable that equals whether the DEAL_II_WORDS_BIGENDIAN preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WORDS_BIGENDIAN
  constexpr const bool words_bigendian = true;
#else
  constexpr const bool         words_bigendian                          = false;
#endif

  // Language features
  /**
   * A variable that equals whether the DEAL_II_HAVE_CXX14 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_CXX14
  constexpr const bool have_cxx14 = true;
#else
  constexpr const bool         have_cxx14                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_CXX17 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_CXX17
  constexpr const bool have_cxx17 = true;
#else
  constexpr const bool         have_cxx17                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_CXX20 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_CXX20
  constexpr const bool have_cxx20 = true;
#else
  constexpr const bool         have_cxx20                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_CXX23 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_CXX23
  constexpr const bool have_cxx23 = true;
#else
  constexpr const bool         have_cxx23                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_FP_EXCEPTIONS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  constexpr const bool have_fp_exceptions = true;
#else
  constexpr const bool         have_fp_exceptions                       = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
  constexpr const bool have_complex_operator_overloads = true;
#else
  constexpr const bool         have_complex_operator_overloads          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
  constexpr const bool have_cxx17_bessel_functions = true;
#else
  constexpr const bool         have_cxx17_bessel_functions              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_CXX14_CONSTEXPR_BUG preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_CXX14_CONSTEXPR_BUG
  constexpr const bool cxx14_constexpr_bug = true;
#else
  constexpr const bool         cxx14_constexpr_bug                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_CXX11 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_CXX11
  constexpr const bool with_cxx11 = true;
#else
  constexpr const bool         with_cxx11                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_CXX14 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_CXX14
  constexpr const bool with_cxx14 = true;
#else
  constexpr const bool         with_cxx14                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_WITH_CXX17 preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_WITH_CXX17
  constexpr const bool with_cxx17 = true;
#else
  constexpr const bool         with_cxx17                               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_EARLY_DEPRECATIONS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_EARLY_DEPRECATIONS
  constexpr const bool early_deprecations = true;
#else
  constexpr const bool         early_deprecations                       = false;
#endif

  // System features
  /**
   * A variable that equals whether the DEAL_II_HAVE_SYS_RESOURCE_H preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_SYS_RESOURCE_H
  constexpr const bool have_sys_resource_h = true;
#else
  constexpr const bool         have_sys_resource_h                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_UNISTD_H preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_UNISTD_H
  constexpr const bool have_unistd_h = true;
#else
  constexpr const bool         have_unistd_h                            = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_GETHOSTNAME preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_GETHOSTNAME
  constexpr const bool have_gethostname = true;
#else
  constexpr const bool         have_gethostname                         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_GETPID preprocessor symbol
   * is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_GETPID
  constexpr const bool have_getpid = true;
#else
  constexpr const bool         have_getpid                              = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_HAVE_JN preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_HAVE_JN
  constexpr const bool have_jn = true;
#else
  constexpr const bool         have_jn                                  = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_MSVC preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_MSVC
  constexpr const bool msvc = true;
#else
  constexpr const bool         msvc                                     = false;
#endif

  // Feature configuration
  /**
   * A variable that equals whether the DEAL_II_ADOLC_WITH_ATRIG_ERF
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_ADOLC_WITH_ATRIG_ERF
  constexpr const bool adolc_with_atrig_erf = true;
#else
  constexpr const bool         adolc_with_atrig_erf                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING
  constexpr const bool adolc_with_advanced_branching = true;
#else
  constexpr const bool         adolc_with_advanced_branching            = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_ADOLC_WITH_TAPELESS_REFCOUNTING
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_ADOLC_WITH_TAPELESS_REFCOUNTING
  constexpr const bool adolc_with_tapeless_refcounting = true;
#else
  constexpr const bool         adolc_with_tapeless_refcounting          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_ARBORX_WITH_MPI preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_ARBORX_WITH_MPI
  constexpr const bool arborx_with_mpi = true;
#else
  constexpr const bool         arborx_with_mpi                          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_ARPACK_WITH_PARPACK preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_ARPACK_WITH_PARPACK
  constexpr const bool arpack_with_parpack = true;
#else
  constexpr const bool         arpack_with_parpack                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_GMSH_WITH_API preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_GMSH_WITH_API
  constexpr const bool gmsh_with_api = true;
#else
  constexpr const bool         gmsh_with_api                            = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_PETSC_WITH_COMPLEX preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_PETSC_WITH_COMPLEX
  constexpr const bool petsc_with_complex = true;
#else
  constexpr const bool         petsc_with_complex                       = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_PETSC_WITH_HYPRE preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_PETSC_WITH_HYPRE
  constexpr const bool petsc_with_hypre = true;
#else
  constexpr const bool         petsc_with_hypre                         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_PETSC_WITH_MUMPS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_PETSC_WITH_MUMPS
  constexpr const bool petsc_with_mumps = true;
#else
  constexpr const bool         petsc_with_mumps                         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_PETSC_WITH_KOKKOS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_PETSC_WITH_KOKKOS
  constexpr const bool petsc_with_kokkos = true;
#else
  constexpr const bool         petsc_with_kokkos                        = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_SUNDIALS_WITH_IDAS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_SUNDIALS_WITH_IDAS
  constexpr const bool sundials_with_idas = true;
#else
  constexpr const bool         sundials_with_idas                       = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_SYMENGINE_WITH_LLVM preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_SYMENGINE_WITH_LLVM
  constexpr const bool symengine_with_llvm = true;
#else
  constexpr const bool         symengine_with_llvm                      = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS preprocessor symbol is defined
   * as determined by CMake and stored in the include/deal.II/base/config.h
   * file.
   */
#ifdef DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS
  constexpr const bool boost_has_broken_header_deprecations = true;
#else
  constexpr const bool         boost_has_broken_header_deprecations     = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_CGAL_HAS_DEPRECATED_BOOST_INCLUDES preprocessor symbol is defined
   * as determined by CMake and stored in the include/deal.II/base/config.h
   * file.
   */
#ifdef DEAL_II_CGAL_HAS_DEPRECATED_BOOST_INCLUDES
  constexpr const bool cgal_has_deprecated_boost_includes = true;
#else
  constexpr const bool         cgal_has_deprecated_boost_includes       = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
  constexpr const bool trilinos_cxx_supports_sacado_complex_rad = true;
#else
  constexpr const bool         trilinos_cxx_supports_sacado_complex_rad = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_AMESOS
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_AMESOS
  constexpr const bool trilinos_with_amesos = true;
#else
  constexpr const bool         trilinos_with_amesos                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_AMESOS2
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_AMESOS2
  constexpr const bool trilinos_with_amesos2 = true;
#else
  constexpr const bool         trilinos_with_amesos2                    = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_AZTECOO
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_AZTECOO
  constexpr const bool trilinos_with_aztecoo = true;
#else
  constexpr const bool         trilinos_with_aztecoo                    = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_BELOS preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_BELOS
  constexpr const bool trilinos_with_belos = true;
#else
  constexpr const bool         trilinos_with_belos                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_EPETRA
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_EPETRA
  constexpr const bool trilinos_with_epetra = true;
#else
  constexpr const bool         trilinos_with_epetra                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_EPETRAEXT
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_EPETRAEXT
  constexpr const bool trilinos_with_epetraext = true;
#else
  constexpr const bool         trilinos_with_epetraext                  = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_IFPACK
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_IFPACK
  constexpr const bool trilinos_with_ifpack = true;
#else
  constexpr const bool         trilinos_with_ifpack                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_IFPACK2
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_IFPACK2
  constexpr const bool trilinos_with_ifpack2 = true;
#else
  constexpr const bool         trilinos_with_ifpack2                    = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_ML preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_ML
  constexpr const bool trilinos_with_ml = true;
#else
  constexpr const bool         trilinos_with_ml                         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_MUELU preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_MUELU
  constexpr const bool trilinos_with_muelu = true;
#else
  constexpr const bool         trilinos_with_muelu                      = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_NOX preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_NOX
  constexpr const bool trilinos_with_nox = true;
#else
  constexpr const bool         trilinos_with_nox                        = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_ROL preprocessor
   * symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_ROL
  constexpr const bool trilinos_with_rol = true;
#else
  constexpr const bool         trilinos_with_rol                        = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_SACADO
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_SACADO
  constexpr const bool trilinos_with_sacado = true;
#else
  constexpr const bool         trilinos_with_sacado                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_SEACAS
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_SEACAS
  constexpr const bool trilinos_with_seacas = true;
#else
  constexpr const bool         trilinos_with_seacas                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_TEUCHOS
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TEUCHOS
  constexpr const bool trilinos_with_teuchos = true;
#else
  constexpr const bool         trilinos_with_teuchos                    = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_TPETRA
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TPETRA
  constexpr const bool trilinos_with_tpetra = true;
#else
  constexpr const bool         trilinos_with_tpetra                     = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_TRILINOS_WITH_TPETRA_INST_COMPLEX_DOUBLE preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_COMPLEX_DOUBLE
  constexpr const bool trilinos_with_tpetra_inst_complex_double = true;
#else
  constexpr const bool         trilinos_with_tpetra_inst_complex_double = false;
#endif

  /**
   * A variable that equals whether the
   * DEAL_II_TRILINOS_WITH_TPETRA_INST_COMPLEX_FLOAT preprocessor symbol is
   * defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_COMPLEX_FLOAT
  constexpr const bool trilinos_with_tpetra_inst_complex_float = true;
#else
  constexpr const bool         trilinos_with_tpetra_inst_complex_float  = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_TPETRA_INST_DOUBLE
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_DOUBLE
  constexpr const bool trilinos_with_tpetra_inst_double = true;
#else
  constexpr const bool         trilinos_with_tpetra_inst_double         = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_TPETRA_INST_FLOAT
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_FLOAT
  constexpr const bool trilinos_with_tpetra_inst_float = true;
#else
  constexpr const bool         trilinos_with_tpetra_inst_float          = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_TPETRA_MUELU
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU
  constexpr const bool trilinos_with_tpetra_muelu = true;
#else
  constexpr const bool         trilinos_with_tpetra_muelu               = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_TRILINOS_WITH_ZOLTAN
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
  constexpr const bool trilinos_with_zoltan = true;
#else
  constexpr const bool         trilinos_with_zoltan                     = false;
#endif

  /**
   * A variable that equals whether the DEAL_II_MPI_WITH_DEVICE_SUPPORT
   * preprocessor symbol is defined as determined by CMake and stored in the
   * include/deal.II/base/config.h file.
   */
#ifdef DEAL_II_MPI_WITH_DEVICE_SUPPORT
  constexpr const bool mpi_with_device_support = true;
#else
  constexpr const bool         mpi_with_device_support                  = false;
#endif

  // Version information
  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_BOOST_VERSION_MAJOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file.
   */
  constexpr const unsigned int boost_version_major =
    DEAL_II_BOOST_VERSION_MAJOR;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_BOOST_VERSION_MINOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file.
   */
  constexpr const unsigned int boost_version_minor =
    DEAL_II_BOOST_VERSION_MINOR;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_BOOST_VERSION_SUBMINOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file.
   */
  constexpr const unsigned int boost_version_subminor =
    DEAL_II_BOOST_VERSION_SUBMINOR;

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_KOKKOS_VERSION_MAJOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_KOKKOS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_KOKKOS
  constexpr const unsigned int kokkos_version_major =
    DEAL_II_KOKKOS_VERSION_MAJOR;
#else
  constexpr const unsigned int kokkos_version_major                     = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_KOKKOS_VERSION_MINOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_KOKKOS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_KOKKOS
  constexpr const unsigned int kokkos_version_minor =
    DEAL_II_KOKKOS_VERSION_MINOR;
#else
  constexpr const unsigned int kokkos_version_minor                     = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_KOKKOS_VERSION_SUBMINOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_KOKKOS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_KOKKOS
  constexpr const unsigned int kokkos_version_subminor =
    DEAL_II_KOKKOS_VERSION_SUBMINOR;
#else
  constexpr const unsigned int kokkos_version_subminor                  = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_OPENCASCADE_VERSION_MAJOR preprocessor symbol as determined by
   * CMake and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_OPENCASCADE preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_OPENCASCADE
  constexpr const unsigned int opencascade_version_major =
    DEAL_II_OPENCASCADE_VERSION_MAJOR;
#else
  constexpr const unsigned int opencascade_version_major                = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_OPENCASCADE_VERSION_MINOR preprocessor symbol as determined by
   * CMake and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_OPENCASCADE preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_OPENCASCADE
  constexpr const unsigned int opencascade_version_minor =
    DEAL_II_OPENCASCADE_VERSION_MINOR;
#else
  constexpr const unsigned int opencascade_version_minor                = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_OPENCASCADE_VERSION_SUBMINOR preprocessor symbol as determined by
   * CMake and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_OPENCASCADE preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_OPENCASCADE
  constexpr const unsigned int opencascade_version_subminor =
    DEAL_II_OPENCASCADE_VERSION_SUBMINOR;
#else
  constexpr const unsigned int opencascade_version_subminor             = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_P4EST_VERSION_MAJOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_P4EST preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_P4EST
  constexpr const unsigned int p4est_version_major =
    DEAL_II_P4EST_VERSION_MAJOR;
#else
  constexpr const unsigned int p4est_version_major                      = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_P4EST_VERSION_MINOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_P4EST preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_P4EST
  constexpr const unsigned int p4est_version_minor =
    DEAL_II_P4EST_VERSION_MINOR;
#else
  constexpr const unsigned int p4est_version_minor                      = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_P4EST_VERSION_SUBMINOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_P4EST preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_P4EST
  constexpr const unsigned int p4est_version_subminor =
    DEAL_II_P4EST_VERSION_SUBMINOR;
#else
  constexpr const unsigned int p4est_version_subminor                   = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_P4EST_VERSION_PATCH preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_P4EST preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_P4EST
  constexpr const unsigned int p4est_version_patch =
    DEAL_II_P4EST_VERSION_PATCH;
#else
  constexpr const unsigned int p4est_version_patch                      = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_SUNDIALS_VERSION_MAJOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_SUNDIALS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_SUNDIALS
  constexpr const unsigned int sundials_version_major =
    DEAL_II_SUNDIALS_VERSION_MAJOR;
#else
  constexpr const unsigned int sundials_version_major                   = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_SUNDIALS_VERSION_MINOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_SUNDIALS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_SUNDIALS
  constexpr const unsigned int sundials_version_minor =
    DEAL_II_SUNDIALS_VERSION_MINOR;
#else
  constexpr const unsigned int sundials_version_minor                   = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_SUNDIALS_VERSION_PATCH preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_SUNDIALS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_SUNDIALS
  constexpr const unsigned int sundials_version_patch =
    DEAL_II_SUNDIALS_VERSION_PATCH;
#else
  constexpr const unsigned int sundials_version_patch                   = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_TRILINOS_VERSION_MAJOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_TRILINOS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_TRILINOS
  constexpr const unsigned int trilinos_version_major =
    DEAL_II_TRILINOS_VERSION_MAJOR;
#else
  constexpr const unsigned int trilinos_version_major                   = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_TRILINOS_VERSION_MINOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_TRILINOS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_TRILINOS
  constexpr const unsigned int trilinos_version_minor =
    DEAL_II_TRILINOS_VERSION_MINOR;
#else
  constexpr const unsigned int trilinos_version_minor                   = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_TRILINOS_VERSION_SUBMINOR preprocessor symbol as determined by
   * CMake and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_TRILINOS preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_TRILINOS
  constexpr const unsigned int trilinos_version_subminor =
    DEAL_II_TRILINOS_VERSION_SUBMINOR;
#else
  constexpr const unsigned int trilinos_version_subminor                = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_CGAL_VERSION_MAJOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_CGAL preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_CGAL
  constexpr const unsigned int cgal_version_major = DEAL_II_CGAL_VERSION_MAJOR;
#else
  constexpr const unsigned int cgal_version_major                       = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_CGAL_VERSION_MINOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_CGAL preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_CGAL
  constexpr const unsigned int cgal_version_minor = DEAL_II_CGAL_VERSION_MINOR;
#else
  constexpr const unsigned int cgal_version_minor                       = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_CGAL_VERSION_SUBMINOR preprocessor symbol as determined by CMake
   * and stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_CGAL preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_CGAL
  constexpr const unsigned int cgal_version_subminor =
    DEAL_II_CGAL_VERSION_SUBMINOR;
#else
  constexpr const unsigned int cgal_version_subminor                    = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_MPI_VERSION_MAJOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_MPI preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_MPI
  constexpr const unsigned int mpi_version_major = DEAL_II_MPI_VERSION_MAJOR;
#else
  constexpr const unsigned int mpi_version_major                        = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_MPI_VERSION_MINOR preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or 0 if the
   * DEAL_II_WITH_MPI preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_MPI
  constexpr const unsigned int mpi_version_minor = DEAL_II_MPI_VERSION_MINOR;
#else
  constexpr const unsigned int mpi_version_minor                        = 0;
#endif

  /**
   * A variable that equals the preprocessor value of the
   * DEAL_II_GMSH_EXECUTABLE_PATH preprocessor symbol as determined by CMake and
   * stored in the include/deal.II/base/config.h file, or the empty string if
   * the DEAL_II_WITH_GMSH preprocessor symbol is not defined there.
   */
#ifdef DEAL_II_WITH_GMSH
  constexpr const char gmsh_executable_path[] = DEAL_II_GMSH_EXECUTABLE_PATH;
#else
  constexpr const char         gmsh_executable_path[]                   = "";
#endif

  /**
   * Equivalent of `DEAL_II_VERSION_GTE(major, minor, subminor)`.
   */
  constexpr bool
  dealii_version_gte(const unsigned int major,
                     const unsigned int minor,
                     const unsigned int subminor)
  {
    return ((dealii_version_major * 10000 + dealii_version_minor * 100 +
             dealii_version_subminor) >=
            (major * 10000 + minor * 100 + subminor));
  }

  /**
   * Equivalent of `DEAL_II_BOOST_VERSION_GTE(major, minor, subminor)`.
   */
  constexpr bool
  boost_version_gte(const unsigned int major,
                    const unsigned int minor,
                    const unsigned int subminor)
  {
    return ((boost_version_major * 100000 + boost_version_minor * 100 +
             boost_version_subminor) >=
            (major * 100000 + minor * 100 + subminor));
  }

  /**
   * Equivalent of `DEAL_II_KOKKOS_VERSION_GTE(major, minor, subminor)`; returns
   * `false` when `DEAL_II_WITH_KOKKOS` is not defined.
   */
  constexpr bool
  kokkos_version_gte(const unsigned int major,
                     const unsigned int minor,
                     const unsigned int subminor)
  {
    return (
      with_kokkos &&
      ((kokkos_version_major * 10000 + kokkos_version_minor * 100 +
        kokkos_version_subminor) >= (major * 10000 + minor * 100 + subminor)));
  }

  /**
   * Equivalent of `DEAL_II_OPENCASCADE_VERSION_GTE(major, minor, subminor)`;
   * returns `false` when `DEAL_II_WITH_OPENCASCADE` is not defined.
   */
  constexpr bool
  opencascade_version_gte(const unsigned int major,
                          const unsigned int minor,
                          const unsigned int subminor)
  {
    return (with_opencascade &&
            ((opencascade_version_major * 10000 +
              opencascade_version_minor * 100 + opencascade_version_subminor) >=
             (major * 10000 + minor * 100 + subminor)));
  }

  /**
   * Equivalent of `DEAL_II_P4EST_VERSION_GTE(major, minor, subminor, patch)`;
   * returns `false` when `DEAL_II_WITH_P4EST` is not defined.
   */
  constexpr bool
  p4est_version_gte(const unsigned int major,
                    const unsigned int minor,
                    const unsigned int subminor,
                    const unsigned int patch)
  {
    return (with_p4est &&
            ((p4est_version_major * 1000000 + p4est_version_minor * 10000 +
              p4est_version_subminor * 100 + p4est_version_patch) >=
             (major * 1000000 + minor * 10000 + subminor * 100 + patch)));
  }

  /**
   * Equivalent of `DEAL_II_SUNDIALS_VERSION_GTE(major, minor, patch)`; returns
   * `false` when `DEAL_II_WITH_SUNDIALS` is not defined.
   */
  constexpr bool
  sundials_version_gte(const unsigned int major,
                       const unsigned int minor,
                       const unsigned int patch)
  {
    return (
      with_sundials &&
      ((sundials_version_major * 10000 + sundials_version_minor * 100 +
        sundials_version_patch) >= (major * 10000 + minor * 100 + patch)));
  }

  /**
   * Equivalent of `DEAL_II_SUNDIALS_VERSION_LT(major, minor, patch)`; returns
   * `false` when `DEAL_II_WITH_SUNDIALS` is not defined.
   */
  constexpr bool
  sundials_version_lt(const unsigned int major,
                      const unsigned int minor,
                      const unsigned int patch)
  {
    return (with_sundials &&
            ((sundials_version_major * 10000 + sundials_version_minor * 100 +
              sundials_version_patch) < (major * 10000 + minor * 100 + patch)));
  }

  /**
   * Equivalent of `DEAL_II_PETSC_VERSION_LT(major, minor, subminor)`; returns
   * `false` when `DEAL_II_WITH_PETSC` is not defined.
   */
  constexpr bool
  petsc_version_lt(const unsigned int major,
                   const unsigned int minor,
                   const unsigned int subminor)
  {
#ifdef DEAL_II_WITH_PETSC
    return PETSC_VERSION_LT(major, minor, subminor);
#else
    (void)major;
    (void)minor;
    (void)subminor;
    return false;
#endif
  }

  /**
   * Equivalent of `DEAL_II_PETSC_VERSION_GTE(major, minor, subminor)`; returns
   * `false` when `DEAL_II_WITH_PETSC` is not defined.
   */
  constexpr bool
  petsc_version_gte(const unsigned int major,
                    const unsigned int minor,
                    const unsigned int subminor)
  {
#ifdef DEAL_II_WITH_PETSC
    return PETSC_VERSION_GE(major, minor, subminor);
#else
    (void)major;
    (void)minor;
    (void)subminor;
    return false;
#endif
  }

  /**
   * Equivalent of `DEAL_II_SLEPC_VERSION_LT(major, minor, subminor)`; returns
   * `false` when `DEAL_II_WITH_SLEPC` is not defined.
   */
  constexpr bool
  slepc_version_lt(const unsigned int major,
                   const unsigned int minor,
                   const unsigned int subminor)
  {
#ifdef DEAL_II_WITH_SLEPC
    return SLEPC_VERSION_LT(major, minor, subminor);
#else
    (void)major;
    (void)minor;
    (void)subminor;
    return false;
#endif
  }

  /**
   * Equivalent of `DEAL_II_SLEPC_VERSION_GTE(major, minor, subminor)`; returns
   * `false` when `DEAL_II_WITH_SLEPC` is not defined.
   */
  constexpr bool
  slepc_version_gte(const unsigned int major,
                    const unsigned int minor,
                    const unsigned int subminor)
  {
#ifdef DEAL_II_WITH_SLEPC
    return SLEPC_VERSION_GE(major, minor, subminor);
#else
    (void)major;
    (void)minor;
    (void)subminor;
    return false;
#endif
  }

  /**
   * Equivalent of `DEAL_II_TRILINOS_VERSION_GTE(major, minor, subminor)`;
   * returns `false` when `DEAL_II_WITH_TRILINOS` is not defined.
   */
  constexpr bool
  trilinos_version_gte(const unsigned int major,
                       const unsigned int minor,
                       const unsigned int subminor)
  {
    return (with_trilinos &&
            ((trilinos_version_major * 10000 + trilinos_version_minor * 100 +
              trilinos_version_subminor) >=
             (major * 10000 + minor * 100 + subminor)));
  }

  /**
   * Equivalent of `DEAL_II_CGAL_VERSION_GTE(major, minor, subminor)`; returns
   * `false` when `DEAL_II_WITH_CGAL` is not defined.
   */
  constexpr bool
  cgal_version_gte(const unsigned int major,
                   const unsigned int minor,
                   const unsigned int subminor)
  {
    return (with_cgal && ((cgal_version_major * 10000 +
                           cgal_version_minor * 100 + cgal_version_subminor) >=
                          (major * 10000 + minor * 100 + subminor)));
  }

  /**
   * Equivalent of `DEAL_II_MPI_VERSION_GTE(major, minor)`; returns `false` when
   * `DEAL_II_WITH_MPI` is not defined.
   */
  constexpr bool
  mpi_version_gte(const unsigned int major, const unsigned int minor)
  {
    return (with_mpi && ((mpi_version_major * 100 + mpi_version_minor) >=
                         (major * 100 + minor)));
  }
} // namespace configuration

DEAL_II_NAMESPACE_CLOSE

#endif
