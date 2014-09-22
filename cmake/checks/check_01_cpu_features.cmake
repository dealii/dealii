## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------


########################################################################
#                                                                      #
#                   Platform and CPU specific tests:                   #
#                                                                      #
########################################################################

#
# This file sets up
#
#   DEAL_II_WORDS_BIGENDIAN
#   DEAL_II_HAVE_SSE2                    *)
#   DEAL_II_HAVE_AVX                     *)
#   DEAL_II_HAVE_AVX512                  *)
#   DEAL_II_HAVE_OPENMP_SIMD             *)
#   DEAL_II_COMPILER_VECTORIZATION_LEVEL
#
# *)
# It is is possible to manually set the above values to their corresponding
# values, when platform introspection is disabled with
# DEAL_II_ALLOW_PLATFORM_INTROSPECTION=OFF,
#


#
# Determine the Endianess of the platform:
#
IF(CMAKE_C_COMPILER_WORKS)
  INCLUDE(TestBigEndian)

  CLEAR_CMAKE_REQUIRED()
  TEST_BIG_ENDIAN(DEAL_II_WORDS_BIGENDIAN)
  RESET_CMAKE_REQUIRED()
ELSE()
  MESSAGE(STATUS
    "No suitable C compiler was found! Assuming little endian platform."
    )
  SET(DEAL_II_WORDS_BIGENDIAN "0")
ENDIF()


#
# Check whether the compiler allows for vectorization and that
# vectorization actually works on the given CPU. For this test, we use
# compiler intrinsics similar to what is used in the deal.II library and
# check whether the arithmetic operations are correctly performed on
# examples where all numbers are exactly represented as floating point
# numbers.
#
# - Matthias Maier, rewritten 2012
#

IF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
  #
  # Take care that the following tests are rerun if CMAKE_REQUIRED_FLAGS
  # changes..
  #
  IF(NOT "${CMAKE_REQUIRED_FLAGS}" STREQUAL "${DEAL_II_CHECK_CPU_FEATURES_SAVED}")
    UNSET(DEAL_II_HAVE_SSE2 CACHE)
    UNSET(DEAL_II_HAVE_AVX CACHE)
    UNSET(DEAL_II_HAVE_AVX512 CACHE)
  ENDIF()
  SET(DEAL_II_CHECK_CPU_FEATURES_SAVED
    "${CMAKE_REQUIRED_FLAGS}" CACHE INTERNAL "" FORCE
    )

  CHECK_CXX_SOURCE_RUNS(
    "
    #include <emmintrin.h>
    int main()
    {
    __m128d a, b;
    const unsigned int vector_bytes = sizeof(__m128d);
    const int n_vectors = vector_bytes/sizeof(double);
    __m128d * data =
      reinterpret_cast<__m128d*>(_mm_malloc (2*vector_bytes, vector_bytes));
    double * ptr = reinterpret_cast<double*>(&a);
    ptr[0] = (volatile double)(1.0);
    for (int i=1; i<n_vectors; ++i)
      ptr[i] = 0.0;
    b = _mm_set1_pd ((volatile double)(2.25));
    data[0] = _mm_add_pd (a, b);
    data[1] = _mm_mul_pd (b, data[0]);
    ptr = reinterpret_cast<double*>(&data[1]);
    unsigned int return_value = 0;
    if (ptr[0] != 7.3125)
      return_value = 1;
    for (int i=1; i<n_vectors; ++i)
      if (ptr[i] != 5.0625)
        return_value = 1;
    _mm_free (data);
    return return_value;
    }
    "
    DEAL_II_HAVE_SSE2)

  CHECK_CXX_SOURCE_RUNS(
    "
    #include <immintrin.h>
    int main()
    {
      __m256d a, b;
      const unsigned int vector_bytes = sizeof(__m256d);
      const int n_vectors = vector_bytes/sizeof(double);
      __m256d * data =
        reinterpret_cast<__m256d*>(_mm_malloc (2*vector_bytes, vector_bytes));
      double * ptr = reinterpret_cast<double*>(&a);
      ptr[0] = (volatile double)(1.0);
      for (int i=1; i<n_vectors; ++i)
        ptr[i] = 0.0;
      b = _mm256_set1_pd ((volatile double)(2.25));
      data[0] = _mm256_add_pd (a, b);
      data[1] = _mm256_mul_pd (b, data[0]);
      ptr = reinterpret_cast<double*>(&data[1]);
      unsigned int return_value = 0;
      if (ptr[0] != 7.3125)
        return_value = 1;
      for (int i=1; i<n_vectors; ++i)
        if (ptr[i] != 5.0625)
          return_value = 1;
      _mm_free (data);
      return return_value;
    }
    "
    DEAL_II_HAVE_AVX)

  CHECK_CXX_SOURCE_RUNS(
    "
    #include <immintrin.h>
    int main()
    {
      __m512d a, b;
      const unsigned int vector_bytes = sizeof(__m512d);
      const int n_vectors = vector_bytes/sizeof(double);
      __m512d * data =
        reinterpret_cast<__m512d*>(_mm_malloc (2*vector_bytes, vector_bytes));
      double * ptr = reinterpret_cast<double*>(&a);
      ptr[0] = (volatile double)(1.0);
      for (int i=1; i<n_vectors; ++i)
        ptr[i] = 0.0;
      const volatile double x = 2.25;
      b = _mm512_set_pd(x, x, x, x, x, x, x, x);
      data[0] = _mm512_add_pd (a, b);
      data[1] = _mm512_mul_pd (b, data[0]);
      ptr = reinterpret_cast<double*>(&data[1]);
      unsigned int return_value = 0;
      if (ptr[0] != 7.3125)
        return_value = 1;
      for (int i=1; i<n_vectors; ++i)
        if (ptr[i] != 5.0625)
          return_value = 1;
      _mm_free (data);
      return return_value;
    }
    "
    DEAL_II_HAVE_AVX512)
ENDIF()

IF(DEAL_II_HAVE_AVX512)
  SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 3)
ELSEIF(DEAL_II_HAVE_AVX)
  SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 2)
ELSEIF(DEAL_II_HAVE_SSE2)
  SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 1)
ELSE()
  SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 0)
ENDIF()


#
# OpenMP 4.0 can be used for vectorization (supported by gcc-4.9.1 and
# later). Only the vectorization instructions
# are allowed, the threading must done through TBB.
#

IF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    CHECK_CXX_COMPILER_FLAG("-openmp-simd" DEAL_II_HAVE_OPENMP_SIMD)
  ELSE()
    CHECK_CXX_COMPILER_FLAG("-fopenmp-simd" DEAL_II_HAVE_OPENMP_SIMD)
  ENDIF()
ENDIF()

IF(DEAL_II_HAVE_OPENMP_SIMD)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    ADD_FLAGS(DEAL_II_CXX_FLAGS "-openmp-simd")
    ADD_FLAGS(DEAL_II_LINKER_FLAGS "-openmp")
  ELSE()
    ADD_FLAGS(DEAL_II_CXX_FLAGS "-fopenmp-simd")
  ENDIF()
  SET(DEAL_II_OPENMP_SIMD_PRAGMA "_Pragma(\"omp simd\")")
ELSE()
  SET(DEAL_II_OPENMP_SIMD_PRAGMA " ")
ENDIF()
