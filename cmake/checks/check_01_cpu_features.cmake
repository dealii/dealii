## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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
#   DEAL_II_COMPILER_VECTORIZATION_LEVEL
#   DEAL_II_HAVE_OPENMP_SIMD             *)
#   DEAL_II_OPENMP_SIMD_PRAGMA
#
# *)
# It is is possible to manually set the above values to their corresponding
# values, when platform introspection is disabled with
# DEAL_II_ALLOW_PLATFORM_INTROSPECTION=OFF,
#


#
# Determine the Endianness of the platform:
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
  # Take care that the following tests are rerun if the
  # CMAKE_REQUIRED_FLAGS changes..
  #
  UNSET_IF_CHANGED(CHECK_CPU_FEATURES_FLAGS_SAVED "${CMAKE_REQUIRED_FLAGS}"
    DEAL_II_HAVE_SSE2 DEAL_II_HAVE_AVX DEAL_II_HAVE_AVX512
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

  #
  # clang-3.6.0 has a bug in operator+ on two VectorizedArray components as
  # defined in deal.II. Therefore, the test for AVX needs to also test for
  # operator+ to be correctly implemented.
  #
  CHECK_CXX_SOURCE_RUNS(
    "
    #ifndef __AVX__
    #error \"__AVX__ flag not set, no support for AVX\"
    #endif
    #include <immintrin.h>
    class VectorizedArray
    {
    public:
      VectorizedArray &
      operator += (const VectorizedArray &vec)
      {
        data = _mm256_add_pd (data, vec.data);
        return *this;
      }
      __m256d data;
    };
    inline
    VectorizedArray
    operator + (const VectorizedArray &u, const VectorizedArray &v)
    {
      VectorizedArray tmp = u;
      return tmp+=v;
    }
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
      VectorizedArray c, d, e;
      c.data = b;
      d.data = b;
      e = c + d;
      ptr = reinterpret_cast<double*>(&e.data);
      for (int i=0; i<n_vectors; ++i)
        if (ptr[i] != 4.5)
          return_value = 1;
      _mm_free (data);
      return return_value;
    }
    "
    DEAL_II_HAVE_AVX)

  CHECK_CXX_SOURCE_RUNS(
    "
    #ifndef __AVX512F__
    #error \"__AVX512F__ flag not set, no support for AVX512\"
    #endif
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
      b = _mm512_set1_pd(x);
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
# OpenMP 4.0 can be used for vectorization. Only the vectorization
# instructions are allowed, the threading must be done through TBB.
#

#
# Choosing the right compiler flag is a bit of a mess:
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  IF("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "15" )
    SET(_keyword "qopenmp")
  ELSEIF("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "14" )
    SET(_keyword "openmp")
  ENDIF()
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(_keyword "openmp")
ELSE()
  SET(_keyword "fopenmp")
ENDIF()

CHECK_CXX_COMPILER_FLAG("-${_keyword}-simd" DEAL_II_HAVE_OPENMP_SIMD)

SET(DEAL_II_OPENMP_SIMD_PRAGMA " ")
IF(DEAL_II_HAVE_OPENMP_SIMD)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "-${_keyword}-simd")
  # Intel is special:
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    ADD_FLAGS(DEAL_II_LINKER_FLAGS "-${_keyword}")
  ENDIF()
  SET(DEAL_II_OPENMP_SIMD_PRAGMA "_Pragma(\"omp simd\")")
ENDIF()
