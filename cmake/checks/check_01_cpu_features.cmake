## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


########################################################################
#                                                                      #
#                   Platform and CPU specific tests:                   #
#                                                                      #
########################################################################

#
# This file sets up
#
#   DEAL_II_WORDS_BIGENDIAN
#   DEAL_II_HAVE_SSE2                    (*)
#   DEAL_II_HAVE_AVX                     (*)
#   DEAL_II_HAVE_AVX512                  (*)
#   DEAL_II_HAVE_ALTIVEC                 (*)
#   DEAL_II_HAVE_ARM_NEON                (*)
#   DEAL_II_HAVE_OPENMP_SIMD             (*)
#   DEAL_II_VECTORIZATION_WIDTH_IN_BITS
#   DEAL_II_OPENMP_SIMD_PRAGMA
#
# (*) It is possible to manually set the above values to their
# corresponding values, when platform introspection is disabled with
# DEAL_II_ALLOW_PLATFORM_INTROSPECTION=OFF,
#


#
# Determine the Endianness of the platform:
#
if(CMAKE_C_COMPILER_WORKS)
  include(TestBigEndian)

  clear_cmake_required()
  TEST_BIG_ENDIAN(DEAL_II_WORDS_BIGENDIAN)
  reset_cmake_required()
else()
  message(STATUS
    "No suitable C compiler was found! Assuming little endian platform."
    )
  set(DEAL_II_WORDS_BIGENDIAN "0")
endif()


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

if(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
  #
  # Take care that the following tests are rerun if the
  # CMAKE_REQUIRED_FLAGS changes..
  #
  unset_if_changed(CHECK_CPU_FEATURES_FLAGS_SAVED "${CMAKE_REQUIRED_FLAGS}"
    DEAL_II_HAVE_SSE2 DEAL_II_HAVE_AVX DEAL_II_HAVE_AVX512 DEAL_II_HAVE_ALTIVEC DEAL_II_HAVE_ARM_NEON
    )

  CHECK_CXX_SOURCE_RUNS(
    "
    #include <x86intrin.h>
    int main()
    {
    __m128d a, b;
    const unsigned int vector_bytes = sizeof(__m128d);
    const int n_vectors = vector_bytes/sizeof(double);
    __m128d * data =
      reinterpret_cast<__m128d*>(_mm_malloc (2*vector_bytes, vector_bytes));
    double * ptr = reinterpret_cast<double*>(&a);
    ptr[0] = static_cast<volatile double>(1.0);
    for (int i=1; i<n_vectors; ++i)
      ptr[i] = 0.0;
    b = _mm_set1_pd (static_cast<volatile double>(2.25));
    data[0] = _mm_add_pd (a, b);
    data[1] = _mm_mul_pd (b, data[0]);
    ptr = reinterpret_cast<double*>(&data[1]);
    int return_value = 0;
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
    #include <x86intrin.h>
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
      ptr[0] = static_cast<volatile double>(1.0);
      for (int i=1; i<n_vectors; ++i)
        ptr[i] = 0.0;
      b = _mm256_set1_pd (static_cast<volatile double>(2.25));
      data[0] = _mm256_add_pd (a, b);
      data[1] = _mm256_mul_pd (b, data[0]);
      ptr = reinterpret_cast<double*>(&data[1]);
      int return_value = 0;
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
    #include <x86intrin.h>
    int main()
    {
      __m512d a, b;
      const unsigned int vector_bytes = sizeof(__m512d);
      const int n_vectors = vector_bytes/sizeof(double);
      __m512d * data =
        reinterpret_cast<__m512d*>(_mm_malloc (2*vector_bytes, vector_bytes));
      double * ptr = reinterpret_cast<double*>(&a);
      ptr[0] = static_cast<volatile double>(1.0);
      for (int i=1; i<n_vectors; ++i)
        ptr[i] = 0.0;
      const volatile double x = 2.25;
      b = _mm512_set1_pd(x);
      data[0] = _mm512_add_pd (a, b);
      data[1] = _mm512_mul_pd (b, data[0]);
      ptr = reinterpret_cast<double*>(&data[1]);
      int return_value = 0;
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

  CHECK_CXX_SOURCE_RUNS(
    "
    #ifndef __ALTIVEC__
    #error \"__ALTIVEC__ flag not set, no support for Altivec\"
    #endif
    #include <altivec.h>
    #undef vector
    #undef pixel
    #undef bool
    int main()
    {
    __vector double a, b, data1, data2;
    const int n_vectors = sizeof(a)/sizeof(double);
    double * ptr = reinterpret_cast<double*>(&a);
    ptr[0] = static_cast<volatile double>(1.0);
    for (int i=1; i<n_vectors; ++i)
      ptr[i] = 0.0;
    b = vec_splats (static_cast<volatile double>(2.25));
    data1 = vec_add (a, b);
    data2 = vec_mul (b, data1);
    ptr = reinterpret_cast<double*>(&data2);
    int return_value = 0;
    if (ptr[0] != 7.3125)
      return_value += 1;
    for (int i=1; i<n_vectors; ++i)
      if (ptr[i] != 5.0625)
        return_value += 2;
    b = vec_splats (static_cast<volatile double>(-1.0));
    data1 = vec_abs(vec_mul (b, data2));
    vec_vsx_st(data1, 0, ptr);
    b = vec_vsx_ld(0, ptr);
    ptr = reinterpret_cast<double*>(&b);
    if (ptr[0] != 7.3125)
      return_value += 4;
    for (int i=1; i<n_vectors; ++i)
      if (ptr[i] != 5.0625)
        return_value += 8;
    return return_value;
    }
    "
    DEAL_II_HAVE_ALTIVEC)

    CHECK_CXX_SOURCE_RUNS(
      "
      #ifndef __ARM_NEON
      #error Preprocessor flag not found
      #endif
      #include <arm_neon.h>
      #include <stdlib.h>
      int main()
      {
      float64x2_t a, b;
      const unsigned int vector_bytes = sizeof(float64x2_t);
      const int n_vectors = vector_bytes/sizeof(double);
      float64x2_t * data =
        reinterpret_cast<float64x2_t*>(aligned_alloc (vector_bytes, 2*vector_bytes));
      double * ptr = reinterpret_cast<double*>(&a);
      ptr[0] = static_cast<volatile double>(1.0);
      for (int i=1; i<n_vectors; ++i)
        ptr[i] = 0.0;
      b = vdupq_n_f64 (static_cast<volatile double>(2.25));
      data[0] = vaddq_f64 (a, b);
      data[1] = vmulq_f64 (b, data[0]);
      ptr = reinterpret_cast<double*>(&data[1]);
      int return_value = 0;
      if (ptr[0] != 7.3125)
        return_value = 1;
      for (int i=1; i<n_vectors; ++i)
        if (ptr[i] != 5.0625)
          return_value = 1;
      free(data);
      return return_value;
      }
      "
      DEAL_II_HAVE_ARM_NEON)

  #
  # OpenMP 4.0 can be used for vectorization. Only the vectorization
  # instructions are allowed, the threading must be done through TBB.
  #

  #
  # Choosing the right compiler flag is a bit of a mess:
  #
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "15" )
      set(_keyword "qopenmp")
    elseif("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "14" )
      set(_keyword "openmp")
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(_keyword "openmp")
  else()
    set(_keyword "fopenmp")
  endif()

  CHECK_CXX_COMPILER_FLAG("-${_keyword}-simd" DEAL_II_HAVE_OPENMP_SIMD)

endif() # IF DEAL_II_ALLOW_PLATFORM_INTROSPECTION


#
# Choose DEAL_II_COMPILER_VECTORIZATION level depending on AVX support
# (that was autodetected or manually specified).
#

if(DEAL_II_HAVE_AVX512)
  set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 512)
elseif(DEAL_II_HAVE_AVX)
  set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 256)
elseif(DEAL_II_HAVE_SSE2)
  set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 128)
else()
  set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 0)
endif()

if(DEAL_II_HAVE_ALTIVEC)
  set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 128)
endif()

if(DEAL_II_HAVE_ARM_NEON)
  set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 128)
endif()

#
# If we have OpenMP SIMD support (i.e. DEAL_II_HAVE_OPENMP_SIMD is true)
# populate DEAL_II_OPENMP_SIMD_PRAGMA.
#

set(DEAL_II_OPENMP_SIMD_PRAGMA " ")
if(DEAL_II_HAVE_OPENMP_SIMD)
  add_flags(DEAL_II_CXX_FLAGS "-${_keyword}-simd")
  # Intel is special:
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_flags(DEAL_II_LINKER_FLAGS "-${_keyword}")
  endif()
  set(DEAL_II_OPENMP_SIMD_PRAGMA "_Pragma(\"omp simd\")")
endif()
