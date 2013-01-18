#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####


###########################################################################
#                                                                         #
#                    Platform and CPU specific tests:                     #
#                                                                         #
###########################################################################

#
# This file sets up
#
#   DEAL_II_WORDS_BIGENDIAN
#   DEAL_II_HAVE_SSE2                    *)
#   DEAL_II_HAVE_AVX                     *)
#   DEAL_II_COMPILER_VECTORIZATION_LEVEL
#
# Note: It is possible to disable the platform introspection tests (e.g.
# for cross compiling or for packaging) by defining the variables (*) in
# the cache prior to configure
#


#
# Determine the Endianess of the platform:
#
INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(DEAL_II_WORDS_BIGENDIAN)


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
CHECK_CXX_SOURCE_RUNS(
  "
  #include <emmintrin.h>
  #include <mm_malloc.h>
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
  #include <mm_malloc.h>
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


IF(DEAL_II_HAVE_SSE2)
  IF(DEAL_II_HAVE_AVX)
    SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 2)
  ELSE()
    SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 1)
  ENDIF()
ELSE()
  SET(DEAL_II_COMPILER_VECTORIZATION_LEVEL 0)
ENDIF()

