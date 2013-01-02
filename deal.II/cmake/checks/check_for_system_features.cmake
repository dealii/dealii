#####
##
## Copyright (C) 2012 by the deal.II authors
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

#
# Check whether the compiler allows to use arithmetic operations
# +-*/ on vectorized data types or whether we need to use
# _mm_add_pd for addition and so on. +-*/ is preferred because
# it allows the compiler to choose other optimizations like
# fused multiply add, whereas _mm_add_pd explicitly enforces the
# assembler command.
#
# - Matthias Maier, rewritten 2012
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <emmintrin.h>
  int main()
  {
    __m128d a, b;
    a = _mm_set_sd (1.0);
    b = _mm_set1_pd (2.1);
    __m128d c = a + b;
    __m128d d = b - c;
    __m128d e = c * a + d;
    __m128d f = e/a;
    (void)f;
  }
  "
  DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS)


###########################################################################
#                                                                         #
#                     POSIX and Linux specific tests:                     #
#                                                                         #
###########################################################################

CHECK_INCLUDE_FILE("sys/resource.h"  HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE("sys/time.h" HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE("sys/times.h" HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)


#
# Check for various posix specific functions. On a posix system they should
# be all defined in unistd.h. On other platforms, most notably
# Windows/MinGW unistd.h is available but not all posix functions. So test
# for each funtion as well.
#
CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)
CHECK_FUNCTION_EXISTS(gethostname HAVE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS(getpid HAVE_GETPID)
CHECK_FUNCTION_EXISTS(rand_r HAVE_RAND_R)
CHECK_FUNCTION_EXISTS(times HAVE_TIMES)


#
# Do we have the Bessel function jn?
#
FIND_LIBRARY(m_lib NAMES m)
MARK_AS_ADVANCED(m_lib)

IF(NOT m_lib MATCHES "-NOTFOUND")
  SET(CMAKE_REQUIRED_LIBRARIES ${m_lib})
  CHECK_FUNCTION_EXISTS(jn HAVE_JN)
  SET(CMAKE_REQUIRED_LIBRARIES)
  IF(HAVE_JN)
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${m_lib})
  ENDIF()
ENDIF()


###########################################################################
#                                                                         #
#                    Windos and CYGWIN specific setup:                    #
#                                                                         #
###########################################################################

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  #
  # Export DEAL_II_MSVC if we are on a Windows platform:
  #
  SET(DEAL_II_MSVC TRUE)

  #
  # Disable -ggdb and -g on Windows/MinGW targets for the moment until the
  # compilation issues with too big files is resolved
  #
  # - Matthias Maier, 2012
  #
  STRIP_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-ggdb")
  STRIP_FLAG(DEAL_II_SHARED_LINKER_FLAGS_DEBUG "-ggdb")
  STRIP_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-g")
  STRIP_FLAG(DEAL_II_SHARED_LINKER_FLAGS_DEBUG "-g")
ENDIF()


#
# Disable shared libraries on CYGWIN and Windows targets for the moment.
# Our support for shared libraries on Windows is a bit buggy atm..
#
# - Matthias Maier, 2012
#
IF( CMAKE_SYSTEM_NAME MATCHES "CYGWIN" OR
    CMAKE_SYSTEM_NAME MATCHES "Windows" )
  MESSAGE(WARNING "\n"
    "BUILD_SHARED_LIBS forced to OFF\n\n"
    )
  SET(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

