
#
# Check for various compiler features.
#



#
# Check whether the std::vector::iterator is just a plain pointer
#
# (Yes. It is not a bug. But the logic is the same.)
CHECK_CXX_COMPILER_BUG(
  "
  #include <vector>
  template <typename T> void f(T) {}
  template void f(int *);
  template void f(std::vector<int>::iterator);
  int main(){return 0;}
  "
  DEAL_II_VECTOR_ITERATOR_IS_POINTER)



#
# Check for existence of the __builtin_expect facility of newer
# gcc compilers. This can be used to hint the compiler's branch
# prediction unit in some cases. We use it in the AssertThrow
# macros.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  bool f() {}
  int main(){ if (__builtin_expect(f(),false)) ; }
  "
  HAVE_BUILTIN_EXPECT)



#
# Newer versions of gcc have a very nice feature: you can set
# a verbose terminate handler, that not only aborts a program
# when an exception is thrown and not caught somewhere, but
# before aborting it prints that an exception has been thrown,
# and possibly what the std::exception::what() function has to
# say. Since many people run into the trap of not having a
# catch clause in main(), they wonder where that abort may be
# coming from.  The terminate handler then at least says what is
# missing in their program.
#
# This test checks whether this feature is available.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <exception>
  namespace __gnu_cxx
  {
    extern void __verbose_terminate_handler ();
  }
  struct preload_terminate_dummy
  {
    preload_terminate_dummy()
    {
      std::set_terminate (__gnu_cxx::__verbose_terminate_handler);
    }
  };
  static preload_terminate_dummy dummy;
  int main() { throw 1; return 0; }
  "
  HAVE_VERBOSE_TERMINATE)



#
# Check whether glibc-like stacktrace information is available
# for the Exception class. If it is, then try to also determine
# whether the compiler accepts the -rdynamic flag, since that is
# recommended for linking if one wants to have meaningful
# backtraces.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <execinfo.h>
  #include <stdlib.h>
  void * array[25];
  int nSize = backtrace(array, 25);
  char ** symbols = backtrace_symbols(array, nSize);
  int main(){ free(symbols); return 0; }
  "
  HAVE_GLIBC_STACKTRACE)

#
# On Mac OS X, -rdynamic is accepted by the compiler (i.e.
# it doesn't produce an error) but we always get a warning
# that it isn't supported. That's pretty stupid because
# we can't test for it. Consequently, only run the test
# if not on OS X.
#
IF(HAVE_GLIBC_STACKTRACE AND NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")

  ENABLE_IF_SUPPORTED(CMAKE_SHARED_LINKER_FLAGS "-rdynamic")

ENDIF()



#
# Check whether the compiler offers a way to demangle symbols
# from within the program. Used inside the exception stacktrace
# mechanism.
#
# The example code is taken from
#   http://gcc.gnu.org/onlinedocs/libstdc++/18_support/howto.html#6
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <exception>
  #include <iostream>
  #include <cxxabi.h>
  #include <cstdlib>

  struct empty { };

  template <typename T, int N>
  struct bar { };

  int     status;
  char   *realname;

  int main()
  {
    // exception classes not in <stdexcept>, thrown by the implementation
    // instead of the user
    std::bad_exception  e;
    realname = abi::__cxa_demangle(e.what(), 0, 0, &status);
    free(realname);


    // typeid
    bar<empty,17>          u;
    const std::type_info  &ti = typeid(u);

    realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
    free(realname);

      return 0;
  }
  "
  HAVE_LIBSTDCXX_DEMANGLER)



#
# Check whether the compiler allows for vectorization and that
# vectorization actually works. For this test, we use compiler
# intrinsics similar to what is used in the deal.II library and
# check whether the arithmetic operations are correctly performed
# on examples where all numbers are exactly represented as
# floating point numbers.
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
  DEAL_II_COMPILER_HAS_SSE2)

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
  DEAL_II_COMPILER_HAS_AVX)

IF(DEAL_II_COMPILER_HAS_SSE2)
  IF(DEAL_II_COMPILER_HAS_AVX)
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



#
# Check if the declared prototype of abort() has a throw()
# specification. We overload abort() in our testsuite, so have
# to make sure that we match the exception specification
# correctly.
#
# (Yes. It is not a bug. But the logic is the same.)
CHECK_CXX_COMPILER_BUG(
  "
  #include <cstdlib>
  extern \"C\" void abort () { for(;;) ; }
  int main(){ return 0; }
  "
  DEAL_II_ABORT_NOTHROW_EXCEPTION)



#
# Gcc and some other compilers have __PRETTY_FUNCTION__, showing
# an unmangled version of the function we are presently in,
# while __FUNCTION__ (or __func__ in ISO C99) simply give the
# function name which would not include the arguments of that
# function, leading to problems in C++ with overloaded function
# names.
#
# If __PRETTY_FUNCTION__ is not available, try to find out whether
# __func__ is available and use the preprocessor to set the first
# thing to the second. If this is also not the case, then set it
# to something indicating non-availability.
#

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <iostream>
  int main()
  {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    return 0;
  }
  "
  DEAL_II_COMPILER_HAS_ATTRIBUTE_PRETTY_FUNCTION)

IF(NOT DEAL_II_COMPILER_HAS_ATTRIBUTE_PRETTY_FUNCTION)
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <iostream>
    int main()
    {
      std::cout << __func__ << std::endl;
      return 0;
    }
    "
    DEAL_II_COMPILER_HAS_ATTRIBUTE_FUNC)

  IF(DEAL_II_COMPILER_HAS_ATTRIBUTE_FUNC)
    SET(__PRETTY_FUNCTION__ "__func")
  ELSE()
    SET(__PRETTY_FUNCTION__ "(not available)")
  ENDIF()

ENDIF()



#
# Check for minimal vector capacity
#
GET_CXX_SOURCE_RETURN_VALUE(
  "
  #include <vector>
  int main () {
    std::vector<int> v(1);
    v.reserve (1);
    v.resize (1);
    return v.capacity();
  }
  "
  DEAL_II_MIN_VECTOR_CAPACITY
  DEAL_II_MIN_VECTOR_CAPACITY_RETURN_VALUE)

IF(NOT DEAL_II_MIN_VECTOR_CAPACITY)
  # We have a problem...
  MESSAGE(WARNING
    "Could not determine DEAL_II_MIN_VECTOR_CAPACITY, "
    "source might not compile..."
    )
ELSE()
  SET(DEAL_II_MIN_VECTOR_CAPACITY
    ${DEAL_II_MIN_VECTOR_CAPACITY_RETURN_VALUE}
    )
ENDIF()



#
# Do same thing with std::vector<bool>
#
GET_CXX_SOURCE_RETURN_VALUE(
  "
  #include <vector>
  int main () {
    std::vector<bool> v(1);
    v.reserve (1);
    v.resize (1);
    return v.capacity();
  }
  "
  DEAL_II_MIN_BOOL_VECTOR_CAPACITY
  DEAL_II_MIN_BOOL_VECTOR_CAPACITY_RETURN_VALUE)

IF(NOT DEAL_II_MIN_BOOL_VECTOR_CAPACITY)
  # We have a problem...
  MESSAGE(WARNING
    "Could not determine DEAL_II_MIN_VECTOR_CAPACITY, "
    "source might not compile..."
    )
ELSE()
  SET(DEAL_II_MIN_BOOL_VECTOR_CAPACITY
    ${DEAL_II_MIN_BOOL_VECTOR_CAPACITY_RETURN_VALUE})
ENDIF()



#
# Newer versions of gcc can pass a flag to the assembler to
# compress debug sections. At the time of writing this test,
# this can save around 230 MB of disk space on the object
# files we produce (810MB down to 570MB for the debug versions
# of object files). Unfortunately, the sections have to be
# unpacked again when they are put into the shared libs, so
# no savings there.
#
# The flag also doesn't appear to be working on Cygwin, as
# per email by John Fowkes on the mailing list in Feb 2012,
# so don't run the test on cygwin.
#
IF(NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-Wa,--compress-debug-sections")
  ENABLE_IF_SUPPORTED(DEAL_II_C_FLAGS_DEBUG "-Wa,--compress-debug-sections")
ENDIF()

