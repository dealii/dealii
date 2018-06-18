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
#                 Check for various compiler features:                 #
#                                                                      #
########################################################################

#
# This file sets up:
#
#   DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
#   DEAL_II_VECTOR_ITERATOR_IS_POINTER
#   DEAL_II_HAVE_BUILTIN_EXPECT
#   DEAL_II_HAVE_GLIBC_STACKTRACE
#   DEAL_II_HAVE_LIBSTDCXX_DEMANGLER
#   DEAL_II_COMPILER_HAS_ATTRIBUTE_PRETTY_FUNCTION
#   DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED
#   DEAL_II_COMPILER_HAS_ATTRIBUTE_ALWAYS_INLINE
#   DEAL_II_DEPRECATED
#   DEAL_II_ALWAYS_INLINE
#   DEAL_II_RESTRICT
#   DEAL_II_COMPILER_HAS_DIAGNOSTIC_PRAGMA
#   DEAL_II_COMPILER_HAS_FUSE_LD_GOLD
#

#
# A couple of test results depend on compiler flags and the C++ mode. Rerun
# these tests if necessary
#

UNSET_IF_CHANGED(CHECK_CXX_FEATURES_FLAGS_SAVED
  "${CMAKE_REQUIRED_FLAGS}${DEAL_II_CXX_VERSION_FLAG}${DEAL_II_WITH_CXX14}${DEAL_II_WITH_CXX17}"
  DEAL_II_COMPILER_HAS_CXX14_ATTRIBUTE_DEPRECATED
  DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED
  )


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
#ifdef __AVX512F__
  #include <immintrin.h>
#endif
  int main()
  {
    __m128d a, b;
    a = _mm_set_sd (1.0);
    b = _mm_set1_pd (2.1);
    __m128d c = a + b;
    __m128d d = b - c;
    __m128d e = c * a + d;
    __m128d f = e/a;
#ifdef __AVX512F__
    __m512d g, h;
    g = _mm512_set1_pd (1.0);
    h = _mm512_set1_pd (2.1);
    __m512d i = g + h;
    g = i - g;
    h *= i;
    i = h/i;
    (void)i;
#endif
    (void)f;
  }
  "
  DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS)


#
# Check whether the std::vector::iterator is just a plain pointer
#
# (Yes. It is not a bug. But the logic is the same.)
#
# - Matthias Maier, rewritten 2012
#
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
# GCC compilers. This can be used to hint the compiler's branch
# prediction unit in some cases. We use it in the AssertThrow
# macros.
#
# - Matthias Maier, rewritten 2012
#
CHECK_CXX_SOURCE_COMPILES(
  "
  bool f() {}
  int main(){ if (__builtin_expect(f(),false)) ; }
  "
  DEAL_II_HAVE_BUILTIN_EXPECT)


#
# Check whether glibc-like stacktrace information is available
# for the Exception class. If it is, then try to also determine
# whether the compiler accepts the -rdynamic flag, since that is
# recommended for linking if one wants to have meaningful
# backtraces.
#
# - Matthias Maier, rewritten 2012
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
  DEAL_II_HAVE_GLIBC_STACKTRACE)

IF(DEAL_II_HAVE_GLIBC_STACKTRACE AND NOT DEAL_II_STATIC_EXECUTABLE)
  ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-rdynamic")
ENDIF()


#
# Check whether the compiler offers a way to demangle symbols
# from within the program. Used inside the exception stacktrace
# mechanism.
#
# The example code is taken from
#   http://gcc.gnu.org/onlinedocs/libstdc++/18_support/howto.html#6
#
# - Matthias Maier, rewritten 2012
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
  DEAL_II_HAVE_LIBSTDCXX_DEMANGLER)


#
# GCC and some other compilers have __PRETTY_FUNCTION__, showing
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
# - Matthias Maier, rewritten 2012
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
    SET(__PRETTY_FUNCTION__ "__func__")
  ELSE()
    SET(__PRETTY_FUNCTION__ "\"(not available)\"")
  ENDIF()

ENDIF()


#
# Newer versions of GCC can pass a flag to the assembler to
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
# Finally, Intel's icpc compiler complains about the flag
# but apparently only if the file to be compiled contains
# particular content. See bug #46 in the Google Code bug
# data base (http://code.google.com/p/dealii/issues/detail?id=46).
# It proved impossible to track down under which circumstances
# this happens, and so it was disabled for icpc.
#
# - Matthias Maier, rewritten 2012, 2013
#
IF( (NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN") AND
    (NOT CMAKE_SYSTEM_NAME MATCHES "Windows") AND
    (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel") )
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-Wa,--compress-debug-sections")
ENDIF()


#
# GCC and some other compilers have an attribute of the form
# __attribute__((deprecated)) that can be used to make the
# compiler warn whenever a deprecated function is used. C++14
# provides a standardized attribute of the form [[deprecated]
# with the exact same functionality.
# See if one of these attribute is available.
#
# If it is, set the variable DEAL_II_DEPRECATED to its value. If
# it isn't, set it to an empty string (actually, to a single
# space, since the empty string causes CMAKE to #undef the
# variable in config.h), i.e., to something the compiler will
# ignore
#
# - Wolfgang Bangerth, 2012
#

# some compilers compile the attributes but they do not work:
# "warning: use of the 'deprecated' attribute is a C++14 extension" (clang in c++11 mode)
# "warning #1292: unknown attribute "deprecated"" (icc)
PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_FLAGS}")
PUSH_CMAKE_REQUIRED("-Werror")
PUSH_CMAKE_REQUIRED("-Wno-deprecated-declarations")
PUSH_CMAKE_REQUIRED("-Wno-deprecated")
PUSH_CMAKE_REQUIRED("-Wno-unused-command-line-argument")

# first see if the compiler accepts the attribute
CHECK_CXX_SOURCE_COMPILES(
  "
          [[deprecated]] int old_fn ();
          int old_fn () { return 0; }
          int (*fn_ptr)() = old_fn;

          struct [[deprecated]] bob
          {
            [[deprecated]] bob(int i);
            [[deprecated]] void test();
          };

          int main () {}
  "
  DEAL_II_COMPILER_HAS_CXX14_ATTRIBUTE_DEPRECATED
  )

CHECK_CXX_SOURCE_COMPILES(
  "
          __attribute__((deprecated)) int old_fn ();
          int old_fn () { return 0; }
          int (*fn_ptr)() = old_fn;

          struct __attribute__((deprecated)) bob
          {
            __attribute__((deprecated)) bob(int i);
            __attribute__((deprecated)) void test();
          };

          int main () {}
  "
  DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED
  )

RESET_CMAKE_REQUIRED()

IF(DEAL_II_COMPILER_HAS_CXX14_ATTRIBUTE_DEPRECATED)
  SET(DEAL_II_DEPRECATED "[[deprecated]]")
ELSEIF(DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED AND NOT DEAL_II_WITH_CUDA)
  SET(DEAL_II_DEPRECATED "__attribute__((deprecated))")
ELSE()
  SET(DEAL_II_DEPRECATED " ")
ENDIF()


#
# Do a similar check with the always_inline attribute on functions.
#
CHECK_CXX_SOURCE_COMPILES(
  "
          __attribute__((always_inline)) int fn () { return 0; }
          int main () { return fn(); }
  "
  DEAL_II_COMPILER_HAS_ATTRIBUTE_ALWAYS_INLINE
  )

IF(DEAL_II_COMPILER_HAS_ATTRIBUTE_ALWAYS_INLINE)
  SET(DEAL_II_ALWAYS_INLINE "__attribute__((always_inline))")
ELSE()
  SET(DEAL_II_ALWAYS_INLINE " ")
ENDIF()

#
# Check whether the compiler understands the __restrict keyword.
#
CHECK_CXX_SOURCE_COMPILES(
  "
          void fn (double *__restrict a, double *__restrict b) { a[0] = b[0]; a[1] = b[0]; }
          int main() { }
  "
  DEAL_II_COMPILER_HAS_RESTRICT_KEYWORD
  )

IF(DEAL_II_COMPILER_HAS_RESTRICT_KEYWORD)
  SET(DEAL_II_RESTRICT "__restrict")
ELSE()
  SET(DEAL_II_RESTRICT " ")
ENDIF()


#
# GCC and Clang allow fine grained control of diagnostics via the "GCC
# diagnostic" pragma. Check whether the compiler supports the "push" and
# "pop" mechanism and the "ignored" toggle. Further, test for the
# alternative "_Pragma(...)" variant (and that it does not emit a warning).
#
# - Matthias Maier, 2015
#
PUSH_CMAKE_REQUIRED("-Werror")
CHECK_CXX_SOURCE_COMPILES(
  "
  _Pragma(\"GCC diagnostic push\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wextra\\\\\\\"\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wunknown-pragmas\\\\\\\"\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wpragmas\\\\\\\"\")
  int main() { return 0; }
  _Pragma(\"GCC diagnostic pop\")
  "
  DEAL_II_COMPILER_HAS_DIAGNOSTIC_PRAGMA)
RESET_CMAKE_REQUIRED()


#
# Use the 'gold' linker if possible, given that it's substantially faster.
#
# We have to try to link a full executable with -fuse-ld=gold to check
# whether "ld.gold" is actually available.
#
# Clang always reports "argument unused during compilation"
# if "-fuse-ld=" is used, but fails at link time for an unsupported linker.
#
# ICC also emits a warning but passes for unsupported linkers
# unless we turn diagnostic warnings into errors.
#
# Wolfgang Bangerth, Matthias Maier, Daniel Arndt, 2015, 2018
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  PUSH_CMAKE_REQUIRED("-Wno-unused-command-line-argument")
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  PUSH_CMAKE_REQUIRED("-diag-error warn")
ENDIF()
PUSH_CMAKE_REQUIRED("-Werror")
PUSH_CMAKE_REQUIRED("-fuse-ld=gold")
CHECK_CXX_SOURCE_COMPILES(
  "
  int main() { return 0; }
  "
  DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)
RESET_CMAKE_REQUIRED()

IF(DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)
  ADD_FLAGS(DEAL_II_LINKER_FLAGS "-fuse-ld=gold")
ENDIF()

