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
#                 Check for various compiler features:                 #
#                                                                      #
########################################################################

#
# This file sets up:
#
#   DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
#   DEAL_II_VECTOR_ITERATOR_IS_POINTER
#   HAVE_BUILTIN_EXPECT
#   HAVE_VERBOSE_TERMINATE
#   HAVE_GLIBC_STACKTRACE
#   HAVE_LIBSTDCXX_DEMANGLER
#   DEAL_II_COMPILER_HAS_ATTRIBUTE_PRETTY_FUNCTION
#   DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED
#   DEAL_II_DEPRECATED
#


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
# gcc compilers. This can be used to hint the compiler's branch
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
# - Matthias Maier, rewritten 2012
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
  HAVE_GLIBC_STACKTRACE)

IF(HAVE_GLIBC_STACKTRACE AND NOT DEAL_II_STATIC_EXECUTABLE)
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-rdynamic")
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
  HAVE_LIBSTDCXX_DEMANGLER)


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
# Gcc and some other compilers have an attribute of the form
# __attribute__((deprecated)) that can be used to make the
# compiler warn whenever a deprecated function is used. See
# if this attribute is available.
#
# If it is, set the variable DEAL_II_DEPRECATED to its value. If
# it isn't, set it to an empty string (actually, to a single
# space, since the empty string causes CMAKE to #undef the
# variable in config.h), i.e., to something the compiler will
# ignore
#
# - Wolfgang Bangerth, 2012
#

# first see if the compiler accepts the attribute
CHECK_CXX_SOURCE_COMPILES(
  "
          int old_fn () __attribute__((deprecated));
          int old_fn () { return 0; }
          int (*fn_ptr)() = old_fn;

          int main () {}
  "
  DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED
  )

IF(DEAL_II_COMPILER_HAS_ATTRIBUTE_DEPRECATED)
  SET(DEAL_II_DEPRECATED "__attribute__((deprecated))")
ELSE()
  SET(DEAL_II_DEPRECATED " ")
ENDIF()

