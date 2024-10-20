## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
#   DEAL_II_COMPILER_HAS_ATTRIBUTE_ALWAYS_INLINE
#   DEAL_II_ALWAYS_INLINE
#   DEAL_II_RESTRICT
#   DEAL_II_COMPILER_HAS_DIAGNOSTIC_PRAGMA
#   DEAL_II_COMPILER_HAS_FUSE_LD_GOLD
#

#
# A couple of test results depend on compiler flags and the C++ mode.
# Nota Bene: If your test depends on the value of compile flags set in
# ${DEAL_II_CXX_FLAGS} it is probably a language feature and should go into
# check_01_cxx_features.cmake
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
  #include <x86intrin.h>
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
check_cxx_compiler_bug(
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
# Intel compilers don't handle __builtin_expect in C++14 constexpr contexts
# properly so we disable this feature in case we are going to use
# DEAL_II_CONSTEXPR with an Intel compiler.
#
# - Matthias Maier, rewritten 2012
#
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  CHECK_CXX_SOURCE_COMPILES(
    "
    bool f() { return true; }
    int main(){ if (__builtin_expect(f(),false)) {} }
    "
    DEAL_II_HAVE_BUILTIN_EXPECT)
endif()


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

if(DEAL_II_HAVE_GLIBC_STACKTRACE)
  enable_if_links(DEAL_II_LINKER_FLAGS "-rdynamic")
endif()


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

if(NOT DEAL_II_COMPILER_HAS_ATTRIBUTE_PRETTY_FUNCTION)
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

  if(DEAL_II_COMPILER_HAS_ATTRIBUTE_FUNC)
    set(__PRETTY_FUNCTION__ "__func__")
  else()
    set(__PRETTY_FUNCTION__ "\"(not available)\"")
  endif()
endif()


#
# Newer versions of GCC can pass a flag to the assembler to compress
# debug sections. At the time of writing this test, this can save
# around 230 MB of disk space on the object files we produce (810MB
# down to 570MB for the debug versions of object files). Regardless of
# the setting, the linker will store the debug information
# uncompressed by default. Fortunately, we can also enable compressed
# output for the linker.
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
if( (NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN") AND
    (NOT CMAKE_SYSTEM_NAME MATCHES "Windows") AND
    (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel") )
  add_flags(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS_DEBUG}")
  enable_if_supported(DEAL_II_CXX_FLAGS_DEBUG "-Wa,--compress-debug-sections")
  enable_if_links(DEAL_II_LINKER_FLAGS_DEBUG "-Wl,--compress-debug-sections=zlib")
  reset_cmake_required()
endif()


#
# Do a similar check with the always_inline attribute on functions.
#
CHECK_CXX_SOURCE_COMPILES(
  "
          __attribute__((always_inline)) inline int fn () { return 0; }
          int main () { return fn(); }
  "
  DEAL_II_COMPILER_HAS_ATTRIBUTE_ALWAYS_INLINE
  )

if(DEAL_II_COMPILER_HAS_ATTRIBUTE_ALWAYS_INLINE)
  set(DEAL_II_ALWAYS_INLINE "__attribute__((always_inline))")
else()
  set(DEAL_II_ALWAYS_INLINE " ")
endif()


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

if(DEAL_II_COMPILER_HAS_RESTRICT_KEYWORD)
  set(DEAL_II_RESTRICT "__restrict")
else()
  set(DEAL_II_RESTRICT " ")
endif()


#
# GCC and Clang allow fine grained control of diagnostics via the "GCC
# diagnostic" pragma. Check whether the compiler supports the "push" and
# "pop" mechanism and the "ignored" toggle. Further, test for the
# alternative "_Pragma(...)" variant (and that it does not emit a warning).
#
# - Matthias Maier, 2015
#
add_flags(CMAKE_REQUIRED_FLAGS "${_werror_flag}")
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
reset_cmake_required()


#
# Use 'mold', 'lld' or the 'gold' linker if possible, given that either of them
# is substantially faster.
#
# We have to try to link a full executable with -fuse-ld=mold, -fuse-ld=lld or
# -fuse-ld=gold to check whether "ld.mold", "ld.lld" or "ld.gold" is actually
# available.
#
# Clang always reports "argument unused during compilation", but fails at link
# time for an unsupported linker.
#
# ICC also emits a warning but passes for unsupported linkers
# unless we turn diagnostic warnings into errors.
#
# We also test linker support with "-shared -fPIC". This catches an
# incompatibility where LLD refuses to produce a shared object from an
# object file compiled by the Intel Compiler:
#
#   ld.lld: error: can't create dynamic relocation R_X86_64_64 against symbol:
#   __gxx_personality_v0 in readonly segment; recompile object files with -fPIC
#   or pass '-Wl,-z,notext' to allow text relocations in the output
#
# even if we actually had -fPIC option present. If we add -Wl,-z,notext, it
# will link, but the produced libdeal_II.so is faulty and will crash randomly.
#
# Wolfgang Bangerth, Matthias Maier, Daniel Arndt, Binrui Dong, 2015, 2018-2020
#

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_flags(CMAKE_REQUIRED_FLAGS "-Wno-unused-command-line-argument")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_flags(CMAKE_REQUIRED_FLAGS "-diag-error warn")
  endif()
  add_flags(CMAKE_REQUIRED_FLAGS "-Werror")
  add_flags(CMAKE_REQUIRED_FLAGS "-shared")
  add_flags(CMAKE_REQUIRED_FLAGS "-fPIC")

  #
  # Check for ld.mold, ld.lld and ld.gold support:
  #
  add_flags(CMAKE_REQUIRED_FLAGS "-fuse-ld=mold")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <iostream>
    void foo() { std::cout << \"Hello, world!\" << std::endl; }
    "
    DEAL_II_COMPILER_HAS_FUSE_LD_MOLD)

  strip_flag(CMAKE_REQUIRED_FLAGS "-fuse-ld=mold")
  add_flags(CMAKE_REQUIRED_FLAGS "-fuse-ld=lld")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <iostream>
    void foo() { std::cout << \"Hello, world!\" << std::endl; }
    "
    DEAL_II_COMPILER_HAS_FUSE_LD_LLD)

  strip_flag(CMAKE_REQUIRED_FLAGS "-fuse-ld=lld")
  add_flags(CMAKE_REQUIRED_FLAGS "-fuse-ld=gold")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <iostream>
    void foo() { std::cout << \"Hello, world!\" << std::endl; }
    "
    DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)

  if(DEAL_II_COMPILER_HAS_FUSE_LD_MOLD)
    add_flags(DEAL_II_LINKER_FLAGS "-fuse-ld=mold")
  elseif(DEAL_II_COMPILER_HAS_FUSE_LD_LLD)
    add_flags(DEAL_II_LINKER_FLAGS "-fuse-ld=lld")
  elseif(DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)
    add_flags(DEAL_II_LINKER_FLAGS "-fuse-ld=gold")
  endif()

  reset_cmake_required()
endif()


#
# Check whether the compiler supports interprocedural and link-time
# optimization. If so, we can set the corresponding property on the
# compile and link targets later on.
#
include(CheckIPOSupported)
check_ipo_supported(RESULT _res OUTPUT _out LANGUAGES CXX)
if (_res)
  message(STATUS "Compiler supports interprocedural/link-time optimizations")
  set(DEAL_II_COMPILER_SUPPORTS_IPO "YES")
else()
  message(STATUS "Compiler does not support interprocedural/link-time optimizations: ${_out}")
  set(DEAL_II_COMPILER_SUPPORTS_IPO "NO")
endif()
