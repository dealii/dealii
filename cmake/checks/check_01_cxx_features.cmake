## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Check for various C++ language features
#
# This file sets up
#
#   DEAL_II_HAVE_CXX17
#   DEAL_II_HAVE_CXX20
#   DEAL_II_HAVE_CXX23
#
#   DEAL_II_HAVE_FP_EXCEPTIONS
#   DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
#   DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
#   DEAL_II_HAVE_CXX17_LEGENDRE_FUNCTIONS
#   DEAL_II_FALLTHROUGH
#   DEAL_II_CONSTEXPR
#


########################################################################
#                                                                      #
#                         C++ Version Support:                         #
#                                                                      #
########################################################################

#
# We need compiler flags specified in ${DEAL_II_CXX_FLAGS} for all the
# tests. Create a small macro to easily set CMAKE_REQUIRED_FLAGS
#
macro(_set_up_cmake_required)
  reset_cmake_required()
  set(CMAKE_REQUIRED_FLAGS "")
  add_flags(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_WARNING_FLAGS}")
  add_flags(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS}")
  add_flags(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS_SAVED}")
endmacro()


#
# Wrap the following checks into a macro to make it easier to rerun them.
#
macro(_test_cxx23_support)
  unset_if_changed(CHECK_CXX23_FEATURES_FLAGS_SAVED
    "${CMAKE_REQUIRED_FLAGS}${CMAKE_CXX_STANDARD}"
    DEAL_II_HAVE_CXX23_FEATURES
    )

  # Strictly speaking "202100L" indicates support for a preliminary version
  # of the C++23 standard (which will have "202302L" when finalized). gcc-13
  # exports this version number when configured with C++23 support.
  add_flags(CMAKE_REQUIRED_FLAGS "-Werror")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <version>

    #if __cplusplus < 202100L && !defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #  error \"insufficient support for C++23\"
    #endif

    struct Dummy {
      static void operator()() {}
    };

    int main() {
      [[assume(2 > 1)]];
    }
    "
    DEAL_II_HAVE_CXX23_FEATURES)
  reset_cmake_required()

  if(DEAL_II_HAVE_CXX23_FEATURES)
    message(STATUS "C++23 support is enabled.")
    set(DEAL_II_HAVE_CXX23 TRUE)
    set(_cxx_standard 23)
  else()
    message(STATUS "C++23 support is disabled.")
    set(DEAL_II_HAVE_CXX23 FALSE)
  endif()
endmacro()


#
# Wrap the following checks into a macro to make it easier to rerun them.
#
macro(_test_cxx20_support)
  unset_if_changed(CHECK_CXX20_FEATURES_FLAGS_SAVED
    "${CMAKE_REQUIRED_FLAGS}${CMAKE_CXX_STANDARD}"
    DEAL_II_HAVE_CXX20_FEATURES
    )

  # Strictly speaking "201709L" indicates support for a preliminary version
  # of the C++20 standard (which will have "202002L" when finalized). gcc-10
  # exports this version number when configured with C++20 support.
  # clang-10 exports the final "202002L" version instead, as does gcc-11.
  #
  # Beyond this, check for some features we actually need.
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <version>

    #if __cplusplus < 201709L && !defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #  error \"insufficient support for C++20\"
    #endif

    #if !(defined __cpp_lib_ranges) || (__cpp_lib_ranges < 201911)
    #  error \"insufficient support for C++20\"
    #endif

    #if !(defined __cpp_lib_type_identity)
    #  error \"insufficient support for C++20: __cpp_lib_type_identity not defined\"
    #endif

    #if !(defined __cpp_lib_type_identity) || (__cpp_lib_type_identity < 201806)
    #  error \"insufficient support for C++20: __cpp_lib_type_identity is too old \"
    #endif

    // Test concepts and requires clauses
    template <int dim, int spacedim>
    concept is_valid_dim_spacedim = (dim >= 1 && spacedim <= 3 &&
                                     dim <= spacedim);

    template <int dim, int spacedim>
    requires is_valid_dim_spacedim<dim,spacedim>
    class Triangulation
    {};

    Triangulation<1,3> t;


    int main()
    {
    }
    "
    DEAL_II_HAVE_CXX20_FEATURES)

  if(DEAL_II_HAVE_CXX20_FEATURES)
    message(STATUS "C++20 support is enabled.")
    set(DEAL_II_HAVE_CXX20 TRUE)
    set(_cxx_standard 20)
  else()
    message(STATUS "C++20 support is disabled.")
    set(DEAL_II_HAVE_CXX20 FALSE)
  endif()
endmacro()


#
# Wrap the following checks into a macro to make it easier to rerun them.
#
macro(_test_cxx17_support)
  unset_if_changed(CHECK_CXX17_FEATURES_FLAGS_SAVED
    "${CMAKE_REQUIRED_FLAGS}${CMAKE_CXX_STANDARD}"
    DEAL_II_HAVE_CXX17_FEATURES
    DEAL_II_HAVE_CXX17_CONSTEXPR_LAMBDA_BUG_OK
    DEAL_II_HAVE_CXX14_FEATURES
    DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK
    DEAL_II_HAVE_CXX11_FEATURES
    DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK
    )

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <iostream>
    #include <optional>
    #include <tuple>

    #if __cplusplus < 201703L && !defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #  error \"insufficient support for C++17\"
    #endif

    //check for some C++17 features that we use in our headers:
    using std::apply;
    using std::optional;

    [[nodiscard]] int test_nodiscard()
    {
      return 1;
    }

    int main()
    {
      const unsigned int n=1;
      switch (n)
      {
        case 1:
          std::cout << n;
          [[fallthrough]];
        case 2:
          std::cout << n;
      }

      [[maybe_unused]] int i = test_nodiscard();

      constexpr bool flag = false;
      if constexpr(flag)
        return 1;
      return 0;
    }
    "
    DEAL_II_HAVE_CXX17_FEATURES)

  # Some compilers treat lambdas as constexpr functions when compiling with
  # C++17 support even if they don't fulfill all the constexpr function
  # requirements. Consequently, these compilers don't allow try-blocks or
  # non-literal return types in lambdas. This is a bug.
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <string>
    int main()
    {
      auto c = []()
      {
        return std::string{};
      }();
      (void) c;

      return []()
      {
        try
        {}
        catch(...)
        {}
        return 0;
      }();
    }
    "
    DEAL_II_HAVE_CXX17_CONSTEXPR_LAMBDA_BUG_OK)

  # Check some generic C++14 features
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <memory>
    #include <algorithm>

    // Check the version language macro, but skip MSVC because
    // MSVC reports 199711 even in MSVC 2017.
    #if __cplusplus < 201402L && !defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #  error \"insufficient support for C++14\"
    #endif

    int main()
    {
      auto ptr = std::make_unique<int>(42);
      constexpr int max = std::max(0, 1);
      (void) ptr;
      (void) max;
      return 0;
    }
    "
    DEAL_II_HAVE_CXX14_FEATURES)

  # Clang-3.5* or older, bail out with a spurious error message in case
  # of an undeduced auto return type.
  #
  # https://llvm.org/bugs/show_bug.cgi?id=16876
  set(_flags "${DEAL_II_CXX_FLAGS_DEBUG}")
  strip_flag(_flags "-Wa,--compress-debug-sections")
  add_flags(CMAKE_REQUIRED_FLAGS "${_flags}")
  CHECK_CXX_SOURCE_COMPILES(
    "
    struct foo
    {
      auto func();
    };
    int main()
    {
      foo bar;
      (void) bar;
    }
    "
    DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK)
  _set_up_cmake_required()

  # Check some generic C++11 features
  CHECK_CXX_SOURCE_COMPILES(
    "
    // common C++11 include files
    #include <array>
    #include <condition_variable>
    #include <type_traits>

    // thread_local storage specification
    static thread_local std::array<int,3> p;

    // Check the version language macro, but skip MSVC because
    // MSVC reports 199711 even in MSVC 2017.
    #if __cplusplus < 201103L && !defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #  error \"insufficient support for C++11\"
    #endif

    int main()
    {
      std::condition_variable c;
      p[0] = 1;
      c.notify_all();

     // type traits functionality
     constexpr auto m0 = std::is_trivial<double>::value;
     (void) m0;
     constexpr auto m1 = std::is_standard_layout<double>::value;
     (void) m1;
    }
    "
    DEAL_II_HAVE_CXX11_FEATURES)

  # clang libc++ bug, see https://llvm.org/bugs/show_bug.cgi?id=20084
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <functional>
    struct A { void foo() const {} };
    int main() { A a; std::bind(&A::foo,a)(); return 0; }
    "
    DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK)

  if(DEAL_II_HAVE_CXX17_FEATURES AND
     DEAL_II_HAVE_CXX17_CONSTEXPR_LAMBDA_BUG_OK AND
     DEAL_II_HAVE_CXX14_FEATURES AND
     DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK AND
     DEAL_II_HAVE_CXX11_FEATURES AND
     DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK)
    message(STATUS "C++17 support is enabled.")
    set(DEAL_II_HAVE_CXX17 TRUE)
    set(_cxx_standard 17)
  else()
    message(STATUS "C++17 support is disabled.")
    set(DEAL_II_HAVE_CXX17 FALSE)
  endif()
endmacro()


#
# Try to find out what we support:
#

_set_up_cmake_required()
_test_cxx17_support()

if(NOT DEAL_II_HAVE_CXX17)
  #
  # We failed to detect C++17 support. Let's make an attempt to set the
  # -std= compiler flag. (But in order to minimize confusion let's not
  # override any manually specified -std= flag or CMAKE_CXX_STANDARD
  # variable set by the user.)
  #
  if(NOT "${DEAL_II_CXX_FLAGS_SAVED}" MATCHES "-std=" AND "${CMAKE_CXX_STANDARD}" STREQUAL "")
    message(STATUS "C++17 support not available. Try to set -std=c++17 explicitly")
    set(CMAKE_CXX_STANDARD 17) # manually set the C++ standard
    _set_up_cmake_required()
    _test_cxx17_support()
  endif()
endif()

if(NOT DEAL_II_HAVE_CXX17)
  message(FATAL_ERROR
    "\nThe current version of deal.II requires a compiler with enabled "
    "C++17 support. Make sure to use a modern enough compiler (GCC version "
    "9 onwards, Clang version 10 onwards, or Microsoft MS VS 2019 onwards) "
    "and check that the compiler flag \"-std=\" is either unset, or set to "
    "at least c++17. Similarly, please make sure that the CMake variable "
    "CMAKE_CXX_STANDARD is either unset, or set at least to 17.\n\n"
    )
endif()

_test_cxx20_support()

_test_cxx23_support()

set_if_empty(CMAKE_CXX_STANDARD "${_cxx_standard}")
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT "${CMAKE_CXX_STANDARD}" STREQUAL "${_cxx_standard}")
  message(FATAL_ERROR
    "\nThe current version of deal.II was configured with CMAKE_CXX_STANDARD "
    "set to \"${CMAKE_CXX_STANDARD}\", but we detected only support for standard "
    "version \"${_cxx_standard}\". Either unset the CMake variable "
    "CMAKE_CXX_STANDARD, or ensure that it is at most set to \"${_cxx_standard}\".\n\n"
    )
endif()


########################################################################
#                                                                      #
#                   Check for various C++ features:                    #
#                                                                      #
########################################################################


#
# Some compilers are too generous in accepting some of the language
# features that we test below and do not issue an error but a warning. Set
# -pedantic -Werror to make the feature detection more reliable.
#
set(_werror_flag "")
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  enable_if_supported(_werror_flag "/WX /EHsc")
else()
  enable_if_supported(_werror_flag "-pedantic")
  enable_if_supported(_werror_flag "-Werror")
  enable_if_supported(_werror_flag "-Wno-unused-command-line-argument")
endif()

#
# We need to reset CMAKE_REQUIRED_* variables again after the above call to
# enable_if_supported()
#
_set_up_cmake_required()
add_flags(CMAKE_REQUIRED_FLAGS "${_werror_flag}")

unset_if_changed(CHECK_CXX_FEATURES_FLAGS_SAVED
  "${CMAKE_REQUIRED_FLAGS}${CMAKE_CXX_STANDARD}"
  DEAL_II_HAVE_FP_EXCEPTIONS
  DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
  DEAL_II_HAVE_CXX17_ATTRIBUTE_FALLTHROUGH
  DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH
  DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
  DEAL_II_HAVE_CXX17_LEGENDRE_FUNCTIONS
  DEAL_II_CXX14_CONSTEXPR_BUG_OK
  )


#
# Check that we can use feenableexcept through the C++11 header file cfenv:
#
# The test is a bit more complicated because we also check that no garbage
# exception is thrown if we convert std::numeric_limits<double>::lowest to a
# string. This sadly happens with some compiler support libraries :-(
#
# - Timo Heister, 2015
#
set(_snippet
  "
  #include <cfenv>
  #include <limits>
  #include <sstream>

  int main()
  {
    feenableexcept(FE_DIVBYZERO|FE_INVALID);
    std::ostringstream description;
    const double lower_bound = std::numeric_limits<double>::lowest();

    description << lower_bound;

    return 0;
  }
  "
  )
if(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
  CHECK_CXX_SOURCE_RUNS("${_snippet}" DEAL_II_HAVE_FP_EXCEPTIONS)
else()
  #
  # If we are not allowed to do platform introspection, just test whether
  # we can compile above code.
  #
  CHECK_CXX_SOURCE_COMPILES("${_snippet}" DEAL_II_HAVE_FP_EXCEPTIONS)
endif()


#
# Check whether the standard library provides operator* overloads for
# mixed floating point multiplication of complex and real valued
# numbers. The C++ standard does not say that such overloads should
# exist, and so they shouldn't. In other words, a standards compliant
# compiler will fail this check.
#
# In practice, the absence of these overloads makes it quite difficult
# to write mixed-precision linear algebra functions. To make this less
# of a pain, we declare such overloads ourselves (in namespace dealii)
# unless the compiler already has them (which, as mentioned above, it
# shouldn't).
#
# - Matthias Maier, 2015
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <complex>

  int main()
  {
    double() * std::complex<float>();
    std::complex<float>() * double();
    float() * std::complex<double>();
    std::complex<double>() * float();
    std::complex<double>() * std::complex<float>();
    std::complex<float>() * std::complex<double>();

    return 0;
  }
  "
  DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS)

#
# Try to enable a fallthrough attribute. This is a language feature in C++17,
# but a compiler extension in earlier language versions.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  int main()
  {
    int i = 42;
    int j = 10;
    switch(i)
      {
      case 1:
        ++j;
        [[fallthrough]];
      case 2:
        ++j;
        [[fallthrough]];
      default:
        break;
      }
   }
   "
   DEAL_II_HAVE_CXX17_ATTRIBUTE_FALLTHROUGH
   )

#
# see if the current compiler configuration supports the GCC extension
# __attribute__((fallthrough)) syntax instead
#
CHECK_CXX_SOURCE_COMPILES(
  "
  int main()
  {
    int i = 42;
    int j = 10;
    switch(i)
      {
      case 1:
        ++j;
        __attribute__((fallthrough));
      case 2:
        ++j;
        __attribute__((fallthrough));
      default:
        break;
      }
  }
  "
  DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH
  )

if(DEAL_II_HAVE_CXX17_ATTRIBUTE_FALLTHROUGH)
  set(DEAL_II_FALLTHROUGH "[[fallthrough]]")
elseif(DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH)
  set(DEAL_II_FALLTHROUGH "__attribute__((fallthrough))")
else()
  set(DEAL_II_FALLTHROUGH " ")
endif()


#
# Check for c++17 Bessel function support. Unfortunately libc++ version 10
# does not have those.
#

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  using std::cyl_bessel_j;
  using std::cyl_bessel_jf;
  int main()
  {
  }
  "
  DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
  )


#
# Check for c++17 Legendre function support.
#

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  using std::legendre;
  using std::legendref;
  int main()
  {
  }
  "
  DEAL_II_HAVE_CXX17_LEGENDRE_FUNCTIONS
  )


#
# Check for correct c++14 constexpr support.
#
# As long as there exists an argument value such that an invocation of the
# function or constructor could be an evaluated subexpression of a core constant
# expression, C++14 allows to call non-constexpr functions from constexpr
# functions.
#
# Unfortunately, not all compilers obey the standard in this regard. In some
# cases, MSVC 2019 crashes with an internal compiler error when we
# declare the respective functions as 'constexpr' even though the test below
# passes, see #9080.
#
# We only run this check if we have CXX14 support, otherwise the use of constexpr
# is limited (non-const constexpr functions for example).
#

# MSVC has considerable problems with "constexpr", disable unconditionally
# for now
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(DEAL_II_CXX14_CONSTEXPR_BUG true)
else()
  check_cxx_compiler_bug(
    "
    #define Assert(x,y) if (!(x)) throw y;
    void bar()
    {}

    constexpr int
    foo(const int n)
    {
      Assert(n>0, \"hello\");
      if(!(n >= 0))
        bar();
      return n;
    }

    int main()
    {
      constexpr unsigned int n=foo(1);
      return n;
    }
    "
    DEAL_II_CXX14_CONSTEXPR_BUG)
endif()

set(DEAL_II_CONSTEXPR "constexpr")
if(DEAL_II_CXX14_CONSTEXPR_BUG)
  set(DEAL_II_CONSTEXPR " ")
endif()


reset_cmake_required()
