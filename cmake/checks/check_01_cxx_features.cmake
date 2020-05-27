## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
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

#
# Check for various C++ language features
#
# This file sets up
#
#   DEAL_II_HAVE_CXX14
#   DEAL_II_HAVE_CXX17
#   DEAL_II_HAVE_CXX20
#
#   DEAL_II_HAVE_FP_EXCEPTIONS
#   DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
#   DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
#   DEAL_II_FALLTHROUGH
#   DEAL_II_DEPRECATED
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
MACRO(_set_up_cmake_required)
  # Let's put the user supplied `DEAL_II_CXX_FLAGS_SAVED` last so that we
  # never override a user supplied -std=c++XY flag in our tests.
  RESET_CMAKE_REQUIRED()
  SET(CMAKE_REQUIRED_FLAGS "")
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS}")
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS_SAVED}")
ENDMACRO()


#
# Wrap the following checks into a macro to make it easier to rerun them.
#
MACRO(_test_cxx20_support)

  UNSET_IF_CHANGED(CHECK_CXX20_FEATURES_FLAGS_SAVED
    "${CMAKE_REQUIRED_FLAGS}"
    DEAL_II_HAVE_CXX20_FEATURES
    )

  # Strictly speaking "201709L" indicates support for a preliminary version
  # of C++20 standard (which will have "202002L" when finalized). gcc-10
  # exports this version number when configured with C++20 support.
  # clang-10 exports the final "202002L" version instead.
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <cmath>
    #include <ranges>

    #if __cplusplus < 201709L && !defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #  error \"insufficient support for C++20\"
    #endif

    #if !(defined __cpp_lib_ranges) || (__cpp_lib_ranges < 201911)
    #  error \"insufficient support for C++20\"
    #endif

    int main()
    {
    }
    "
    DEAL_II_HAVE_CXX20_FEATURES)

  IF(DEAL_II_HAVE_CXX20_FEATURES)
    MESSAGE(STATUS "C++20 support is enabled.")
    SET(DEAL_II_HAVE_CXX20 TRUE)
  ELSE()
    MESSAGE(STATUS "C++20 support is disabled.")
    SET(DEAL_II_HAVE_CXX20 FALSE)
  ENDIF()
ENDMACRO()


#
# Wrap the following checks into a macro to make it easier to rerun them.
#
MACRO(_test_cxx17_support)

  UNSET_IF_CHANGED(CHECK_CXX17_FEATURES_FLAGS_SAVED
    "${CMAKE_REQUIRED_FLAGS}"
    DEAL_II_HAVE_CXX17_FEATURES
    DEAL_II_HAVE_CXX17_CONSTEXPR_LAMBDA_BUG_OK
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

  IF(DEAL_II_HAVE_CXX17_FEATURES AND
     DEAL_II_HAVE_CXX17_CONSTEXPR_LAMBDA_BUG_OK)
    MESSAGE(STATUS "C++17 support is enabled.")
    SET(DEAL_II_HAVE_CXX17 TRUE)
  ELSE()
    MESSAGE(STATUS "C++17 support is disabled.")
    SET(DEAL_II_HAVE_CXX17 FALSE)
  ENDIF()
ENDMACRO()


#
# Wrap the following checks into a macro to make it easier to rerun them.
#
MACRO(_test_cxx14_support)
  UNSET_IF_CHANGED(CHECK_CXX14_FEATURES_FLAGS_SAVED
    "${CMAKE_REQUIRED_FLAGS}"
    DEAL_II_HAVE_CXX14_FEATURES
    DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK
    DEAL_II_HAVE_CXX11_FEATURES
    DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK
    )

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
  SET(_flags "${DEAL_II_CXX_FLAGS_DEBUG}")
  STRIP_FLAG(_flags "-Wa,--compress-debug-sections")
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${_flags}")
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
     constexpr auto m2 = std::is_pod<double>::value;
     (void) m2;
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

  IF(DEAL_II_HAVE_CXX14_FEATURES AND
     DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK AND
     DEAL_II_HAVE_CXX11_FEATURES AND
     DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK)
    MESSAGE(STATUS "C++14 support is enabled.")
    SET(DEAL_II_HAVE_CXX14 TRUE)
  ELSE()
    MESSAGE(STATUS "C++14 support is disabled.")
    SET(DEAL_II_HAVE_CXX14 FALSE)
  ENDIF()
ENDMACRO()


#
# Try to find out what we support:
#

_set_up_cmake_required()
_test_cxx14_support()

IF(NOT DEAL_II_HAVE_CXX14)
  MESSAGE(STATUS "C++14 support not available. Try to set -std=c++14 explicitly")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-std=c++14")
  _set_up_cmake_required()
  _test_cxx14_support()
ENDIF()

IF(NOT DEAL_II_HAVE_CXX14)
  MESSAGE(FATAL_ERROR
    "\nThe current version of deal.II requires a compiler with enabled "
    "C++14 support. Make sure to use a modern enough compiler (GCC version "
    "5 onwards, Clang version 4 onwards, or Microsoft MS VS 2015 onwards) "
    "and check that the compiler flag \"-std=\" is either unset, or set to "
    "at least c++14.\n\n"
    )
ENDIF()

_test_cxx17_support()
_test_cxx20_support()


########################################################################
#                                                                      #
#                   Check for various C++ features:                    #
#                                                                      #
########################################################################


#
# Some compilers are too generous in accepting some of the language
# features that we test below and do not issue an error but a warning. Set
# -Werror to make the feature detection more reliable.
#
SET(_werror_flag "")
IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  ENABLE_IF_SUPPORTED(_werror_flag "/WX /EHsc")
ELSE()
  ENABLE_IF_SUPPORTED(_werror_flag "-Werror")
  ENABLE_IF_SUPPORTED(_werror_flag "-Wno-unused-command-line-argument")
ENDIF()
ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${_werror_flag}")

UNSET_IF_CHANGED(CHECK_CXX_FEATURES_FLAGS_SAVED
  "${CMAKE_REQUIRED_FLAGS}"
  DEAL_II_HAVE_FP_EXCEPTIONS
  DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
  DEAL_II_HAVE_CXX17_ATTRIBUTE_DEPRECATED
  DEAL_II_HAVE_ATTRIBUTE_DEPRECATED
  DEAL_II_HAVE_CXX17_ATTRIBUTE_FALLTHROUGH
  DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH
  DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
  DEAL_II_CXX14_CONSTEXPR_BUG_OK
  )


#
# Check that we can use feenableexcept through the C++11 header file cfenv:
#
# The test is a bit more complicated because we also check that no garbage
# exception is thrown if we convert -std::numeric_limits<double>::max to a
# string. This sadly happens with some compiler support libraries :-(
#
# - Timo Heister, 2015
#
SET(_snippet
  "
  #include <cfenv>
  #include <limits>
  #include <sstream>

  int main()
  {
    feenableexcept(FE_DIVBYZERO|FE_INVALID);
    std::ostringstream description;
    const double lower_bound = -std::numeric_limits<double>::max();

    description << lower_bound;

    return 0;
  }
  "
  )
IF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
  CHECK_CXX_SOURCE_RUNS("${_snippet}" DEAL_II_HAVE_FP_EXCEPTIONS)
ELSE()
  #
  # If we are not allowed to do platform introspection, just test whether
  # we can compile above code.
  #
  CHECK_CXX_SOURCE_COMPILES("${_snippet}" DEAL_II_HAVE_FP_EXCEPTIONS)
ENDIF()


#
# Check whether the standard library provides operator* overloads for mixed
# floating point multiplication of complex and real valued numbers.
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
# Even though [[deprecated]] is a C++14 feature we have to check
# wether we can actually use the [[deprecated]] attribute in all
# cases we care about; some of the following are C++17 features.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  [[deprecated]] int old_fn ();
  int old_fn () { return 0; }

  struct [[deprecated]] bob
  {
    [[deprecated]] bob(int i);
    [[deprecated]] void test();
  };

  enum color
  {
    red [[deprecated]]
  };

  template <int dim>
  struct foo {};
  using bar [[deprecated]] = foo<2>;

  int main () {}
  "
  DEAL_II_HAVE_CXX17_ATTRIBUTE_DEPRECATED
  )

#
# Also test the corresponding GCC extension
#
CHECK_CXX_SOURCE_COMPILES(
  "
  __attribute__((deprecated)) int old_fn ();
  int old_fn () { return 0; }

  struct __attribute__((deprecated)) bob
  {
    __attribute__((deprecated)) bob(int i);
    __attribute__((deprecated)) void test();
  };

  enum color
  {
    red __attribute__((deprecated))
  };

  template <int dim>
  struct foo {};
  using bar __attribute__((deprecated)) = foo<2>;

  int main () {}
  "
  DEAL_II_HAVE_ATTRIBUTE_DEPRECATED
  )

IF(DEAL_II_HAVE_CXX17_ATTRIBUTE_DEPRECATED)
  SET(DEAL_II_DEPRECATED "[[deprecated]]")
ELSEIF(DEAL_II_HAVE_ATTRIBUTE_DEPRECATED AND NOT DEAL_II_WITH_CUDA)
  SET(DEAL_II_DEPRECATED "__attribute__((deprecated))")
ELSE()
  SET(DEAL_II_DEPRECATED " ")
ENDIF()


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

IF(DEAL_II_HAVE_CXX17_ATTRIBUTE_FALLTHROUGH)
  SET(DEAL_II_FALLTHROUGH "[[fallthrough]]")
ELSEIF(DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH)
  SET(DEAL_II_FALLTHROUGH "__attribute__((fallthrough))")
ELSE()
  SET(DEAL_II_FALLTHROUGH " ")
ENDIF()


#
# Check for c++17 bessel function support. Unfortunately libc++ version 10
# does not have those.
#

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  using std::cyl_bessel_j;
  using std::cyl_bessel_jf;
  using std::cyl_bessel_jl;
  int main()
  {
  }
  "
  DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
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
IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  SET(DEAL_II_CXX14_CONSTEXPR_BUG true)
ELSE()
  CHECK_CXX_COMPILER_BUG(
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
ENDIF()

SET(DEAL_II_CONSTEXPR "constexpr")
IF(DEAL_II_CXX14_CONSTEXPR_BUG)
  SET(DEAL_II_CONSTEXPR " ")
ENDIF()


RESET_CMAKE_REQUIRED()
