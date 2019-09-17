## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II authors
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
#   DEAL_II_WITH_CXX14
#   DEAL_II_WITH_CXX17
#
#   DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH
#   DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE
#   DEAL_II_HAVE_CXX14_CONSTEXPR_CAN_CALL_NONCONSTEXPR
#   DEAL_II_HAVE_FP_EXCEPTIONS
#   DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
#
#   DEAL_II_CONSTEXPR
#   DEAL_II_FALLTHROUGH
#


#
# MSVC needs different compiler flags to turn warnings into errors
# additionally a suitable exception handling model is required
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  SET(_werror_flag "/WX /EHsc")
ELSE()
  SET(_werror_flag "-Werror")
ENDIF()


########################################################################
#                                                                      #
#                         C++ Version Support:                         #
#                                                                      #
########################################################################


#
#                              BIG DISCLAIMER
#
# Frankly speaking - this whole file is a mess and has to be rewritten and
# restructured at some point. So, if you are the lucky one who wants to add
# DEAL_II_WITH_CXX2X support, do not add yet another 300 lines of spaghetti
# code but restructure the whole thing to something sane.
#


#
IF(DEAL_II_WITH_CXX17 AND DEFINED DEAL_II_WITH_CXX14 AND NOT DEAL_II_WITH_CXX14)
  MESSAGE(FATAL_ERROR
    "Compiling deal.II with C++17 support (i.e., DEAL_II_WITH_CXX17=ON) requires"
    " that C++14 support not be explicitly disabled (i.e., DEAL_II_WITH_CXX14 may"
    " not be set to a logically false value)."
    )
ENDIF()

IF(DEFINED DEAL_II_WITH_CXX14 AND NOT DEAL_II_WITH_CXX14)
  SET(DEAL_II_WITH_CXX17 OFF CACHE STRING "" FORCE)
ENDIF()

#
# Check the user supplied DEAL_II_CXX_VERSION_FLAG
#

IF(NOT "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
  CHECK_CXX_COMPILER_FLAG(${DEAL_II_CXX_VERSION_FLAG} DEAL_II_CXX_VERSION_FLAG_VALID)
  IF(NOT DEAL_II_CXX_VERSION_FLAG_VALID)
    MESSAGE(FATAL_ERROR
      "The supplied flag \"${DEAL_II_CXX_VERSION_FLAG}\" was not recognized "
      "by the compiler."
      )
  ENDIF()

  #
  # A quick check that may, hopefully, give a more useful error message
  #
  IF("${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "-std=c++98" OR
     "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "-std=c++03")
    MESSAGE(FATAL_ERROR
      "\ndeal.II no longer supports compilation under the C++98 or C++03 "
      "standards: please either specify a value for DEAL_II_CXX_VERSION_FLAG "
      "corresponding to C++11 or leave the field blank (so that CMake may "
      "automatically determine a valid version flag).\n"
      )
  ENDIF()

  SET(_user_provided_cxx_version_flag TRUE)
ELSE()
  SET(_user_provided_cxx_version_flag FALSE)
ENDIF()

#
# A macro to check for various C++ version flags
#

MACRO(_check_cxx_flag _suffix)
  IF("${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
    IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
      CHECK_CXX_COMPILER_FLAG("/std:c++${_suffix}" DEAL_II_HAVE_FLAG_stdcxx${_suffix})
      IF(DEAL_II_HAVE_FLAG_stdcxx${_suffix})
        SET(DEAL_II_CXX_VERSION_FLAG "/std:c++${_suffix}")
      ENDIF()
    ELSE()
      CHECK_CXX_COMPILER_FLAG("-std=c++${_suffix}" DEAL_II_HAVE_FLAG_stdcxx${_suffix})
      IF(DEAL_II_HAVE_FLAG_stdcxx${_suffix})
        SET(DEAL_II_CXX_VERSION_FLAG "-std=c++${_suffix}")
      ENDIF()
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO(_check_version _version _symbolic)
  _check_cxx_flag("${_version}")
  _check_cxx_flag("${_symbolic}")
  IF(DEAL_II_WITH_CXX${_version} AND "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
    MESSAGE(FATAL_ERROR
      "C++${_version} support was requested but CMake was not able to find a valid C++${_version}"
      " flag. Try to manually specify DEAL_II_CXX_VERSION_FLAG and rerun CMake."
      )
  ENDIF()
ENDMACRO()

#
# Check for proper C++17 support and set up DEAL_II_HAVE_CXX17:
#
IF(NOT DEFINED DEAL_II_WITH_CXX17 OR DEAL_II_WITH_CXX17)
  _check_version("17" "1z")

  IF(NOT "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
    # Set CMAKE_REQUIRED_FLAGS for the unit tests
    MESSAGE(STATUS "Using C++ version flag \"${DEAL_II_CXX_VERSION_FLAG}\"")
    ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG} ${_werror_flag}")

    UNSET_IF_CHANGED(CHECK_CXX_FEATURES_FLAGS_CXX17_SAVED
      "${CMAKE_REQUIRED_FLAGS}${DEAL_II_CXX_VERSION_FLAG}"
      DEAL_II_HAVE_CXX17_ATTRIBUTES
      DEAL_II_HAVE_CXX17_IF_CONSTEXPR
      )

    #
    # Test that the c++17 attributes are supported.
    #
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <iostream>

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
      }
      "
      DEAL_II_HAVE_CXX17_ATTRIBUTES)

    #
    # Test that the c++17 if constexpr is supported.
    #
    CHECK_CXX_SOURCE_COMPILES(
      "
      int main()
      {
        constexpr bool flag = false;
        if constexpr(flag)
          return 1;
        return 0;
      }
      "
      DEAL_II_HAVE_CXX17_IF_CONSTEXPR)

    #
    # Some compilers treat lambdas as constexpr functions when compiling
    # with C++17 support even if they don't fulfill all the constexpr
    # function requirements. Consequently, these compilers don't allow
    # try-blocks or non-literal return types in lambdas. This is a bug.
    #
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
      DEAL_II_NON_CONSTEXPR_LAMBDA)

    RESET_CMAKE_REQUIRED()
  ENDIF()

  IF( DEAL_II_HAVE_CXX17_ATTRIBUTES AND
      DEAL_II_HAVE_CXX17_IF_CONSTEXPR AND
      DEAL_II_NON_CONSTEXPR_LAMBDA)
    SET(DEAL_II_HAVE_CXX17 TRUE)
  ELSE()
    IF(NOT _user_provided_cxx_version_flag)
      SET(DEAL_II_CXX_VERSION_FLAG "")
    ENDIF()
  ENDIF()
ENDIF()



#
# Check for proper C++14 support and set up DEAL_II_HAVE_CXX14:
#
IF(NOT DEFINED DEAL_II_WITH_CXX14 OR DEAL_II_WITH_CXX14)
  _check_version("14" "1y")

  IF(NOT "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")

    UNSET_IF_CHANGED(CHECK_CXX_FEATURES_FLAGS_CXX14_SAVED
      "${CMAKE_REQUIRED_FLAGS}${DEAL_II_CXX_VERSION_FLAG}"
      DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK
      DEAL_II_HAVE_CXX14_CONSTEXPR_STDMAXMIN
      DEAL_II_HAVE_CXX14_MAKE_UNIQUE
      )

    # Set CMAKE_REQUIRED_FLAGS for the unit tests
    ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

    #
    # We assume std::make_unique works
    #
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <memory>

      int main()
      {
        auto ptr = std::make_unique<int>(42);
        return 0;
      }
      "
      DEAL_II_HAVE_CXX14_MAKE_UNIQUE)

    #
    # This test checks constexpr std::max/min support. Unfortunately,
    # gcc-4.9 claims to support C++14 but fails to provide a constexpr
    # compatible std::max/min. Disable C++14 support in this case.
    #
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <algorithm>
      int main()
      {
          constexpr int max = std::max(0,1);
          (void) max;
      }
      "
      DEAL_II_HAVE_CXX14_CONSTEXPR_STDMAXMIN)


    #
    # Clang-3.5* or older, bail out with a spurious error message in case
    # of an undeduced auto return type.
    #
    # https://llvm.org/bugs/show_bug.cgi?id=16876
    #
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

    RESET_CMAKE_REQUIRED()
  ENDIF()

  IF( DEAL_II_HAVE_CXX14_MAKE_UNIQUE AND
      DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK AND
      DEAL_II_HAVE_CXX14_CONSTEXPR_STDMAXMIN)
    SET(DEAL_II_HAVE_CXX14 TRUE)
  ELSE()
    IF(NOT _user_provided_cxx_version_flag)
      SET(DEAL_II_CXX_VERSION_FLAG "")
    ENDIF()
  ENDIF()
ENDIF()

#
# Check for proper C++11 support.
#

#
# Set up a default C++11 flag in case both C++14 detection failed and the user
# did not specify a flag:
#
IF("${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
  _check_version("11" "0x")
ENDIF()

UNSET_IF_CHANGED(CHECK_CXX_FEATURES_FLAGS_CXX11_SAVED
  "${CMAKE_REQUIRED_FLAGS}${DEAL_II_CXX_VERSION_FLAG}"
  DEAL_II_HAVE_CXX11_FEATURES
  DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK
  DEAL_II_HAVE_CXX11_ICCLIBSTDCPP47CXX11BUG_OK
  DEAL_II_HAVE_CXX11_ICCNUMERICLIMITSBUG_OK
  DEAL_II_HAVE_CXX11_MACOSXC99BUG_OK
  )

ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
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
    p[0];
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

#
# On Mac OS-X 10.9 with recent gcc compilers in C++11 mode linking to
# some standard C library functions, notably toupper and tolower, fail
# due to unresolved references to these functions.
#
# Thanks to Denis Davydov for the testcase.
#
# Matthias Maier, 2013
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <ctype.h>
  int main ()
  {
    int c = toupper('a');
    (void) c;
  }
  "
  DEAL_II_HAVE_CXX11_MACOSXC99BUG_OK)


#
# icc-13 triggers an internal compiler error when compiling
# std::numeric_limits<...>::min() with -std=c++0x [1].
#
# Reported by Ted Kord.
#
# - Matthias Maier, 2013
#
# [1] http://software.intel.com/en-us/forums/topic/328902
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <limits>
  struct Integer
  {
    static const int min_int_value;
    static const int max_int_value;
  };
  const int Integer::min_int_value = std::numeric_limits<int>::min();
  const int Integer::max_int_value = std::numeric_limits<int>::max();
  int main() { return 0; }
  "
  DEAL_II_HAVE_CXX11_ICCNUMERICLIMITSBUG_OK)


#
# icc-14.0.0 has an astonishing bug [1] where it hits an internal compiler
# error when run in C++11 mode with libstdc++-4.7 (from gcc).
#
# [1] http://software.intel.com/en-us/forums/topic/472385
#
# - Matthias Maier, 2013
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <vector>
  template<typename T> void foo()
  {
    std::vector<double> data(100);
  }
  int main()
  {
    foo<int>();
  }
  "
  DEAL_II_HAVE_CXX11_ICCLIBSTDCPP47CXX11BUG_OK)
RESET_CMAKE_REQUIRED()

IF( DEAL_II_HAVE_CXX11_FEATURES AND
    DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK AND
    DEAL_II_HAVE_CXX11_MACOSXC99BUG_OK AND
    DEAL_II_HAVE_CXX11_ICCNUMERICLIMITSBUG_OK AND
    DEAL_II_HAVE_CXX11_ICCLIBSTDCPP47CXX11BUG_OK )
  # we require C++11 so this is fine
ELSE()
  MESSAGE(FATAL_ERROR "\n"
    "The current combination of compiler and C++ version flag is missing some "
    "features of the C++11 version of the language necessary for compiling "
    "deal.II. Please ensure that the CMake variable DEAL_II_CXX_VERSION_FLAG is "
    "set to the correct flag to enable C++11 support; if it is already set and "
    "you still see this message then you will need to upgrade your compiler.\n")
ENDIF()


#
# Set up a configuration option for C++14 and C++17 support:
#

IF (DEAL_II_HAVE_CXX17 AND NOT DEAL_II_HAVE_CXX14)
  MESSAGE(STATUS "Disabling CXX17 support because CXX14 detection failed.")
  SET(DEAL_II_HAVE_CXX17 FALSE)
ENDIF()

OPTION(DEAL_II_WITH_CXX14
  "Compile deal.II using C++14 language standard."
  ${DEAL_II_HAVE_CXX14}
  )
LIST(APPEND DEAL_II_FEATURES CXX14)
SET(FEATURE_CXX14_PROCESSED TRUE)

OPTION(DEAL_II_WITH_CXX17
  "Compile deal.II using C++17 language standard."
  ${DEAL_II_HAVE_CXX17}
  )
LIST(APPEND DEAL_II_FEATURES CXX17)
SET(FEATURE_CXX17_PROCESSED TRUE)

#
# Bail out if user requested support for a certain C++ version (e.g.,
# DEAL_II_WITH_CXX14) but support is not available due to above tests
#

MACRO(_bailout _version _symbolic)
  IF(DEAL_II_WITH_CXX${_version} AND NOT DEAL_II_HAVE_CXX${_version})
    MESSAGE(FATAL_ERROR "\n"
      "C++${_version} support was requested (DEAL_II_WITH_CXX${_version}=${DEAL_II_WITH_CXX${_version}}) "
      "but it is not supported by the current compiler.\n"
      "Please disable C++${_version} support, i.e. configure with\n"
      "    -DDEAL_II_WITH_CXX${_version}=FALSE,\n"
      "or use a different compiler, instead. (If the compiler flag for C++${_version} "
      "support differs from \"-std=c++${_symbolic}\" or \"-std=c++${_version}\", "
      "a suitable compiler flag has to be specified manually via\n"
      "    -DDEAL_II_CXX_VERSION_FLAG=\"...\"\n\n"
      )
  ENDIF()
ENDMACRO()

_bailout("14" "1y")
_bailout("17" "1z")

#
# Some compilers (notably GCC6 and up) default to C++14 (more exactly, GNU++14,
# which contains some extensions) but most compilers default to C++98. Hence,
# try to avoid adding an extra flag by doing one last test:
#
RESET_CMAKE_REQUIRED()
IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "/Zc:__cplusplus")
ENDIF()
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <memory>

  #if __cplusplus < 201103L
  #  error \"The compiler does not default to C++11 or newer.\"
  #endif

  auto main() -> int
  {
    auto p0 = std::unique_ptr<int>();
    auto p1 = std::move(p0);
  }
  "
  DEAL_II_COMPILER_DEFAULTS_TO_CXX11_OR_NEWER)
RESET_CMAKE_REQUIRED()

IF(_user_provided_cxx_version_flag OR
    NOT DEAL_II_COMPILER_DEFAULTS_TO_CXX11_OR_NEWER OR
    (DEAL_II_COMPILER_DEFAULTS_TO_CXX11_OR_NEWER AND NOT DEAL_II_WITH_CXX14))
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "Using C++ version flag \"${DEAL_II_CXX_VERSION_FLAG}\"")
ENDIF ()

IF (DEAL_II_WITH_CXX17)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX14 successfully set up")
  MESSAGE(STATUS "DEAL_II_WITH_CXX17 successfully set up")
ELSEIF (DEAL_II_WITH_CXX14)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX14 successfully set up")
ELSE()
  MESSAGE(STATUS
    "DEAL_II_WITH_CXX17 and DEAL_II_WITH_CXX14 are both disabled")
ENDIF()


########################################################################
#                                                                      #
#                   Check for various C++ features:                    #
#                                                                      #
########################################################################

UNSET_IF_CHANGED(CHECK_CXX_FEATURES_FLAGS_SAVED
  "${CMAKE_REQUIRED_FLAGS}${DEAL_II_CXX_VERSION_FLAG}${DEAL_II_WITH_CXX14}${DEAL_II_WITH_CXX17}"
  DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH
  DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE
  DEAL_II_HAVE_CXX14_CONSTEXPR_CAN_CALL_NONCONSTEXPR
  DEAL_II_HAVE_FP_EXCEPTIONS
  DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
  )

#
# Try to enable a fallthrough attribute. This is a language feature in C++17,
# but a compiler extension in earlier language versions: check both
# possibilities here.
#
ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS}")
ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${_werror_flag}")
IF(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "-Wno-unused-command-line-argument")
ENDIF()
#
# first try the attribute [[fallthrough]]
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

RESET_CMAKE_REQUIRED()
IF(DEAL_II_HAVE_CXX17_ATTRIBUTE_FALLTHROUGH)
  SET(DEAL_II_FALLTHROUGH "[[fallthrough]]")
ELSEIF(DEAL_II_HAVE_ATTRIBUTE_FALLTHROUGH)
  SET(DEAL_II_FALLTHROUGH "__attribute__((fallthrough))")
ELSE()
  SET(DEAL_II_FALLTHROUGH " ")
ENDIF()

ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <type_traits>
  int main(){ std::is_trivially_copyable<int> bob; }
  "
  DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE)

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
# As long as there exists an argument value such that an invocation of the
# function or constructor could be an evaluated subexpression of a core constant
# expression, C++14 allows to call non-constexpr functions from constexpr
# functions. Unfortunately, not all compilers obey the standard in this regard.
#
CHECK_CXX_SOURCE_COMPILES(
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
  DEAL_II_HAVE_CXX14_CONSTEXPR_CAN_CALL_NONCONSTEXPR)

#
# The macro DEAL_II_CONSTEXPR allows using c++ constexpr features in a portable way.
# Here we enable it only when a constexpr function can call simple non-constexpr
# functions. This requirement is probabely very conservative in most cases, but
# it will prevent breaking builds with certain compilers.
#
IF (DEAL_II_HAVE_CXX14_CONSTEXPR_CAN_CALL_NONCONSTEXPR)
  SET(DEAL_II_CONSTEXPR "constexpr")
ELSE()
  SET(DEAL_II_CONSTEXPR " ")
ENDIF()

RESET_CMAKE_REQUIRED()
