## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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

#
# Check for various C++ language features
#
# This file sets up
#
#   DEAL_II_WITH_CXX14
#   DEAL_II_WITH_CXX17
#
#   DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE
#   DEAL_II_HAVE_ISNAN
#   DEAL_II_HAVE_UNDERSCORE_ISNAN
#   DEAL_II_HAVE_ISFINITE
#   DEAL_II_HAVE_FP_EXCEPTIONS
#   DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
#
#   DEAL_II_FALLTHROUGH
#


########################################################################
#                                                                      #
#                         C++ Version Support:                         #
#                                                                      #
########################################################################

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
    CHECK_CXX_COMPILER_FLAG("-std=c++${_suffix}" DEAL_II_HAVE_FLAG_stdcxx${_suffix})
    IF(DEAL_II_HAVE_FLAG_stdcxx${_suffix})
      SET(DEAL_II_CXX_VERSION_FLAG "-std=c++${_suffix}")
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
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
    PUSH_CMAKE_REQUIRED("-Werror")

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

    RESET_CMAKE_REQUIRED()
  ENDIF()

  IF( DEAL_II_HAVE_CXX17_ATTRIBUTES )
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
    # Set CMAKE_REQUIRED_FLAGS for the unit tests
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")

    #
    # This test does not guarantee full C++14 support, but virtually every
    # compiler with some C++14 support implements this.
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
    # Clang-3.5* or older, bail out with a spurious error message in case
    # of an undeduced auto return type.
    #
    # https://llvm.org/bugs/show_bug.cgi?id=16876
    #
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_FLAGS_DEBUG}")
    CHECK_CXX_SOURCE_COMPILES(
      "
      struct foo
      {
        auto func();
      };
      int main()
      {
        foo bar;
      }
      "
      DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK)

    RESET_CMAKE_REQUIRED()
  ENDIF()

  IF( DEAL_II_HAVE_CXX14_MAKE_UNIQUE AND
      DEAL_II_HAVE_CXX14_CLANGAUTODEBUG_BUG_OK )
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

PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
CHECK_CXX_SOURCE_COMPILES(
  "
  // common C++11 include files
  #include <array>
  #include <condition_variable>
  #include <type_traits>

  // type traits functionality
  constexpr auto m0 = std::is_trivial<double>::value;
  constexpr auto m1 = std::is_standard_layout<double>::value;
  constexpr auto m2 = std::is_pod<double>::value;

  // thread_local storage specification
  thread_local std::array<int,3> p;
  std::condition_variable c;

  // Check the version language macro, but skip MSVC because
  // MSVC reports 199711 even in MSVC 2017.
  #if __cplusplus < 201103L && !defined(_MSC_VER)
  #  error \"insufficient support for C++11\"
  #endif

  int main()
  {
    p[0];
    c.notify_all();
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
    char c = toupper('a');
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

#
# Enable c++17's [[fallthrough]] attribute.
#
IF(DEAL_II_WITH_CXX17)
  SET(DEAL_II_FALLTHROUGH "[[fallthrough]]")
ELSE()
  SET(DEAL_II_FALLTHROUGH " ")
ENDIF()


########################################################################
#                                                                      #
#                   Check for various C++ features:                    #
#                                                                      #
########################################################################

PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <type_traits>
  int main(){ std::is_trivially_copyable<int> bob; }
  "
  DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  int main(){ double d=0; std::isnan (d); return 0; }
  "
  DEAL_II_HAVE_STD_ISNAN)


CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  int main(){ double d=0; isnan (d); return 0; }
  "
  DEAL_II_HAVE_ISNAN)


CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  int main(){ double d=0; _isnan (d); return 0; }
  "
  DEAL_II_HAVE_UNDERSCORE_ISNAN)


CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  int main(){ double d=0; std::isfinite (d); return 0; }
  "
  DEAL_II_HAVE_ISFINITE)


#
# Check that we can use feenableexcept through the C++11 header file cfenv:
#
# The test is a bit more complicated because we also check that no garbage
# exception is thrown if we convert -std::numeric_limits<double>::max to a
# string. This sadly happens with some compiler support libraries :-(
#
# - Timo Heister, 2015
#

IF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
  CHECK_CXX_SOURCE_RUNS(
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
    DEAL_II_HAVE_FP_EXCEPTIONS)
ELSE()
  #
  # If we are not allowed to do platform introspection, just test whether
  # we can compile above code.
  #
  CHECK_CXX_SOURCE_COMPILES(
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
    DEAL_II_HAVE_FP_EXCEPTIONS)
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
RESET_CMAKE_REQUIRED()
