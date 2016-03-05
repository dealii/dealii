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
#   DEAL_II_WITH_CXX11
#   DEAL_II_WITH_CXX14
#
#   DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE
#   DEAL_II_HAVE_ISNAN
#   DEAL_II_HAVE_UNDERSCORE_ISNAN
#   DEAL_II_HAVE_ISFINITE
#   DEAL_II_HAVE_FP_EXCEPTIONS
#   DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
#


########################################################################
#                                                                      #
#                         C++ Version Support:                         #
#                                                                      #
########################################################################

#
# backwards compatibility with the old DEAL_II_CXX11_FLAG option
#
SET_IF_EMPTY(DEAL_II_CXX_VERSION_FLAG "${DEAL_II_CXX11_FLAG}")

IF(DEAL_II_WITH_CXX14 AND DEFINED DEAL_II_WITH_CXX11 AND NOT DEAL_II_WITH_CXX11)
  MESSAGE(FATAL_ERROR
    "Compiling deal.II with C++14 support (i.e., DEAL_II_WITH_CXX14=ON) requires"
    " that C++11 support not be explicitly disabled (i.e., DEAL_II_WITH_CXX11 may"
    " not be set to a logically false value)."
    )
ENDIF()

IF(DEFINED DEAL_II_WITH_CXX11 AND NOT DEAL_II_WITH_CXX11)
  SET(DEAL_II_WITH_CXX14 OFF CACHE STRING "" FORCE)
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

  SET(_user_provided_cxx_version_flag TRUE)
ENDIF()

#
# A macro to check for various C++11 and C++14 flags
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
# Check for proper C++14 support and set up DEAL_II_HAVE_CXX14:
#
IF(NOT DEFINED DEAL_II_WITH_CXX14 OR DEAL_II_WITH_CXX14)
  _check_version("14" "1y")

  IF(NOT "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
    # Set CMAKE_REQUIRED_FLAGS for the unit tests
    MESSAGE(STATUS "Using C++ version flag \"${DEAL_II_CXX_VERSION_FLAG}\"")
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
# Check for proper C++11 support and set up DEAL_II_HAVE_CXX11:
#
IF(NOT DEFINED DEAL_II_WITH_CXX11 OR DEAL_II_WITH_CXX11)

  IF("${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
    _check_version("11" "0x")
  ENDIF()

  IF(NOT "${DEAL_II_CXX_VERSION_FLAG}" STREQUAL "")
    # Set CMAKE_REQUIRED_FLAGS for the unit tests
    MESSAGE(STATUS "Using C++ version flag \"${DEAL_II_CXX_VERSION_FLAG}\"")
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <array>
      std::array<int,3> p;
      int main(){  p[0]; return 0; }
      "
      DEAL_II_HAVE_CXX11_ARRAY)

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <condition_variable>
      std::condition_variable c;
      int main(){ c.notify_all(); return 0; }
      "
      DEAL_II_HAVE_CXX11_CONDITION_VARIABLE)

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <functional>
      void f(int, double){}
      std::function<void (int)> g = std::bind (f, std::placeholders::_1,1.1);
      int main(){ return 0; }
      "
      DEAL_II_HAVE_CXX11_FUNCTIONAL)

    # Make sure we don't run into GCC bug 35569
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <functional>
      void f(int){}
      using namespace std;
      using namespace std::placeholders;
      int main(){ bind(multiplies<int>(),4,_1)(5); return 0; }
      "
      DEAL_II_HAVE_CXX11_FUNCTIONAL_GCCBUG35569_OK)

    # clang libc++ bug, see https://llvm.org/bugs/show_bug.cgi?id=20084
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <functional>
      struct A { void foo() const {} };
      int main() { A a; std::bind(&A::foo,a)(); return 0; }
      "
      DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK)

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <memory>
      std::shared_ptr<int> p(new int(3));
      int main(){ return 0; }
      "
      DEAL_II_HAVE_CXX11_SHARED_PTR)

    PUSH_CMAKE_REQUIRED("-pthread")
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <thread>
      void f(int){}
      int main(){ std::thread t(f,1); t.join(); return 0; }
      "
      DEAL_II_HAVE_CXX11_THREAD)
    RESET_CMAKE_REQUIRED()
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <mutex>
      std::mutex m;
      int main(){ m.lock(); return 0; }
      "
      DEAL_II_HAVE_CXX11_MUTEX)

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <tuple>],
      std::tuple<int,double,char> p(1,1.1,'a');
      int main(){ return 0; }
      "
      DEAL_II_HAVE_CXX11_TUPLE)

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <type_traits>
      const bool m0 = std::is_trivial<double>::value;
      const bool m1 = std::is_standard_layout<double>::value;
      const bool m2 = std::is_pod<double>::value;
      int main(){ return 0; }
      "
      DEAL_II_HAVE_CXX11_TYPE_TRAITS)

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
    # We just disable C++11 mode in this case
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
  ENDIF()

  IF( DEAL_II_HAVE_CXX11_ARRAY AND
      DEAL_II_HAVE_CXX11_CONDITION_VARIABLE AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL_GCCBUG35569_OK AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL_LLVMBUG20084_OK AND
      DEAL_II_HAVE_CXX11_SHARED_PTR AND
      DEAL_II_HAVE_CXX11_THREAD AND
      DEAL_II_HAVE_CXX11_MUTEX AND
      DEAL_II_HAVE_CXX11_TUPLE AND
      DEAL_II_HAVE_CXX11_TYPE_TRAITS AND
      DEAL_II_HAVE_CXX11_MACOSXC99BUG_OK AND
      DEAL_II_HAVE_CXX11_ICCNUMERICLIMITSBUG_OK AND
      DEAL_II_HAVE_CXX11_ICCLIBSTDCPP47CXX11BUG_OK )
      SET(DEAL_II_HAVE_CXX11 TRUE)
  ENDIF()
ENDIF()


#
# Finally disable cxx14 if cxx11 detection failed for whatever reason. This
# can happen if any of our compile checks above fails, for example threading
# support.
# 
IF (DEAL_II_HAVE_CXX14 AND NOT DEAL_II_HAVE_CXX11)
  MESSAGE(STATUS "Disabling CXX14 support because CXX11 detection failed.")
  SET(DEAL_II_HAVE_CXX14 FALSE)
ENDIF()

#
# Set up a configuration options for C++11 and C++14 support:
#

OPTION(DEAL_II_WITH_CXX11
  "Compile deal.II using C++11 language standard."
  ${DEAL_II_HAVE_CXX11}
  )
LIST(APPEND DEAL_II_FEATURES CXX11)

OPTION(DEAL_II_WITH_CXX14
  "Compile deal.II using C++14 language standard."
  ${DEAL_II_HAVE_CXX14}
  )
LIST(APPEND DEAL_II_FEATURES CXX14)

#
# Bail out if user requested C++11 support (DEAL_II_WITH_CXX11) but support
# is not available due to above tests (DEAL_II_HAVE_CXX11):
#

MACRO(_bailout _version)
  IF(DEAL_II_WITH_CXX${_version} AND NOT DEAL_II_HAVE_CXX${_version})
    MESSAGE(FATAL_ERROR "\n"
      "C++${_version} support was requested (DEAL_II_WITH_CXX${_version}=${DEAL_II_WITH_CXX${_version}}) but is not "
      "supported by the current compiler.\n"
      "Please disable C++${_version} support, i.e. configure with\n"
      "    -DDEAL_II_WITH_CXX${_version}=FALSE,\n"
      "or use a different compiler, instead. (If the compiler flag for C++${_version} "
      "support differs from \"-std=c++0x\" or \"-std=c++${_version}\", a suitable "
      "compiler flag has to be specified manually via\n"
      "    -DDEAL_II_CXX_VERSION_FLAG=\"...\"\n\n"
      )
  ENDIF()
ENDMACRO()

_bailout("11")
_bailout("14")

IF (DEAL_II_WITH_CXX14)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX11 successfully set up")
  MESSAGE(STATUS "DEAL_II_WITH_CXX14 successfully set up")
ELSEIF(DEAL_II_WITH_CXX11)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX11 successfully set up")
ELSE()
  MESSAGE(STATUS "DEAL_II_WITH_CXX14 and DEAL_II_WITH_CXX11 are both disabled")
ENDIF()


########################################################################
#                                                                      #
#                   Check for various C++ features:                    #
#                                                                      #
########################################################################

#
# Some compilers (such as Intel 15.3 and GCC 4.9.2) support the flags
# "-std=c++11" and "-std=c++14" but do not support
# 'std::is_trivially_copyable', so check for support in C++11 or newer.
#
IF(DEAL_II_WITH_CXX11)
  PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
  CHECK_CXX_SOURCE_COMPILES(
    "
  #include <type_traits>
  int main(){ std::is_trivially_copyable<int> bob; }
  "
    DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE)
  RESET_CMAKE_REQUIRED()
ELSE()
  SET(DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE FALSE)
ENDIF()

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

# This test requires C++11
IF(DEAL_II_WITH_CXX11)
  PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
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
  RESET_CMAKE_REQUIRED()
ELSE()
  SET(DEAL_II_HAVE_FP_EXCEPTIONS FALSE)
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

