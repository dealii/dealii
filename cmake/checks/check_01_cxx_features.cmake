## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
#   DEAL_II_HAVE_ISNAN
#   DEAL_II_HAVE_UNDERSCORE_ISNAN
#   DEAL_II_HAVE_ISFINITE
#


########################################################################
#                                                                      #
#                         C++ Version Support:                         #
#                                                                      #
########################################################################

#
# backwards compatability with the old DEAL_II_CXX11_FLAG option
#

IF(DEFINED DEAL_II_CXX11_FLAG AND NOT DEFINED DEAL_II_CXX_VERSION_FLAG)
  SET(DEAL_II_CXX_VERSION_FLAG ${DEAL_II_CXX11_FLAG})
ENDIF()

IF(DEAL_II_WITH_CXX11 AND DEAL_II_WITH_CXX14)
  MESSAGE(STATUS
          "Since DEAL_II_WITH_CXX11 and DEAL_II_WITH_CXX14 were both specified,\n"
          "DEAL_II_WITH_CXX14 is chosen as the standard.")
ENDIF()

IF(DEAL_II_WITH_CXX14 AND DEFINED DEAL_II_WITH_CXX11 AND NOT DEAL_II_WITH_CXX11)
  MESSAGE(FATAL_ERROR
    "Compiling deal.II with C++14 support (i.e., DEAL_II_WITH_CXX14=ON) requires"
    " that C++11 support not be explicitly disabled (i.e., DEAL_II_WITH_CXX11 may"
    " not be set to a logically false value).")
ENDIF()

#
# Only run these tests if C++11 support (or newer) should actually be set up:
#

IF(DEAL_II_WITH_CXX11 OR DEAL_II_WITH_CXX14)
  IF(NOT DEAL_II_CXX_VERSION_FLAG)
    IF(DEAL_II_WITH_CXX14)
      CHECK_CXX_COMPILER_FLAG("-std=c++14" DEAL_II_HAVE_FLAG_stdcxx14)
      IF(DEAL_II_HAVE_FLAG_stdcxx14)
        SET(DEAL_II_CXX_VERSION_FLAG "-std=c++14")
      ELSE()
        MESSAGE(FATAL_ERROR
          "CMake was not able to find a valid C++14 flag. Try manually"
          " specifying\nDEAL_II_CXX_VERSION_FLAG\nand reruning CMake.")
      ENDIF()
    ELSEIF(DEAL_II_WITH_CXX11)
      CHECK_CXX_COMPILER_FLAG("-std=c++0x" DEAL_II_HAVE_FLAG_stdcxx0x)
      CHECK_CXX_COMPILER_FLAG("-std=c++11" DEAL_II_HAVE_FLAG_stdcxx11)
      IF(DEAL_II_HAVE_FLAG_stdcxx11)
        SET(DEAL_II_CXX_VERSION_FLAG "-std=c++11")
      ELSEIF(DEAL_II_HAVE_FLAG_stdcxx0x)
        SET(DEAL_II_CXX_VERSION_FLAG "-std=c++0x")
      ELSE()
        MESSAGE(FATAL_ERROR
          "CMake was not able to find a valid C++11 flag. Try manually"
          " specifying\nDEAL_II_CXX_VERSION_FLAG and reruning CMake.")
      ENDIF()
    ENDIF()
  ELSE()
    CHECK_CXX_COMPILER_FLAG(${DEAL_II_CXX_VERSION_FLAG} DEAL_II_CXX_VERSION_FLAG_VALID)
    IF(NOT DEAL_II_CXX_VERSION_FLAG_VALID)
      MESSAGE(FATAL_ERROR
        "The supplied flag \"${DEAL_II_CXX_VERSION_FLAG}\" was not recognized "
        "by the detected compiler.")
    ENDIF()
  ENDIF()

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

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <memory>
    std::shared_ptr<int> p(new int(3));
    int main(){ return 0; }
    "
    DEAL_II_HAVE_CXX11_SHARED_PTR)

  #
  # On some systems with gcc 4.5.0, we can compile the code
  # below but it will throw an exception when run. So test
  # that as well.
  #
  IF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
    PUSH_CMAKE_REQUIRED("-pthread")
    CHECK_CXX_SOURCE_RUNS(
      "
      #include <thread>
      void f(int){}
      int main(){ std::thread t(f,1); t.join(); return 0; }
      "
      DEAL_II_HAVE_CXX11_THREAD)
    RESET_CMAKE_REQUIRED()
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
  ELSE()
    # Just export it ;-)
    SET(DEAL_II_HAVE_CXX11_THREAD TRUE CACHE BOOL "")
  ENDIF()

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


  IF( DEAL_II_HAVE_CXX11_ARRAY AND
      DEAL_II_HAVE_CXX11_CONDITION_VARIABLE AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL_GCCBUG35569_OK AND
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
  IF( DEAL_II_HAVE_CXX11 AND
      DEAL_II_HAVE_CXX14_MAKE_UNIQUE AND
      DEAL_II_WITH_CXX14)
    SET(DEAL_II_HAVE_CXX14 TRUE)
  ENDIF()

  RESET_CMAKE_REQUIRED()
ELSE()
  SET(DEAL_II_WITH_CXX11 FALSE)
  SET(DEAL_II_WITH_CXX14 FALSE)
ENDIF()

#
# Set up a configuration option for C++11 support:
#

OPTION(DEAL_II_WITH_CXX11
  "Compile deal.II using C++11 language standard."
  ${DEAL_II_HAVE_CXX11}
  )

#
# Set up a configuration option for C++14 support:
#

OPTION(DEAL_II_WITH_CXX14
  "Compile deal.II using C++14 language standard."
  ${DEAL_II_HAVE_CXX14}
  )

#
# Bail out if user requested C++11 support (DEAL_II_WITH_CXX11) but support
# is not available due to above tests (DEAL_II_HAVE_CXX11):
#

IF(DEAL_II_WITH_CXX11 AND NOT DEAL_II_HAVE_CXX11)
  MESSAGE(FATAL_ERROR "\n"
    "C++11 support was requested (DEAL_II_WITH_CXX11=TRUE) but is not "
    "supported by the current compiler.\n"
    "Please disable C++11 support, i.e. configure with\n"
    "    -DDEAL_II_WITH_CXX11=FALSE,\n"
    "or use a different compiler, instead. (If the compiler flag for C++11 "
    "support differs from \"-std=c++0x\" or \"-std=c++11\", a suitable "
    "compiler flag has to be specified manually via\n"
    "    -DDEAL_II_CXX_VERSION_FLAG=\"...\"\n\n"
    )
ENDIF()

#
# similar for C++14
#

IF(DEAL_II_WITH_CXX14 AND NOT DEAL_II_HAVE_CXX14)
  MESSAGE(FATAL_ERROR "\n"
    "C++14 support was requested (DEAL_II_WITH_CXX14=${DEAL_II_WITH_CXX14}) "
    "but is not supported by the current compiler and flag combination.\n"
    "Please disable C++14 support, i.e. configure with\n"
    "    -DDEAL_II_WITH_CXX14=FALSE,\n"
    "or use a different compiler, instead. (If the compiler flag for C++14 "
    "support differs from \"-std=c++14\", a suitable "
    "compiler flag has to be specified manually via\n"
    "    -DDEAL_II_CXX_VERSION_FLAG=\"...\"\n\n"
    )
ENDIF()

#
# Set up support for newer versions of C++:
#

IF (DEAL_II_WITH_CXX14)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX14 successfully set up")
ELSEIF(DEAL_II_WITH_CXX11)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX11 successfully set up")

  PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <type_traits>
    int main(){ std::is_trivially_copyable<int> bob; }
    "
    DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE)
  RESET_CMAKE_REQUIRED()
ELSE()
  MESSAGE(STATUS "DEAL_II_WITH_CXX14 and DEAL_II_WITH_CXX11 are both disabled")
ENDIF()


########################################################################
#                                                                      #
#                   Check for various C++ features:                    #
#                                                                      #
########################################################################

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

