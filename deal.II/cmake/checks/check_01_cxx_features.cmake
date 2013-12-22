## ---------------------------------------------------------------------
## $Id$
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

#
# Check for various C++ language features
#
# This file sets up
#
#   DEAL_II_WITH_CXX11
#
#   HAVE_ISNAN
#   HAVE_UNDERSCORE_ISNAN
#   DEAL_II_HAVE_ISFINITE
#


########################################################################
#                                                                      #
#                             C++11 Support:                           #
#                                                                      #
########################################################################


#
# Only run these tests if C++11 support should actually be set up:
#
IF(NOT DEFINED DEAL_II_WITH_CXX11 OR DEAL_II_WITH_CXX11)

  IF("${DEAL_II_CXX11_FLAG}" STREQUAL "")
    CHECK_CXX_COMPILER_FLAG("-std=c++11" DEAL_II_HAVE_FLAG_stdcxx11)
    IF(DEAL_II_HAVE_FLAG_stdcxx11)
      SET(DEAL_II_CXX11_FLAG "-std=c++11")
    ELSEIF(DEAL_II_HAVE_FLAG_stdcxx11)
      CHECK_CXX_COMPILER_FLAG("-std=c++0x" DEAL_II_HAVE_FLAG_stdcxx0x)
      SET(DEAL_II_CXX11_FLAG "-std=x++0x")
    ENDIF()
  ENDIF()

  # Set CMAKE_REQUIRED_FLAGS for the unit tests
  MESSAGE(STATUS "Using C++11 flag \"${DEAL_II_CXX11_FLAG}\"")
  PUSH_TEST_FLAG("${DEAL_II_CXX11_FLAG}")

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
    PUSH_TEST_FLAG("-pthread")
    CHECK_CXX_SOURCE_RUNS(
      "
      #include <thread>
      void f(int){}
      int main(){ std::thread t(f,1); t.join(); return 0; }
      "
      DEAL_II_HAVE_CXX11_THREAD)
    POP_TEST_FLAG()
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

  IF( DEAL_II_HAVE_CXX11_ARRAY AND
      DEAL_II_HAVE_CXX11_CONDITION_VARIABLE AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL AND
      DEAL_II_HAVE_CXX11_FUNCTIONAL_GCCBUG35569_OK AND
      DEAL_II_HAVE_CXX11_SHARED_PTR AND
      DEAL_II_HAVE_CXX11_THREAD AND
      DEAL_II_HAVE_CXX11_MUTEX AND
      DEAL_II_HAVE_CXX11_TUPLE AND
      DEAL_II_HAVE_CXX11_TYPE_TRAITS )
    SET(DEAL_II_HAVE_CXX11 TRUE)
  ENDIF()

  RESET_CMAKE_REQUIRED()
ENDIF()

#
# Set up a configuration option for C++11 support:
#

OPTION(DEAL_II_WITH_CXX11
  "Compile deal.II using C++11 language standard."
  ${DEAL_II_HAVE_CXX11}
  )

#
# Bail out if user requested C++11 support (DEAL_II_WITH_CXX11) but support
# is not available due to above tests (DEAL_II_HAVE_CXX11):
#

IF(DEAL_II_WITH_CXX11 AND NOT DEAL_II_HAVE_CXX11)
  MESSAGE(FATAL_ERROR "\n"
    "C++11 support was requested (DEAL_II_WITH_CXX11=TRUE) but is not
    supported by the current compiler.\n"
    "Please disable C++11 support, i.e. configure with\n"
    "    -DDEAL_II_WITH_CXX11=FALSE,\n"
    "or use a different compiler, instead. (If the compiler flag for C++11 "
    "support differs from \"-std=c++0x\" or \"-std=c++11\", a suitable "
    "compiler flag has to be specified manually.\n\n"
    )
ENDIF()

#
# Set up C++11 support:
#

IF(DEAL_II_WITH_CXX11)
  ADD_FLAGS(DEAL_II_CXX_FLAGS "${DEAL_II_CXX11_FLAG}")
  MESSAGE(STATUS "DEAL_II_WITH_CXX11 successfully set up")

  PUSH_TEST_FLAG("${DEAL_II_CXX11_FLAG}")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <type_traits>
    int main(){ std::is_trivially_copyable<int> bob; }
    "
    DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE)
  POP_TEST_FLAG()
ELSE()
  MESSAGE(STATUS "DEAL_II_WITH_CXX11 disabled")
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
  HAVE_ISNAN)


CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  int main(){ double d=0; _isnan (d); return 0; }
  "
  HAVE_UNDERSCORE_ISNAN)


CHECK_CXX_SOURCE_COMPILES(
  "
  #include <cmath>
  int main(){ double d=0; std::isfinite (d); return 0; }
  "
  DEAL_II_HAVE_ISFINITE)

