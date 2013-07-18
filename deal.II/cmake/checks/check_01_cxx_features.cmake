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
#   DEAL_II_HAVE_CXX11_FLAG
#   DEAL_II_CXX11_FLAG
#   DEAL_II_CAN_USE_CXX11
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
# See if there is a compiler flag to enable C++11 features
#
IF(NOT DEFINED DEAL_II_HAVE_CXX11_FLAG)
  FOREACH(_test_flag
      "-std=c++11"
      "-std=c++0x"
      )
    CHECK_CXX_COMPILER_FLAG("${_test_flag}" DEAL_II_HAVE_CXX11_FLAG)

    IF(DEAL_II_HAVE_CXX11_FLAG)
      # We have found a CXX11_FLAG that the compiler understands
      SET(DEAL_II_CXX11_FLAG "${_test_flag}" CACHE INTERNAL "")
      BREAK()
    ELSE()
      # Remove test result from cache and try the next flag in the list
      UNSET(DEAL_II_HAVE_CXX11_FLAG CACHE)
    ENDIF()
  ENDFOREACH()
ENDIF()


IF(DEAL_II_HAVE_CXX11_FLAG)

  # Set CMAKE_REQUIRED_FLAGS for the unit tests
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
  # TODO: This test will only succeed on platforms where "-pthread" is
  #       recognized. But this isn't easily fixable:
  #       configure_threads.cmake which will determine and setup threads
  #       has to be called later...
  #
  IF(NOT CMAKE_CROSSCOMPILING) # Todo: Is it better to use DEAL_II_ALLOW_PLATFORM_INTROSPECTION here?
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
    SET_IF_EMPTY(DEAL_II_HAVE_CXX11_THREAD TRUE)
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

    MESSAGE(STATUS "Sufficient C++11 support. Enabling ${DEAL_II_CXX11_FLAG}.")

    SET(DEAL_II_CAN_USE_CXX11 TRUE)

    ADD_FLAGS(CMAKE_CXX_FLAGS "${DEAL_II_CXX11_FLAG}")

  ELSE()
    MESSAGE(STATUS "Insufficient C++11 support. Disabling ${DEAL_II_CXX11_FLAG}.")
  ENDIF()

#
# Currently unused
#
#  IF(DEAL_II_CAN_USE_CXX11)
#    #
#    # Also test for a couple of C++11 things that we don't use in the
#    # library but that users may want to use in their applications and that
#    # we might want to test in the testsuite
#    #
#    # TODO: Actually we have to export the test results somehow. :-]
#    #
#
#    CHECK_CXX_SOURCE_COMPILES(
#      "
#      #include <vector>
#      std::vector<int> v;
#      int main(){ auto i = v.begin(); *i; return 0;}
#      "
#      DEAL_II_HAVE_CXX11_AUTO_TYPE)
#
#    CHECK_CXX_SOURCE_COMPILES(
#      "
#      #include <vector>],
#      std::vector<int> v;
#      int main(){ for (std::vector<int>::iterator i : v) *i; return 0;}
#      "
#      DEAL_II_HAVE_CXX11_RANGE_BASED_FOR)
#
#    IF( DEAL_II_HAVE_CXX11_AUTO_TYPE AND
#        DEAL_II_HAVE_CXX11_RANGE_BASED_FOR )
#
#      MESSAGE(STATUS "Additional C++11 support available.")
#
#      SET(DEAL_II_CAN_USE_ADDITIONAL_CXX1X_FEATURES)
#    ENDIF()
#
#  ENDIF()

  POP_TEST_FLAG()

ELSE()
    MESSAGE(STATUS "Insufficient C++11 support. Disabling ${DEAL_II_CXX11_FLAG}.")
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

