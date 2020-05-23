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
# Configuration for the boost library:
#

IF(NOT FEATURE_ZLIB_PROCESSED)
  MESSAGE(FATAL_ERROR "\n"
    "Internal build system error: The configuration of "
    "DEAL_II_WITH_BOOST depends on "
    "DEAL_II_WITH_ZLIB, but CONFIGURE_FEATURE(BOOST) "
    "was called before CONFIGURE_FEATURE(ZLIB).\n\n"
    )
ENDIF()


SET(DEAL_II_WITH_BOOST ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for boost." FORCE
  )


MACRO(FEATURE_BOOST_CONFIGURE_COMMON)
  #
  # Boost version 1.62 - 1.63 checks for the availability of "emplace_hint"
  # incorrectly: It tests for the preprocessor define
  # BOOST_NO_CXX11_HDR_UNORDERED_MAP in .../boost/serialization/map.h
  # thinking that that this define is characteristic for the presence of
  # std::(multi)map::emplace_hint. This is generally correct, except for
  # GCC before 4.8, for which the preprocessor variable is defined, but the
  # function does not exist [1].
  #
  # Thus, simply define a BOOST_NO_CXX11_HDR_UNORDERED_MAP if the gcc
  # compiler version is less than 4.8.
  #
  # [1] https://svn.boost.org/trac/boost/ticket/12755
  #
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
      CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
    LIST(APPEND BOOST_DEFINITIONS "BOOST_NO_CXX11_HDR_UNORDERED_MAP")
    LIST(APPEND BOOST_USER_DEFINITIONS "BOOST_NO_CXX11_HDR_UNORDERED_MAP")
  ENDIF()

  # Some standard library implementations do not implement std::auto_ptr
  # (anymore) which was deprecated for C++11 and removed in the C++17 standard.
  # Older boost versions can't know about this but provide a possibility to
  # circumvent the issue. Hence, we just check ourselves.
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
  IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    ADD_FLAGS(CMAKE_REQUIRED_FLAGS "/WX /EHsc")
  ELSE()
    ADD_FLAGS(CMAKE_REQUIRED_FLAGS "-Werror")
  ENDIF()
  # The configure function is called only once. In case an externally provided
  # boost library is detected, BOOST_INCLUDE_DIRS contains the include paths to
  # be used and BOOST_BUNDLED_INCLUDE_DIRS is empty. For the bundled library, it
  # is the other way around.
  LIST(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIRS} ${BOOST_BUNDLED_INCLUDE_DIRS})

  # In case, the boost library already sets BOOST_NO_AUTO_PTR we report
  # DEAL_II_HAS_AUTO_PTR to be true to avoid redefining the macro.
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <memory>
    #include <boost/config.hpp>

    int main()
    {
    #ifndef BOOST_NO_AUTO_PTR
      int *i = new int;
      std::auto_ptr<int> x(i);
    #endif
      return 0;
    }
    "
    DEAL_II_HAS_AUTO_PTR)

  RESET_CMAKE_REQUIRED()

  IF( NOT DEAL_II_HAS_AUTO_PTR )
    LIST(APPEND BOOST_DEFINITIONS "BOOST_NO_AUTO_PTR")
    LIST(APPEND BOOST_USER_DEFINITIONS "BOOST_NO_AUTO_PTR")
  ENDIF()

  ENABLE_IF_SUPPORTED(BOOST_CXX_FLAGS "-Wno-unused-local-typedefs")
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_BUNDLED)
  #
  # Add rt to the link interface as well, boost/chrono needs it.
  #
  IF(NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    FIND_SYSTEM_LIBRARY(rt_LIBRARY NAMES rt)
    MARK_AS_ADVANCED(rt_LIBRARY)
    IF(NOT rt_LIBRARY MATCHES "-NOTFOUND")
      SET(BOOST_LIBRARIES ${rt_LIBRARY})
    ENDIF()
  ENDIF()

  # We need to set this path before calling the configure function
  # to be able to use the include paths in the checks.
  SET(BOOST_BUNDLED_INCLUDE_DIRS ${BOOST_FOLDER}/include)

  FEATURE_BOOST_CONFIGURE_COMMON()

  IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
    #
    # Bundled boost tries to (dl)open itself as a dynamic library on
    # Windows. Disable this undesired behavior by exporting
    # BOOST_ALL_NO_LIB on Windows platforms (for bundled boost).
    #
    LIST(APPEND BOOST_DEFINITIONS "BOOST_ALL_NO_LIB")
    LIST(APPEND BOOST_USER_DEFINITIONS "BOOST_ALL_NO_LIB")
  ENDIF()
ENDMACRO()

MACRO(FEATURE_BOOST_FIND_EXTERNAL var)
  FIND_PACKAGE(BOOST)

  IF(BOOST_FOUND)
    SET(${var} TRUE)

    IF(DEAL_II_WITH_ZLIB)
      #
      # Test that Boost.Iostreams is usable.
      #
      RESET_CMAKE_REQUIRED()
      ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")
      LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${BOOST_LIBRARIES})
      LIST(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIRS})

      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <string>
        #include <boost/iostreams/device/back_inserter.hpp>
        #include <boost/iostreams/filter/gzip.hpp>
        #include <boost/iostreams/filtering_stream.hpp>

        int main()
        {
          std::string decompressed_buffer;
          char test[1] = {'c'};

          boost::iostreams::filtering_ostream decompressing_stream;
          decompressing_stream.push(boost::iostreams::gzip_decompressor());
          decompressing_stream.push(boost::iostreams::back_inserter(decompressed_buffer));
          decompressing_stream.write (test, 1);
        }
        "
        BOOST_IOSTREAMS_USABLE
        )
      IF(NOT ${BOOST_IOSTREAMS_USABLE})
        MESSAGE(STATUS
          "DEAL_II_WITH_ZLIB=ON requires Boost.Iostreams to be compiled "
          "with zlib support but a simple test failed! "
          "Therefore, the bundled boost package is used."
          )
        SET(BOOST_ADDITIONAL_ERROR_STRING
          "DEAL_II_WITH_ZLIB=ON requires Boost.Iostreams to be compiled "
          "with zlib support but a simple test failed! "
          )
        SET(${var} FALSE)
      ENDIF()
      RESET_CMAKE_REQUIRED()
    ENDIF() # DEAL_II_WITH_ZLIB

    IF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
      #
      # Test that Boost.Serialization is usable.
      #
      IF(NOT DEFINED BOOST_SERIALIZATION_USABLE OR NOT ${BOOST_SERIALIZATION_USABLE})
        # Only run this check if it hasn't successfully run previously.
        MESSAGE(STATUS "Performing Test BOOST_SERIALIZATION_USABLE")

        SET(_binary_test_dir ${CMAKE_CURRENT_BINARY_DIR}/cmake/configure/TestBoostBugWorkdir)

        FILE(REMOVE_RECURSE ${_binary_test_dir})
        FILE(MAKE_DIRECTORY ${_binary_test_dir})
        EXECUTE_PROCESS(
          COMMAND ${CMAKE_COMMAND}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            "-DCMAKE_CXX_FLAGS=${DEAL_II_CXX_FLAGS}"
            "-DCMAKE_EXE_LINKER_FLAGS=${DEAL_II_LINKER_FLAGS}"
            "-DCMAKE_SHARED_LINKER_FLAGS=${DEAL_II_LINKER_FLAGS}"
            "-DBOOST_INCLUDE_DIRS=${BOOST_INCLUDE_DIRS}"
            "-DBOOST_LIBRARIES=${BOOST_LIBRARIES}"
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure/TestBoostBug
          WORKING_DIRECTORY ${_binary_test_dir}
          RESULT_VARIABLE _result
          OUTPUT_QUIET
          ERROR_QUIET
          )
        IF(${_result} EQUAL 0)
          EXECUTE_PROCESS(
            COMMAND ${CMAKE_COMMAND} --build . --target run
            WORKING_DIRECTORY ${_binary_test_dir}
            RESULT_VARIABLE _result
            OUTPUT_QUIET
            ERROR_QUIET
            )
        ENDIF()
        IF(${_result} EQUAL 0)
          MESSAGE(STATUS "Performing Test BOOST_SERIALIZATION_USABLE - Success")
          SET(BOOST_SERIALIZATION_USABLE TRUE CACHE INTERNAL "")
        ELSE()
          MESSAGE(STATUS "Performing Test BOOST_SERIALIZATION_USABLE - Failed")
          SET(BOOST_SERIALIZATION_USABLE FALSE)
        ENDIF()
      ENDIF()

      IF(NOT ${BOOST_SERIALIZATION_USABLE})
        MESSAGE(STATUS
          "The externally provided Boost.Serialization library "
          "failed to pass a crucial test. \n"
          "Therefore, the bundled boost package is used. \n"
          "The configured testing project can be found at \n"
          "${_binary_test_dir}"
          )
        SET(BOOST_ADDITIONAL_ERROR_STRING
          "The externally provided Boost.Serialization library "
          "failed to pass a crucial test."
          )
        SET(${var} FALSE)
      ENDIF()
    ENDIF() # DEAL_II_ALLOW_PLATFORM_INTROSPECTION
  ENDIF()
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_EXTERNAL)
  FEATURE_BOOST_CONFIGURE_COMMON()
ENDMACRO()


CONFIGURE_FEATURE(BOOST)


#
# DEAL_II_WITH_BOOST is always required.
#
IF(NOT DEAL_II_WITH_BOOST)
  IF(DEAL_II_FEATURE_AUTODETECTION)
    FEATURE_ERROR_MESSAGE("BOOST")
  ELSE()
    MESSAGE(FATAL_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_BOOST required, but set to OFF!.\n\n"
      )
  ENDIF()
ENDIF()
