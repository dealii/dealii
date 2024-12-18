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

#
# Configuration for the boost library:
#

if(NOT FEATURE_ZLIB_PROCESSED)
  message(FATAL_ERROR "\n"
    "Internal build system error: The configuration of "
    "DEAL_II_WITH_BOOST depends on "
    "DEAL_II_WITH_ZLIB, but configure_feature(BOOST) "
    "was called before configure_feature(ZLIB).\n\n"
    )
endif()


set(DEAL_II_WITH_BOOST ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for boost." FORCE
  )


macro(feature_boost_find_external var)
  find_package(DEAL_II_BOOST)

  if(BOOST_FOUND)
    set(${var} TRUE)

    if(DEAL_II_WITH_ZLIB)
      #
      # Test that Boost.Iostreams is usable.
      #
      reset_cmake_required()
      list(APPEND CMAKE_REQUIRED_LIBRARIES ${BOOST_TARGETS})
      list(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIRS})

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
      if(NOT ${BOOST_IOSTREAMS_USABLE})
        message(STATUS
          "DEAL_II_WITH_ZLIB=ON requires Boost.Iostreams to be compiled "
          "with zlib support but a simple test failed! "
          "Therefore, the bundled boost package is used."
          )
        set(BOOST_ADDITIONAL_ERROR_STRING
          "DEAL_II_WITH_ZLIB=ON requires Boost.Iostreams to be compiled "
          "with zlib support but a simple test failed! "
          )
        set(${var} FALSE)
      endif()
      reset_cmake_required()
    endif() # DEAL_II_WITH_ZLIB

    if(${BOOST_VERSION} VERSION_LESS 1.74.0 AND DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
      #
      # Test that Boost.Serialization is usable.
      #
      if(NOT DEFINED BOOST_SERIALIZATION_USABLE OR NOT ${BOOST_SERIALIZATION_USABLE})
        # Only run this check if it hasn't successfully run previously.
        message(STATUS "Performing Test BOOST_SERIALIZATION_USABLE")

        set(_binary_test_dir ${CMAKE_CURRENT_BINARY_DIR}/cmake/configure/TestBoostBugWorkdir)

        set(_flags "${DEAL_II_CXX_FLAGS}")
        strip_flag(_flags "-Werror")

        file(REMOVE_RECURSE ${_binary_test_dir})
        file(MAKE_DIRECTORY ${_binary_test_dir})
        execute_process(
          COMMAND ${CMAKE_COMMAND}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            "-DCMAKE_CXX_FLAGS=${_flags}"
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
        if(${_result} EQUAL 0)
          execute_process(
            COMMAND ${CMAKE_COMMAND} --build . --target run
            WORKING_DIRECTORY ${_binary_test_dir}
            RESULT_VARIABLE _result
            OUTPUT_QUIET
            ERROR_QUIET
            )
        endif()
        if(${_result} EQUAL 0)
          message(STATUS "Performing Test BOOST_SERIALIZATION_USABLE - Success")
          set(BOOST_SERIALIZATION_USABLE TRUE CACHE INTERNAL "")
        else()
          message(STATUS "Performing Test BOOST_SERIALIZATION_USABLE - Failed")
          set(BOOST_SERIALIZATION_USABLE FALSE)
        endif()
      endif()

      if(NOT ${BOOST_SERIALIZATION_USABLE})
        message(STATUS
          "The externally provided Boost.Serialization library "
          "failed to pass a crucial test. \n"
          "Therefore, the bundled boost package is used. \n"
          "The configured testing project can be found at \n"
          "${_binary_test_dir}"
          )
        set(BOOST_ADDITIONAL_ERROR_STRING
          "The externally provided Boost.Serialization library "
          "failed to pass a crucial test."
          )
        set(${var} FALSE)
      endif()
    endif() # DEAL_II_ALLOW_PLATFORM_INTROSPECTION
  endif()
endmacro()


macro(feature_boost_configure_external)
  #
  # Avoid a number of warnings:
  #
  enable_if_supported(DEAL_II_WARNING_FLAGS "-Wno-unused-local-typedefs")

  #
  # At least BOOST 1.74 has the problem that some of the BOOST headers
  # include other BOOST headers that are deprecated, and this then leads to
  # warnings. That's rather annoying.
  #
  list(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIRS})

  check_cxx_compiler_bug(
    "
    #define BOOST_CONFIG_HEADER_DEPRECATED_HPP_INCLUDED
    #define BOOST_HEADER_DEPRECATED(a) _Pragma(\"GCC error \\\"stop compilation\\\"\");
    #include <boost/geometry/index/rtree.hpp>
    int main() { return 0; }
    "
    DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS)
  reset_cmake_required()
endmacro()


configure_feature(BOOST)


#
# DEAL_II_WITH_BOOST is always required.
#
if(NOT DEAL_II_WITH_BOOST)
  if(DEAL_II_FEATURE_AUTODETECTION)
    feature_error_message("BOOST")
  else()
    message(FATAL_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_BOOST required, but set to OFF!\n\n"
      )
  endif()
endif()
