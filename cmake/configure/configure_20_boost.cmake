## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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


macro(FEATURE_BOOST_CONFIGURE_COMMON)
  # Some standard library implementations do not implement std::auto_ptr
  # (anymore) which was deprecated for C++11 and removed in the C++17 standard.
  # Older boost versions can't know about this but provide a possibility to
  # circumvent the issue. Hence, we just check ourselves.
  if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_flags(CMAKE_REQUIRED_FLAGS "/WX /EHsc")
  else()
    add_flags(CMAKE_REQUIRED_FLAGS "-Werror")
  endif()
  # The configure function is called only once. In case an externally provided
  # boost library is detected, BOOST_INCLUDE_DIRS contains the include paths to
  # be used and BOOST_BUNDLED_INCLUDE_DIRS is empty. For the bundled library, it
  # is the other way around.
  list(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIRS} ${BOOST_BUNDLED_INCLUDE_DIRS})

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

  reset_cmake_required()

  # Fix some problems by defining some additional preprocessor symbols.
  # Ultimately these are added into DEAL_II_DEFINITIONS. They are separate
  # here so that they show up in detailed.log under DEAL_II_WITH_BOOST as,
  # logically, they are part of our boost configuration.
  if(NOT DEAL_II_HAS_AUTO_PTR)
    list(APPEND BOOST_DEFINITIONS "BOOST_NO_AUTO_PTR")
  endif()

  enable_if_supported(BOOST_CXX_FLAGS "-Wno-unused-local-typedefs")

  # At least BOOST 1.74 has the problem that some of the BOOST headers
  # include other BOOST headers that are deprecated, and this then leads to
  # warnings. That's rather annoying.

  # The configure function is called only once. In case an externally provided
  # boost library is detected, BOOST_INCLUDE_DIRS contains the include paths to
  # be used and BOOST_BUNDLED_INCLUDE_DIRS is empty. For the bundled library, it
  # is the other way around.
  list(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIRS} ${BOOST_BUNDLED_INCLUDE_DIRS})

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


macro(FEATURE_BOOST_CONFIGURE_BUNDLED)
  #
  # Add rt to the link interface as well, boost/chrono needs it.
  #
  if(NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    find_system_library(rt_LIBRARY NAMES rt)
    mark_as_advanced(rt_LIBRARY)
    if(NOT rt_LIBRARY MATCHES "-NOTFOUND")
      set(BOOST_LIBRARIES ${rt_LIBRARY})
    endif()
  endif()

  # We need to set this path before calling the configure function
  # to be able to use the include paths in the checks.
  set(BOOST_BUNDLED_INCLUDE_DIRS ${BOOST_FOLDER}/include)
  #
  # We still need the version information, which is set up in the FindBoost
  # module in the non-bundled case:
  #
  file(STRINGS "${BOOST_BUNDLED_INCLUDE_DIRS}/boost/version.hpp"
    BOOST_VERSION_STRING
    REGEX "#define.*BOOST_VERSION")

  string(REGEX REPLACE "^.*BOOST_VERSION.* ([0-9]+).*" "\\1"
    BOOST_VERSION_NUMBER "${BOOST_VERSION_STRING}"
    )
  math(EXPR Boost_MAJOR_VERSION "${BOOST_VERSION_NUMBER} / 100000")
  math(EXPR Boost_MINOR_VERSION "${BOOST_VERSION_NUMBER} / 100 % 1000")
  math(EXPR Boost_SUBMINOR_VERSION "${BOOST_VERSION_NUMBER} % 100")

  FEATURE_BOOST_CONFIGURE_COMMON()

  if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    #
    # Bundled boost tries to (dl)open itself as a dynamic library on
    # Windows. Disable this undesired behavior by exporting
    # BOOST_ALL_NO_LIB on Windows platforms (for bundled boost).
    #
    list(APPEND BOOST_DEFINITIONS "BOOST_ALL_NO_LIB")
  endif()
endmacro()

macro(FEATURE_BOOST_FIND_EXTERNAL var)
  find_package(DEAL_II_BOOST)

  if(BOOST_FOUND)
    set(${var} TRUE)

    if(DEAL_II_WITH_ZLIB)
      #
      # Test that Boost.Iostreams is usable.
      #
      reset_cmake_required()
      list(APPEND CMAKE_REQUIRED_LIBRARIES ${BOOST_LIBRARIES})
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

        file(REMOVE_RECURSE ${_binary_test_dir})
        file(MAKE_DIRECTORY ${_binary_test_dir})
        execute_process(
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


macro(FEATURE_BOOST_CONFIGURE_EXTERNAL)
  FEATURE_BOOST_CONFIGURE_COMMON()
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
