cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

#
# Setting up example tests is a bit tricky because we need to patch the
# current source file of the examples to make them suitable for the
# testsuite (sanitize output, reduce run time).
#
# We thus create a set of "source" and "binary" sub directories in our
# CMAKE_CURRENT_BINARY_DIR that we can manipulate at will. We then populate
# the source directory with all (dynamically) patched step-*.cc files,
# create symlinks to all output files, and put the CMakelists.txt.in file
# in place.
#
# As a last step we call into the subdirectory and process all tests:
#

set(DEAL_II_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
if (EXISTS "${DEAL_II_SOURCE_DIR}/examples")
  #
  # Create CMake file structure in the current binary dir:
  #
  set(_base_directory   "${CMAKE_CURRENT_SOURCE_DIR}")
  set(_source_directory "${CMAKE_CURRENT_BINARY_DIR}/source")
  set(_binary_directory "${CMAKE_CURRENT_BINARY_DIR}/binary")
  message(STATUS "Temporary source directory: ${_source_directory}")
  message(STATUS "Temporary binary directory: ${_binary_directory}")
  file(MAKE_DIRECTORY "${_source_directory}")

  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/example_test.h" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt.in" "${_source_directory}/CMakeLists.txt"
    COPYONLY
    )

  enable_testing()
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/CTestFile.cmake" "subdirs(binary)")

  #
  # and call into the freshly created "source" directory:
  #
  add_subdirectory(${_source_directory} ${_binary_directory})
endif()
