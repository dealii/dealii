cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

#
# Make sure that we have been called with a proper deal.II source directory
# layout:
#
set(DEAL_II_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../..")
if (EXISTS "${DEAL_II_SOURCE_DIR}/examples")
  #
  # Create CMake file structure in the current binary dir:
  #
  set(_source_directory "${CMAKE_CURRENT_BINARY_DIR}/source")
  set(_binary_directory "${CMAKE_CURRENT_BINARY_DIR}/binary")
  message(STATUS "Temporary source directory: ${_source_directory}")
  message(STATUS "Temporary binary directory: ${_binary_directory}")
  file(MAKE_DIRECTORY "${_source_directory}")
  file(WRITE "${_source_directory}/CMakeLists.txt" "deal_ii_pickup_tests(examples)")
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/example_test.h" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

  enable_testing()
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/CTestFile.cmake" "subdirs(binary)")

  #
  # Create tests from diffs:
  #
  execute_process(
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/update.sh" "${DEAL_II_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}"
    WORKING_DIRECTORY "${_source_directory}"
    )

  add_subdirectory(${_source_directory} ${_binary_directory})

  #
  # Create a top-level diff target:
  #
  add_custom_target(update_diffs
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/update_diffs.sh" "${DEAL_II_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}"
    WORKING_DIRECTORY "${_source_directory}"
    )
endif()
