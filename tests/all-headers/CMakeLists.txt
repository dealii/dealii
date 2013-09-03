#
# Header tests are special:
#
# Construct a list of all header files and build up a test that just tries
# to compile a simple worker (test_header.cc) that only includes the given
# header file. We omit linking to safe some time.
#
SET(_category all-headers)

FILE(GLOB_RECURSE _header
  ${CMAKE_SOURCE_DIR}/include/*.h
  )

FOREACH(_full_file ${_header})
  GET_FILENAME_COMPONENT(_file ${_full_file} NAME)
  GET_FILENAME_COMPONENT(_path ${_full_file} PATH)
  GET_FILENAME_COMPONENT(_path ${_path} NAME)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    SET(_test ${_category}-${_path}-${_file}.${_build_lowercase})

    #
    # Add an object library for each header file and build configuration:
    #
    ADD_LIBRARY(${_test} OBJECT EXCLUDE_FROM_ALL test_header.cc)
    SET_TARGET_PROPERTIES(${_test} PROPERTIES
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
      COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${_build}}"
      )
    SET_PROPERTY(TARGET ${_test} APPEND PROPERTY
      INCLUDE_DIRECTORIES
        "${CMAKE_BINARY_DIR}/include"
        "${CMAKE_SOURCE_DIR}/include"
        "${CMAKE_SOURCE_DIR}/include/deal.II/"
      )
    SET_PROPERTY(TARGET ${_test} APPEND PROPERTY
      COMPILE_DEFINITIONS
        HEADER=<deal.II/${_path}/${_file}>
      )

    #
    # And finally add the test:
    #
    ADD_TEST(NAME ${_category}/${_path}/${_file}.${_build_lowercase}
      COMMAND ${CMAKE_COMMAND}
        -DTEST=${_test}
        -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_SOURCE_DIR}/cmake/scripts/run_test.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    SET_TESTS_PROPERTIES(${_category}/${_path}/${_file}.${_build_lowercase} PROPERTIES
      LABEL "${_category}"
      TIMEOUT ${TEST_TIME_LIMIT}
      )
  ENDFOREACH()
ENDFOREACH()
