cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

add_executable(dummy dummy.cc)
set(TEST_TARGET dummy)

deal_ii_pickup_tests()


#
# And a configure test:
#

foreach(_build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${_build} _build_lowercase)

  set(_target a-framework-configure.${_build_lowercase})
  set(_test a-framework/configure.${_build_lowercase})

  # Respect TEST_PICKUP_REGEX:
  if( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
      _test MATCHES "${TEST_PICKUP_REGEX}"  )
    add_custom_target(${_target}
      COMMAND echo "${_test}: CONFIGURE failed." && exit 1
      )
    add_test(NAME ${_test}
      COMMAND ${CMAKE_COMMAND} -DTRGT=${_target} -DTEST=${_test}
        -DBINARY_DIR=${CMAKE_BINARY_DIR} -DEXPECT=CONFIGURE
        -P ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    set_tests_properties(${_test} PROPERTIES
      LABEL "${_category}"
      TIMEOUT ${TEST_TIME_LIMIT}
      )
  endif()
endforeach()
