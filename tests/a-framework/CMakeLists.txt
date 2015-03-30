CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)
INCLUDE(../setup_testsubproject.cmake)
PROJECT(testsuite CXX)
DEAL_II_PICKUP_TESTS()

#
# And a configure test:
#

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build_lowercase)

  SET(_target a-framework-configure.${_build_lowercase})
  SET(_test a-framework/configure.${_build_lowercase})

  # Respect TEST_PICKUP_REGEX:
  IF( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
      _test MATCHES "${TEST_PICKUP_REGEX}"  )
    ADD_CUSTOM_TARGET(${_target}
      COMMAND echo "${_test}: CONFIGURE failed." && exit 1
      )
    ADD_TEST(NAME ${_test}
      COMMAND ${CMAKE_COMMAND} -DTRGT=${_target} -DTEST=${_test}
        -DBINARY_DIR=${CMAKE_BINARY_DIR} -DEXPECT=CONFIGURE
        -P ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    SET_TESTS_PROPERTIES(${_test} PROPERTIES
      LABEL "${_category}"
      TIMEOUT ${TEST_TIME_LIMIT}
      )
  ENDIF()
ENDFOREACH()
