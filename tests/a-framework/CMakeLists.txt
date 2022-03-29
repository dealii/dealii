CMAKE_MINIMUM_REQUIRED(VERSION 3.3.0)
INCLUDE(../setup_testsubproject.cmake)
PROJECT(testsuite CXX)

ADD_EXECUTABLE(dummy.release dummy.cc)
ADD_EXECUTABLE(dummy.debug dummy.cc)
SET(TEST_TARGET_RELEASE dummy.release)
SET(TEST_TARGET_DEBUG dummy.debug)

DEAL_II_PICKUP_TESTS()

#
# Limit concurrency between the two parameter file tests.
#
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build)
  IF(TARGET parameter_file_2.${_build}.diff AND
     TARGET parameter_file_1.${_build}.diff )
    SET_TESTS_PROPERTIES(a-framework/parameter_file_2.${_build} PROPERTIES
      DEPENDS a-framework/parameter_file_1.${_build}
      )
  ENDIF()
ENDFOREACH()


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
