CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
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
# We only do this for CMake-3.4 or newer that has policy CMP0064 and
# supports the IF(TEST ...) operator. For older CMake versions that do not
# support the IF(TEST ...) operator we cannot much do (and people have to
# live with spurious test failures.)
#
IF(POLICY CMP0064)
  CMAKE_POLICY(SET CMP0064 NEW)
  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build)
    IF( TEST a-framework/parameter_file_2.${_build} AND
        TEST a-framework/parameter_file_1.${_build} )
      SET_TESTS_PROPERTIES(a-framework/parameter_file_2.${_build} PROPERTIES
        DEPENDS a-framework/parameter_file_1.${_build}
        )
    ENDIF()
  ENDFOREACH()
ENDIF()


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
