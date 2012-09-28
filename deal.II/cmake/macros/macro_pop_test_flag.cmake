#
# A small macro used in the platform checks to remove the right most flag in
# CMAKE_REQUIRED_FLAGS
#
# We assume that the flags in CMAKE_REQUIRED_FLAGS are space separated
#
# Usage:
#     POP_TEST_FLAG()
#

MACRO(POP_TEST_FLAG)
  SET(CMAKE_REQUIRED_FLAGS " ${CMAKE_REQUIRED_FLAGS}")
  STRING(REGEX REPLACE " -[^ ]+$" ""
    CMAKE_REQUIRED_FLAGS
    "${CMAKE_REQUIRED_FLAGS}"
    )
  STRING(STRIP "${CMAKE_REQUIRED_FLAGS}" CMAKE_REQUIRED_FLAGS)
ENDMACRO()

