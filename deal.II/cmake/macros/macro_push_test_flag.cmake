#
# A small macro used in the platform checks to easily add a flag to
# CMAKE_REQUIRED_FLAGS
#
# Usage:
#     PUSH_TEST_FLAG("flag")
#

MACRO(PUSH_TEST_FLAG flag)

  SET(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${flag}")
  STRING(STRIP "${CMAKE_REQUIRED_FLAGS}" CMAKE_REQUIRED_FLAGS)

ENDMACRO()

