#
# Check for a compiler bug.
#
# Usage:
#     CHECK_CXX_COMPILER_BUG(source var),
#
# where source is a snipped of source code and var is a variable that will
# be set to true if the source could not be compiled and linked successfully.
# (This just inverts the logic of CHECK_CXX_SOURCE_COMPILES.)
#

MACRO(CHECK_CXX_COMPILER_BUG source var)
  CHECK_CXX_SOURCE_COMPILES(
    "${source}"
    ${var}_OK)

  IF(${var}_OK)
    MESSAGE(STATUS "Test successful, do not define ${var}")
  ELSE()
    MESSAGE(STATUS "Test unsuccessful, define ${var}")
    SET(${var} 1)
  ENDIF()
ENDMACRO()

