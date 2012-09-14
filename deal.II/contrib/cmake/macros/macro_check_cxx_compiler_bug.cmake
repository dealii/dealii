INCLUDE(CheckCXXSourceCompiles)

MACRO(CHECK_CXX_COMPILER_BUG source var)

  #
  # Check for a compiler bug, i.e. if source does not compile, define var
  # This just inverts the logic of CHECK_CXX_SOURCE_COMPILES.
  #

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
