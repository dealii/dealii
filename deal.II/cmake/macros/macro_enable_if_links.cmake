#
# Tests whether it is possible to compile and link a dummy program with a
# given flag.
# If so, add it to variable.
#
# Usage:
#     ENABLE_IF_LINKS(variable flag)
#

MACRO(ENABLE_IF_LINKS variable flag)
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${flag}")
  CHECK_CXX_SOURCE_COMPILES(
  "
  int main() { return 0; }
  "
  DEAL_II_HAVE_FLAG_${flag}
  )
  STRIP_FLAG(CMAKE_REQUIRED_FLAGS "${flag}")
  IF(DEAL_II_HAVE_FLAG_${flag})
    SET(${variable} "${${variable}} ${flag}")
  ENDIF()
ENDMACRO()

