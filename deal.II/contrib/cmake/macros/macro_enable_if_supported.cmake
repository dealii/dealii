MACRO(ENABLE_IF_SUPPORTED variable flag)
  CHECK_CXX_COMPILER_FLAG(
    "${flag}"
    DEAL_II_HAVE_FLAG_${flag})
  IF(DEAL_II_HAVE_FLAG_${flag})
    SET(${variable} "${${variable}} ${flag}")
  ENDIF()
ENDMACRO()
