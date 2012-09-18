FOREACH(flags ${deal_ii_used_flags})
  #
  # Append the saved cache variable ${flags}_SAVED at the end of ${flags}
  #
  SET(${flags} "${${flags}} ${${flags}_SAVED}")
ENDFOREACH()

MESSAGE("


deal.II successfully configured!


Compiler Flags:

Flags used by the compiler during all build types:

    CMAKE_C_FLAGS:   ${CMAKE_C_FLAGS}
    CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}

Additional flags used by the compiler during all release builds:

    CMAKE_C_FLAGS_RELEASE:   ${CMAKE_C_FLAGS_RELEASE}
    CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}

Additional flags used by the compiler during all debug builds:

    CMAKE_C_FLAGS_DEBUG:   ${CMAKE_C_FLAGS_DEBUG}
    CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}

(Note: Flags set with ccmake or the command line will be appended at the end
of the default configuration)


Configured Features:
")

GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    MESSAGE("    ${var} = ${${var}}")
  ENDIF()
ENDFOREACH()

MESSAGE("
TODO: Tell something about contrib/external

")
