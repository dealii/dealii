#
# Configuration for the zlib library:
#

MACRO(FIND_FEATURE_ZLIB_EXTERNAL var)

  FIND_PACKAGE(ZLIB)

  IF(ZLIB_FOUND)
    SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_ZLIB_EXTERNAL var)

  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
  LIST(APPEND deal_ii_external_libraries ${ZLIB_LIBRARIES})
  SET(HAVE_LIBZ TRUE)

  SET(${var} TRUE)

ENDMACRO()


MACRO(CONFIGURE_FEATURE_ZLIB_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "
Could not find the zlib library!

Please ensure that the zlib library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

")
ENDMACRO()


CONFIGURE_FEATURE(ZLIB)
