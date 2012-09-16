#
# Configuration for the netcdf library:
#

MACRO(FIND_FEATURE_METIS_EXTERNAL var)

  FIND_PACKAGE(METIS)

  IF(METIS_FOUND AND METIS_MAJOR GREATER 4)
    SET(${var} TRUE)
  ELSE()
    MESSAGE(STATUS
      "Could not find a sufficient modern metis installation: Version 5.x
      required."
      )
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_METIS_EXTERNAL var)

  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  LIST(APPEND deal_ii_external_libraries ${METIS_LIBRARY})

  SET(DEAL_II_USE_METIS TRUE)

  SET(${var} TRUE)
ENDMACRO()


MACRO(CONFIGURE_FEATURE_METIS_ERROR_MESSAGE)

  MESSAGE(SEND_ERROR "
Could not find the metis library!

Please ensure that the metis version 5.x is installed on your computer.
If the library is not at a default location, either provide some hints
via environment variables:
METIS_LIBRARY_DIR METIS_INCLUDE_DIR
Or set the relevant variables by hand in ccmake.

")

ENDMACRO()


CONFIGURE_FEATURE(METIS)
