#
# Configuration for the netcdf library:
#

MACRO(FEATURE_METIS_FIND_EXTERNAL var)

  FIND_PACKAGE(METIS)

  IF(METIS_FOUND AND METIS_MAJOR GREATER 4)
    SET(${var} TRUE)
  ELSE()
    MESSAGE(WARNING "
Could not find a sufficient modern metis installation:
Version 5.x required!

Please ensure that a sufficiently modern version of metis is installed on
your computer. If the library is not at a default location, either provide
some hints for the autodetection, or set the relevant variables by hand in
ccmake.

"
      )
  ENDIF()

ENDMACRO()


MACRO(FEATURE_METIS_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  LIST(APPEND deal_ii_external_libraries ${METIS_LIBRARY})

  SET(DEAL_II_USE_METIS TRUE)

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(METIS)
