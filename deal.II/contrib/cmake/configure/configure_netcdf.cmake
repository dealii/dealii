#
# Configuration for the netcdf library:
#

MACRO(FIND_FEATURE_NETCDF_EXTERNAL var)

  FIND_PACKAGE(NETCDF)

  IF(NETCDF_FOUND)
    SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_NETCDF_EXTERNAL var)

  INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIR})
  LIST(APPEND deal_ii_external_libraries ${NETCDF_LIBRARY})
  SET(HAVE_LIBNETCDF TRUE)

  SET(${var} TRUE)
ENDMACRO()


MACRO(CONFIGURE_FEATURE_NETCDF_ERROR_MESSAGE)

  MESSAGE(SEND_ERROR "
Could not find the netcdf library!

Please ensure that the netcdf library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

")

ENDMACRO()


CONFIGURE_FEATURE(NETCDF)
