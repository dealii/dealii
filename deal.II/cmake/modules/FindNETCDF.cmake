#
# Try to find the NETCDF library
#

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(NETCDF_INCLUDE_DIR netcdf.hh
)

#
# TODO: netcdf might externally depend on hdf5. Check and fix this.
#

FIND_LIBRARY(NETCDF_LIBRARY
  NAMES netcdf_c++ netcdf_cpp
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(NETCDF DEFAULT_MSG NETCDF_LIBRARY NETCDF_INCLUDE_DIR)

IF(NETCDF_FOUND)
  MARK_AS_ADVANCED(
    NETCDF_LIBRARY
    NETCDF_INCLUDE_DIR
  )
ENDIF()

