FIND_PACKAGE(Netcdf REQUIRED)

INCLUDE_DIRECTORIES(${Netcdf_INCLUDE_DIR})

LIST(APPEND deal_ii_external_libraries ${Netcdf_LIBRARY})
LIST(APPEND deal_ii_external_debug_libraries ${Netcdf_LIBRARY})

SET(HAVE_LIBNETCDF TRUE)
