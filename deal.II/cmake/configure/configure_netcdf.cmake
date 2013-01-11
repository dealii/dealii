#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Configuration for the netcdf library:
#


MACRO(FEATURE_NETCDF_FIND_EXTERNAL var)
  FIND_PACKAGE(NETCDF)

  IF(NETCDF_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_NETCDF_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIRS})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${NETCDF_LIBRARIES})
  SET(HAVE_LIBNETCDF TRUE)

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(NETCDF)
