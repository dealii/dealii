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
# Configuration for the MUMPS library:
#

OPTION(DEAL_II_WITH_MUMPS
  "Build deal.II with support for MUMPS."
  OFF)

SET(FEATURE_MUMPS_DEPENDS DEAL_II_WITH_MPI)


MACRO(FEATURE_MUMPS_FIND_EXTERNAL var)
  FIND_PACKAGE(MUMPS)

  IF(MUMPS_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MUMPS_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${MUMPS_INCLUDE_DIRS})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${MUMPS_LIBRARIES})
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${MUMPS_LINKER_FLAGS}")

  SET(DEAL_II_USE_MUMPS TRUE)

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(MUMPS)

