#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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
# Configuration for the p4est and sc libraries:
#

SET(FEATURE_P4EST_DEPENDS DEAL_II_WITH_MPI)

MACRO(FEATURE_P4EST_FIND_EXTERNAL var)
  FIND_PACKAGE(P4EST)

  IF(P4EST_FOUND)
    #
    # Check whether p4est supports mpi:
    #
    IF(NOT P4EST_WITH_MPI)
      MESSAGE(STATUS "Insufficient p4est installation found: "
        "p4est has to be configured with MPI enabled."
        )
      SET(HDF5_ADDITIONAL_ERROR_STRING
        "Insufficient p4est installation found!\n"
        "p4est has to be configured with MPI enabled.\n"
        )

      UNSET(P4EST_LIBRARY CACHE)
      UNSET(P4EST_INCLUDE_DIR CACHE)
      UNSET(SC_LIBRARY CACHE)
      UNSET(SC_INCLUDE_DIR CACHE)
      SET(P4EST_DIR "" CACHE PATH
        "An optional hint to a p4est installation/directory"
        )
    ELSE()
      SET(${var} TRUE)
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO(FEATURE_P4EST_CONFIGURE_EXTERNAL)
  INCLUDE_DIRECTORIES(${P4EST_INCLUDE_DIRS})

  # The user has to know the location of the P4EST headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${P4EST_INCLUDE_DIRS})

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${P4EST_LIBRARIES})
ENDMACRO()

CONFIGURE_FEATURE(P4EST)
