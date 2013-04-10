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
# Configuration for the MUMPS library:
#

SET(FEATURE_MUMPS_DEPENDS DEAL_II_WITH_MPI DEAL_II_WITH_LAPACK)

MACRO(FEATURE_MUMPS_CONFIGURE_EXTERNAL)
  INCLUDE_DIRECTORIES(${MUMPS_INCLUDE_DIRS})

  # The user has to know the location of the MUMPS headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${MUMPS_INCLUDE_DIRS})

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${MUMPS_LIBRARIES})
ENDMACRO()

CONFIGURE_FEATURE(MUMPS)
