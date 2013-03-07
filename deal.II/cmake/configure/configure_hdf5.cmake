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
# Configuration for the hdf5 library:
#

MACRO(FEATURE_HDF5_FIND_EXTERNAL var)
  FIND_PACKAGE(HDF5)

  IF(HDF5_FOUND)
    IF( (HDF5_WITH_MPI AND DEAL_II_WITH_MPI)
         OR
         (NOT HDF5_WITH_MPI AND NOT DEAL_II_WITH_MPI))
      SET(${var} TRUE)
    ELSE()
      SET(HDF5_ADDITIONAL_WARNING_STRING
        "Insufficient hdf5 installation found!\n"
        "hdf5 has to be configured with the same MPI configuration as deal.II, but found:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
        "  HDF5_WITH_MPI    = ${HDF5_WITH_MPI}\n"
        )
      MESSAGE(WARNING "\n" ${HDF5_ADDITIONAL_WARNING_STRING} "\n")
    ENDIF()
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(HDF5)
