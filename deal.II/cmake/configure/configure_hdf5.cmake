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


MACRO(FEATURE_HDF5_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${HDF5_LIBRARIES})

  SET(${var} TRUE)
ENDMACRO()


MACRO(FEATURE_HDF5_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find the hdf5 library!\n"
    "Please ensure that the hdf5 library is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ HDF5_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DHDF5_DIR=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(HDF5)
