## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Configuration for the Mesquite library:
#

MACRO(FEATURE_MESQUITE_CONFIGURE_BUNDLED)
  SET(MESQUITE_BUNDLED_INCLUDE_DIRS
    ${MESQUITE_FOLDER}/src/include
    )
ENDMACRO()

MACRO(FEATURE_MESQUITE_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find Mesquite and supporting libraries!\n"
    "Please ensure that the libraries are installed on your computer.\n"
    "If the libraries are not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ MESQUITE_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DMESQUITE_DIR=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    "Alternatively you may choose to compile the bundled libraries\n"
    "by setting DEAL_II_ALLOW_BUNDLED=ON or DEAL_II_FORCE_BUNDLED_MESQUITE=ON.\n\n"
    )
ENDMACRO()

CONFIGURE_FEATURE(MESQUITE)