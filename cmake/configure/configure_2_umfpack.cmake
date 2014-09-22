## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Configuration for the umfpack library:
#

SET(FEATURE_UMFPACK_DEPENDS LAPACK)

MACRO(FEATURE_UMFPACK_CONFIGURE_BUNDLED)
  SET(UMFPACK_BUNDLED_INCLUDE_DIRS
    ${UMFPACK_FOLDER}/UMFPACK/Include
    ${UMFPACK_FOLDER}/AMD/Include
    )
ENDMACRO()

MACRO(FEATURE_UMFPACK_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find umfpack and supporting libraries!\n"
    "Please ensure that the libraries are installed on your computer.\n"
    "If the libraries are not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ UMFPACK_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DUMFPACK_DIR=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    "Relevant hints for UMFPACK are SUITESPARSE_DIR, UMFPACK_DIR\n"
    "(AMD_DIR, CHOLMOD_DIR, COLAMD_DIR, SUITESPARSECONFIG_DIR.)\n"
    "Alternatively you may choose to compile the bundled libraries\n"
    "by setting DEAL_II_ALLOW_BUNDLED=ON or DEAL_II_FORCE_BUNDLED_UMFPACK=ON.\n"
    "(BLAS and LAPACK have to be installed for bundled UMFPACK to be available)\n\n"
    )
ENDMACRO()

CONFIGURE_FEATURE(UMFPACK)
