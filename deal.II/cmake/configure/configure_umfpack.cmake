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
# Configuration for the umfpack library:
#

SET(FEATURE_UMFPACK_DEPENDS DEAL_II_WITH_LAPACK)

MACRO(FEATURE_UMFPACK_CONFIGURE_BUNDLED)
  INCLUDE_DIRECTORIES(
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
