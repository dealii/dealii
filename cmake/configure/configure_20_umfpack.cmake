## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Configuration for the umfpack library:
#

set(FEATURE_UMFPACK_DEPENDS LAPACK)


macro(feature_umfpack_error_message)
  message(FATAL_ERROR "\n"
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
endmacro()


configure_feature(UMFPACK)
