## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Configuration for the Ginkgo library:
#

MACRO(FEATURE_GINKGO_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find Ginkgo and supporting libraries!\n"
    "Please ensure that the libraries are installed on your computer.\n"
    "If the libraries are not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ GINKGO_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DGINKGO_DIR=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    "Relevant hints for GINKGO are GINKGO_DIR.\n"
    )
ENDMACRO()

MACRO(FEATURE_GINKGO_CONFIGURE_EXTERNAL)
  SET(DEAL_II_GINKGO_BUILT_REFERENCE ${GINKGO_BUILT_REFERENCE})
  SET(DEAL_II_GINKGO_BUILT_OPENMP ${GINKGO_BUILT_OMP})
  SET(DEAL_II_GINKGO_BUILT_CUDA ${GINKGO_BUILT_CUDA})
ENDMACRO()

CONFIGURE_FEATURE(GINKGO)
