## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2014 by the deal.II authors
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
# Try to find the MATHJAX library
#
# This module exports
#
#   MATHJAX_PATH
#

# SET(MATHJAX_DIR "" CACHE PATH "An optional hint to a MATHJAX installation")
# SET_IF_EMPTY(MATHJAX_DIR "$ENV{MATHJAX_DIR}")

# DEAL_II_FIND_LIBRARY(MATHJAX_LIBRARY
#   NAMES mathjax
#   HINTS ${MATHJAX_DIR}
#   )

# DEAL_II_FIND_PATH(MATHJAX_DIR
#   NAMES MathJax.js
#   PATHS "$ENV{MATHJAX_DIR}" "/usr/share/javascript/mathjax/"
#   DOC "Path to local MathJax.js"
#   )

DEAL_II_FIND_PATH(MATHJAX_PATH MathJax.js
  PATHS ${MATHJAX_ROOT}
  $ENV{MATHJAX_ROOT}
  "${CMAKE_PREFIX_PATH}/share/javascript/mathjax"
  "${CMAKE_INSTALL_DATADIR}/javascript/mathjax"
  "/usr/share/javascript/mathjax"
  DOC "Root path of MathJax.js"
  )

DEAL_II_PACKAGE_HANDLE(MATHJAX
  # LIBRARIES REQUIRED MATHJAX_LIBRARY
  REQUIRED_VARS  MATHJAX_PATH
  # CLEAR MATHJAX_DIR
  )
