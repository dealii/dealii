## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
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
# Configuration for the OpenCASCADE library:
#

MACRO(FEATURE_OPENCASCADE_CONFIGURE_EXTERNAL)
  #
  # Disable a bunch of warnings caused by OpenCascade headers:
  #
  ENABLE_IF_SUPPORTED(TRILINOS_CXX_FLAGS "-Wno-extra")
ENDMACRO()


CONFIGURE_FEATURE(OPENCASCADE)
