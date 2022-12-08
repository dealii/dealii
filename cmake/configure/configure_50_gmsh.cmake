## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2021 by the deal.II authors
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
# Configuration for the gmsh executable:
#

macro(feature_gmsh_configure_external)
  set(DEAL_II_GMSH_WITH_API ${GMSH_WITH_API})
endmacro()

configure_feature(GMSH)
