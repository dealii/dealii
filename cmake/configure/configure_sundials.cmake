## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
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
# Configuration for the SUNDIALS library:
#

MACRO(FEATURE_SUNDIALS_CONFIGURE_EXTERNAL)
  SET(DEAL_II_SUNDIALS_WITH_IDAS ${SUNDIALS_WITH_IDAS})
ENDMACRO()

CONFIGURE_FEATURE(SUNDIALS)
