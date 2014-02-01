## ---------------------------------------------------------------------
## $Id$
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
# Configuration for functionparser
#

MACRO(FEATURE_FUNCTIONPARSER_FIND_EXTERNAL var)
  MESSAGE(STATUS
    "No module available for finding functionparser externally."
    )
ENDMACRO()


MACRO(FEATURE_FUNCTIONPARSER_CONFIGURE_BUNDLED)
  SET(FUNCTIONPARSER_BUNDLED_INCLUDE_DIRS ${FUNCTIONPARSER_FOLDER})
ENDMACRO()


MACRO(FEATURE_FUNCTIONPARSER_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "No module available for finding functionparser externally.\n"
    "Disable DEAL_II_WITH_FUNCTIONPARSER, or enable DEAL_II_ALLOW_BUNDLED.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(FUNCTIONPARSER)
