#####
##
## Copyright (C) 2012 by the deal.II authors
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
# A small macro used in the platform checks to remove the right most flag in
# CMAKE_REQUIRED_FLAGS
#
# We assume that the flags in CMAKE_REQUIRED_FLAGS are space separated
#
# Usage:
#     POP_TEST_FLAG()
#

MACRO(POP_TEST_FLAG)
  SET(CMAKE_REQUIRED_FLAGS " ${CMAKE_REQUIRED_FLAGS}")
  STRING(REGEX REPLACE " -[^ ]+$" ""
    CMAKE_REQUIRED_FLAGS
    "${CMAKE_REQUIRED_FLAGS}"
    )
  STRING(STRIP "${CMAKE_REQUIRED_FLAGS}" CMAKE_REQUIRED_FLAGS)
ENDMACRO()

