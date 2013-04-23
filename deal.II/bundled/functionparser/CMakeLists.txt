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
# Remove -Wall and -pedantic from CMAKE_CXX_FLAGS (in directory scope) to
# avoid some annoying warnings...
#
STRIP_FLAG(CMAKE_CXX_FLAGS "-Wall")
STRIP_FLAG(CMAKE_CXX_FLAGS "-pedantic")

DEAL_II_ADD_LIBRARY(obj_functionparser OBJECT
  fparser.cc
  fpoptimizer.cc
  )

