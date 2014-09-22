## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# avoid some annoying warnings...
#
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wall")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-switch)
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-nested-anon-types)

INCLUDE_DIRECTORIES(
  include/
)

DEAL_II_ADD_LIBRARY(obj_muparser OBJECT
src/muParserBase.cpp
src/muParserBytecode.cpp
src/muParserCallback.cpp
src/muParser.cpp
src/muParserError.cpp
src/muParserInt.cpp
src/muParserTest.cpp
src/muParserTokenReader.cpp
  )

#IF(DEAL_II_WITH_THREADS)
#  DEAL_II_ADD_DEFINITIONS(obj_muparser "FP_USE_THREAD_SAFE_EVAL")
#ENDIF()
