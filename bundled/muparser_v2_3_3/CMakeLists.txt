## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
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
# avoid some annoying warnings...
#
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wall")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-switch)
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-nested-anon-types)
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-implicit-fallthrough)
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-cast-function-type)
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-float-conversion)

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

# We do not support (yet) compiling DLLs with MSVC. By default, muParser will
# try to set itself up correctly with __declspec(dllexport). Get around this by
# instructing it to use static linkage, which will define (see muParserFixes.h)
# AP_EXPORT_CXX to nothing.
DEAL_II_ADD_DEFINITIONS(obj_muparser "MUPARSER_STATIC")
