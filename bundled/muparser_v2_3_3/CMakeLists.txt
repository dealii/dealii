## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2023 by the deal.II authors
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
strip_flag(DEAL_II_WARNING_FLAGS "-Wall")
strip_flag(DEAL_II_WARNING_FLAGS "-pedantic")
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-switch)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-nested-anon-types)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-implicit-fallthrough)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-cast-function-type)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-float-conversion)

include_directories(
  include/
)

define_object_library(bundled_muparser OBJECT
  src/muParserBase.cpp
  src/muParserBytecode.cpp
  src/muParserCallback.cpp
  src/muParser.cpp
  src/muParserError.cpp
  src/muParserInt.cpp
  src/muParserTest.cpp
  src/muParserTokenReader.cpp
  )
