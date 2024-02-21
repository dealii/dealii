## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/)

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
