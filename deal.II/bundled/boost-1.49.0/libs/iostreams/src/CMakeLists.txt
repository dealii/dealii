## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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

INCLUDE_DIRECTORIES(${BOOST_BUNDLED_INCLUDE_DIRS})

SET(src_boost_iostreams
    file_descriptor.cpp
    gzip.cpp
    mapped_file.cpp
    zlib.cpp
    )

IF(DEALII_WITH_BZIP2)
  LIST(APPEND src_boost_iostreams bzip2.cpp)
ELSE()
  MESSAGE(STATUS "BOOST::Iostreams will not support bz2'ed streams because libbz2 or its header files could not be found")
ENDIF()

DEAL_II_ADD_LIBRARY(obj_boost_iostreams OBJECT ${src_boost_iostreams})
