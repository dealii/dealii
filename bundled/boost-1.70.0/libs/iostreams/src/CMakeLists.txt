## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 - 2020 by the deal.II authors
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

set(src_boost_iostreams
    file_descriptor.cpp
    gzip.cpp
    mapped_file.cpp
    zlib.cpp
    )
enable_if_supported(DEAL_II_CXX_FLAGS -Wno-c11-extensions)
enable_if_supported(DEAL_II_CXX_FLAGS -Wno-deprecated-copy)
enable_if_supported(DEAL_II_CXX_FLAGS -Wno-uninitialized)

deal_ii_add_library(obj_boost_iostreams OBJECT ${src_boost_iostreams})
