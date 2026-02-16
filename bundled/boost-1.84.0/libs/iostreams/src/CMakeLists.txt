## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2013 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

set(src_boost_iostreams
    file_descriptor.cpp
    gzip.cpp
    mapped_file.cpp
    zlib.cpp
    )
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-c11-extensions)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-deprecated-copy)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-uninitialized)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-deprecated-builtins)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-deprecated-declarations)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
define_object_library(bundled_boost_iostreams OBJECT ${src_boost_iostreams})
