## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Keep track of all of the header files in this directory
#
set(_header_files
    mg_base.h
    mg_base.templates.h
    mg_block_smoother.h
    mg_coarse.h
    mg_constrained_dofs.h
    mg_matrix.h
    mg_smoother.h
    mg_tools.h
    mg_transfer_block.h
    mg_transfer_block.templates.h
    mg_transfer_component.h
    mg_transfer_component.templates.h
    mg_transfer_global_coarsening.h
    mg_transfer_global_coarsening.templates.h
    mg_transfer.h
    mg_transfer_internal.h
    mg_transfer_matrix_free.h
    mg_transfer_matrix_free.templates.h
    mg_transfer.templates.h
    multigrid.h
    multigrid.templates.h
    sparse_matrix_collection.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
