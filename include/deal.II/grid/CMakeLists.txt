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
    cell_data.h
    cell_id.h
    cell_id_translator.h
    cell_status.h
    composition_manifold.h
    connectivity.h
    filtered_iterator.h
    grid_generator.h
    grid_in.h
    grid_out.h
    grid_refinement.h
    grid_tools_cache.h
    grid_tools_cache_update_flags.h
    grid_tools_geometry.h
    grid_tools.h
    grid_tools_topology.h
    intergrid_map.h
    magic_numbers.h
    manifold.h
    manifold_lib.h
    persistent_tria.h
    reference_cell.h
    tensor_product_manifold.h
    tria_accessor.h
    tria_accessor.templates.h
    tria_description.h
    tria_faces.h
    tria.h
    tria_iterator_base.h
    tria_iterator.h
    tria_iterator_selector.h
    tria_iterator.templates.h
    tria_levels.h
    tria_objects.h
    tria_objects_orientations.h
    tria_orientation.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
