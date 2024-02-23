// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_grid_tria_info_cache_update_flags_h
#define dealii_grid_tria_info_cache_update_flags_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  /**
   * The enum type given to the Cache class to select what
   * information to update.
   *
   * You can select more than one flag by concatenation using the bitwise or
   * <code>operator|(CacheUpdateFlags,CacheUpdateFlags)</code>.
   */
  enum CacheUpdateFlags
  {
    /**
     * Update Nothing.
     */
    update_nothing = 0x000,

    /**
     * Update vertex_to_cell_map, as returned by
     * GridTools::vertex_to_cell_map().
     */
    update_vertex_to_cell_map = 0x001,

    /**
     * Update vertex_to_cell_centers_directions, as returned by
     * GridTools::vertex_to_cell_centers_directions()
     */
    update_vertex_to_cell_centers_directions =
      update_vertex_to_cell_map | 0x002,

    /**
     * Update a mapping of used vertices.
     */
    update_used_vertices = 0x008,

    /**
     * Update an RTree of the used vertices.
     */
    update_used_vertices_rtree = 0x010,

    /**
     * Update an RTree of the cell bounding boxes.
     */
    update_cell_bounding_boxes_rtree = 0x020,

    /**
     * Update the covering rtree object, initialized with pairs
     * of a bounding box and an unsigned int. The bounding
     * boxes are used to describe approximately which portion
     * of the mesh contains locally owned cells by the
     * process of rank the second element of the pair.
     */
    update_covering_rtree = 0x040,

    /**
     * Update an RTree of locally owned cell bounding boxes.
     */
    update_locally_owned_cell_bounding_boxes_rtree = 0x080,

    /**
     * Update vertex to neighbor subdomain
     */
    update_vertex_to_neighbor_subdomain = 0x100,

    /**
     * Update the information about which subdomains are connected to each
     * vertex.
     */
    update_vertex_with_ghost_neighbors = 0x200,

    /**
     * Update all objects.
     */
    update_all = 0xFFF,
  };


  /**
   * Output operator which outputs assemble flags as a set of or'd text values.
   *
   * @ref CacheUpdateFlags
   */
  template <typename StreamType>
  inline StreamType &
  operator<<(StreamType &s, const CacheUpdateFlags u)
  {
    s << " CacheUpdateFlags";
    if (u & update_vertex_to_cell_map)
      s << "|vertex_to_cell_map";
    if (u & update_vertex_to_cell_centers_directions)
      s << "|vertex_to_cells_centers_directions";
    if (u & update_covering_rtree)
      s << "|covering_rtree";
    return s;
  }


  /**
   * Global operator which returns an object in which all bits are set which are
   * either set in the first or the second argument. This operator exists since
   * if it did not then the result of the bit-or <tt>operator |</tt> would be an
   * integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type CacheUpdateFlags.
   *
   * @ref CacheUpdateFlags
   */
  inline CacheUpdateFlags
  operator|(const CacheUpdateFlags f1, const CacheUpdateFlags f2)
  {
    return static_cast<CacheUpdateFlags>(static_cast<unsigned int>(f1) |
                                         static_cast<unsigned int>(f2));
  }

  /**
   * Global operator which returns an object in which all bits are set which are
   * not set in the argument. This operator exists since
   * if it did not then the result of the bit-negation <tt>operator ~</tt> would
   * be an integer which would in turn trigger a compiler warning when we tried
   * to assign it to an object of type CacheUpdateFlags.
   *
   * @ref CacheUpdateFlags
   */
  inline CacheUpdateFlags
  operator~(const CacheUpdateFlags f1)
  {
    return static_cast<CacheUpdateFlags>(static_cast<unsigned int>(f1) ^
                                         static_cast<unsigned int>(update_all));
  }



  /**
   * Global operator which sets the bits from the second argument also in the
   * first one.
   *
   * @ref CacheUpdateFlags
   */
  inline CacheUpdateFlags &
  operator|=(CacheUpdateFlags &f1, const CacheUpdateFlags f2)
  {
    f1 = f1 | f2;
    return f1;
  }


  /**
   * Global operator which returns an object in which all bits are set which are
   * set in the first as well as the second argument. This operator exists since
   * if it did not then the result of the bit-and <tt>operator &</tt> would be
   * an integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type CacheUpdateFlags.
   *
   * @ref CacheUpdateFlags
   */
  inline CacheUpdateFlags
  operator&(const CacheUpdateFlags f1, const CacheUpdateFlags f2)
  {
    return static_cast<CacheUpdateFlags>(static_cast<unsigned int>(f1) &
                                         static_cast<unsigned int>(f2));
  }


  /**
   * Global operator which clears all the bits in the first argument if they are
   * not also set in the second argument.
   *
   * @ref CacheUpdateFlags
   */
  inline CacheUpdateFlags &
  operator&=(CacheUpdateFlags &f1, const CacheUpdateFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }

} // namespace GridTools
DEAL_II_NAMESPACE_CLOSE

#endif
