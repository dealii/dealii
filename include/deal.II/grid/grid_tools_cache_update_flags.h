// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

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
   *
   * @author Luca Heltai, 2017.
   */
  enum CacheUpdateFlags
  {
    /**
     * Update Nothing.
     */
    update_nothing = 0x00,

    /**
     * Update vertex_to_cell_map, as returned by GridTools::vertex_to_cell_map().
     */
    update_vertex_to_cell_map = 0x01,

    /**
     * Update vertex_to_cell_centers_directions, as returned by
     * GridTools::vertex_to_cell_centers_directions()
     */
    update_vertex_to_cell_centers_directions = update_vertex_to_cell_map | 0x0002,

    /**
     * Update a KDTree object, initialized with the vertices of the Triangulation.
     */
    update_vertex_kdtree = 0x04,

    /**
     * Update a mapping of used vertices.
     */
    update_used_vertices = 0x08,

    /**
     * Update the global bounding boxes describing the locally owned part, for each
     * process, of the mesh
     */
    update_global_bounding_boxes = 0x10,

    /**
     * Update all objects.
     */
    update_all = 0xFF,
  };


  /**
   * Output operator which outputs assemble flags as a set of or'd text values.
   *
   * @ref CacheUpdateFlags
   */
  template <class StreamType>
  inline
  StreamType &operator << (StreamType &s,
                           const CacheUpdateFlags u)
  {
    s << " CacheUpdateFlags";
    if (u & update_vertex_to_cell_map)                 s << "|vertex_to_cell_map";
    if (u & update_vertex_to_cell_centers_directions)  s << "|vertex_to_cells_centers_directions";
#ifdef DEAL_II_WITH_NANOFLANN
    if (u & update_vertex_kdtree)                      s << "|vertex_kdtree";
#endif
#ifdef DEAL_II_WITH_MPI
    if (u & update_global_bounding_boxes)              s << "|global_bounding_boxes";
#endif
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
  inline
  CacheUpdateFlags
  operator | (const CacheUpdateFlags f1,
              const CacheUpdateFlags f2)
  {
    return static_cast<CacheUpdateFlags> (
             static_cast<unsigned int> (f1) |
             static_cast<unsigned int> (f2));
  }

  /**
   * Global operator which returns an object in which all bits are set which are
   * not set in the argument. This operator exists since
   * if it did not then the result of the bit-negation <tt>operator ~</tt> would be an
   * integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type CacheUpdateFlags.
   *
   * @ref CacheUpdateFlags
   */
  inline
  CacheUpdateFlags
  operator ~ (const CacheUpdateFlags f1)
  {
    return static_cast<CacheUpdateFlags> (
             ~static_cast<unsigned int> (f1)
           );
  }



  /**
   * Global operator which sets the bits from the second argument also in the
   * first one.
   *
   * @ref CacheUpdateFlags
   */
  inline
  CacheUpdateFlags &
  operator |= (CacheUpdateFlags &f1,
               const CacheUpdateFlags f2)
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
  inline
  CacheUpdateFlags
  operator & (const CacheUpdateFlags f1,
              const CacheUpdateFlags f2)
  {
    return static_cast<CacheUpdateFlags> (
             static_cast<unsigned int> (f1) &
             static_cast<unsigned int> (f2));
  }


  /**
   * Global operator which clears all the bits in the first argument if they are
   * not also set in the second argument.
   *
   * @ref CacheUpdateFlags
   */
  inline
  CacheUpdateFlags &
  operator &= (CacheUpdateFlags &f1,
               const CacheUpdateFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }

}
DEAL_II_NAMESPACE_CLOSE

#endif
