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

#ifndef dealii_grid_tria_info_cache_flags_h
#define dealii_grid_tria_info_cache_flags_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * The enum type given to the TriangulationInfoCache class to select what
 * information to cache.
 *
 * You can select more than one flag by concatenation using the bitwise or
 * <code>operator|(TriangulationInfoCacheFlags,TriangulationInfoCacheFlags)</code>.
 *
 * @author Luca Heltai, 2017.
 */
enum TriangulationInfoCacheFlags
{
  /**
   * Cache Nothing.
   */
  cache_nothing = 0,

  /**
   * Cache vertex_to_cell_map, as returned by GridTools::vertex_to_cell_map().
   */
  cache_vertex_to_cell_map = 0x0001,
};


/**
 * Output operator which outputs assemble flags as a set of or'd text values.
 *
 * @ref TriangulationInfoCacheFlags
 */
template <class StreamType>
inline
StreamType &operator << (StreamType &s, TriangulationInfoCacheFlags u)
{
  s << " TriangulationInfoCacheFlags";
  if (u & cache_vertex_to_cell_map) s << "|vertex_to_cell_map"        ;
  return s;
}


/**
 * Global operator which returns an object in which all bits are set which are
 * either set in the first or the second argument. This operator exists since
 * if it did not then the result of the bit-or <tt>operator |</tt> would be an
 * integer which would in turn trigger a compiler warning when we tried to
 * assign it to an object of type TriangulationInfoCacheFlags.
 *
 * @ref TriangulationInfoCacheFlags
 */
inline
TriangulationInfoCacheFlags
operator | (TriangulationInfoCacheFlags f1, TriangulationInfoCacheFlags f2)
{
  return static_cast<TriangulationInfoCacheFlags> (
           static_cast<unsigned int> (f1) |
           static_cast<unsigned int> (f2));
}




/**
 * Global operator which sets the bits from the second argument also in the
 * first one.
 *
 * @ref TriangulationInfoCacheFlags
 */
inline
TriangulationInfoCacheFlags &
operator |= (TriangulationInfoCacheFlags &f1, TriangulationInfoCacheFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set which are
 * set in the first as well as the second argument. This operator exists since
 * if it did not then the result of the bit-and <tt>operator &</tt> would be
 * an integer which would in turn trigger a compiler warning when we tried to
 * assign it to an object of type TriangulationInfoCacheFlags.
 *
 * @ref TriangulationInfoCacheFlags
 */
inline
TriangulationInfoCacheFlags
operator & (TriangulationInfoCacheFlags f1, TriangulationInfoCacheFlags f2)
{
  return static_cast<TriangulationInfoCacheFlags> (
           static_cast<unsigned int> (f1) &
           static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if they are
 * not also set in the second argument.
 *
 * @ref TriangulationInfoCacheFlags
 */
inline
TriangulationInfoCacheFlags &
operator &= (TriangulationInfoCacheFlags &f1, TriangulationInfoCacheFlags f2)
{
  f1 = f1 & f2;
  return f1;
}

DEAL_II_NAMESPACE_CLOSE

#endif
