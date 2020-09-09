// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_simplex_grid_generator_h
#define dealii_simplex_grid_generator_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SIMPLEX_SUPPORT

namespace GridGenerator
{
  /**
   * Create a coordinate-parallel brick from the two diagonally opposite
   * corner points @p p1 and @p p2. The number of vertices in coordinate
   * direction @p i is given by <tt>repetitions[i]+1</tt>.
   *
   * @note This function connects internally 4/8 vertices to
   *   quadrilateral/hexahedral cells and subdivides these into 2/5
   * triangular/tetrahedral cells.
   *
   * @note Currently, this function only works for `dimDEAL_II_EQUALS spacedim`.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_rectangle_with_simplices(
    Triangulation<dim, spacedim> &   tria,
    const std::vector<unsigned int> &repetitions,
    const Point<dim> &               p1,
    const Point<dim> &               p2,
    const bool                       colorize = false);

  /**
   * Initialize the given triangulation with a hypercube (square in 2D and
   * cube in 3D) consisting of @p repetitions cells in each direction.
   * The hypercube volume is the tensor product interval
   * $[left,right]^{\text{dim}}$ in the present number of dimensions, where
   * the limits are given as arguments. They default to zero and unity, then
   * producing the unit hypercube.
   *
   * @note This function connects internally 4/8 vertices to
   * quadrilateral/hexahedral cells and subdivides these into 2/5
   * triangular/tetrahedral cells.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_cube_with_simplices(Triangulation<dim, spacedim> &tria,
                                       const unsigned int repetitions,
                                       const double       p1       = 0.0,
                                       const double       p2       = 1.0,
                                       const bool         colorize = false);

} // namespace GridGenerator

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
