// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_cgal_intersections_h
#define dealii_cgal_intersections_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN
namespace CGALWrappers
{
  /**
   * Given two deal.II affine cells, compute their intersection and return a
   * vector of simplices, each one identified by an array of deal.II Points. All
   * the simplices together are a subdivision of the intersection. If cells are
   * non-affine, a geometrical error is introduced. If the
   * measure of one of the simplices is below a certain treshold which defaults
   * to 1e-9, then it is discarded. In case the two cells are disjoint, an empty
   * array is returned.
   *
   * @note If the two cells have different intrinsic dimensions, than the
   * first iterator is assumed to have with higher intrinsic dimension.
   * The same holds for the two Mapping objects.
   *
   *
   * @param cell0 Iterator to the first cell.
   * @param cell1 Iterator to the second cell.
   * @param mapping0 Mapping for the first cell.
   * @param mapping1 Mapping for the second cell.
   * @param discard_tolerance Treshold to decide whether or not a found simplex is included.
   * @param vertex_tolerance Tolerance to find an intersection. 
   * In some cases intersections should only be found if cells are intersecting up to machine 
   * precision. This is the default case.
   * In some cases however, we want to find an intersection with a certain tolerance. This is 
   * for example the case if we constructed a triangulation from two non-fitting 
   * triangulations to realize required jumps in element sizes. To find intersections up to a 
   * certain tolerance we slightly modify the vertex positions before asking CGAL about 
   * intersections. This is done by enforcing the same floating point representation for 
   * vertices that are equivalent up to a certain tolerance using @c Utilities::round().
   * @return Vector of arrays, where each array identify a simplex by its vertices.
   * @note Rounding with vertex_tolerance does not have any effect on the triangulation, it is 
   * only applied to a copy of vertices handed to CGAL.
   */
  template <int dim0, int dim1, int spacedim>
  std::vector<std::array<Point<spacedim>, dim1 + 1>>
  compute_intersection_of_cells(
    const typename Triangulation<dim0, spacedim>::cell_iterator &cell0,
    const typename Triangulation<dim1, spacedim>::cell_iterator &cell1,
    const Mapping<dim0, spacedim> &                              mapping0,
    const Mapping<dim1, spacedim> &                              mapping1,
    const double discard_tolerance = 1e-9,
    const double vertex_tolerance  = std::numeric_limits<double>::min());


  /**
   * Same function as above, but working directly with vertices.
   */
  template <int dim0, int dim1, int spacedim, int n_vertices0, int n_vertices1>
  std::vector<std::array<Point<spacedim>, dim1 + 1>>
  compute_intersection_of_cells(
    const std::array<Point<spacedim>, n_vertices0> &vertices0,
    const std::array<Point<spacedim>, n_vertices1> &vertices1,
    const double                                    discard_tolerance = 1e-9)
  {
    (void)vertices0;
    (void)vertices1;
    (void)discard_tolerance;
    Assert(false, ExcMessage("No explicit template instantiation available"));
    return {};
  }
} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
