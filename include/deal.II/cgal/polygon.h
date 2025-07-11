// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cgal_polygon_h
#define dealii_cgal_polygon_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CGAL
#  include <deal.II/cgal/point_conversion.h>
#  include <deal.II/cgal/utilities.h>

#  include <CGAL/version.h>
#  if CGAL_VERSION_MAJOR >= 6
#    include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#  endif
#  include <CGAL/Polygon_2.h>
#  include <CGAL/Polygon_with_holes_2.h>


DEAL_II_NAMESPACE_OPEN

namespace CGALWrappers
{
  /**
   * Build a CGAL::Polygon_2 from a deal.II cell.
   *
   * The class Polygon_2 is a wrapper around a container of points that can
   * be used to represent polygons.
   * The points are added in counterclockwise order to a Polygon_2
   *
   * More information on this class is available at
   * https://doc.cgal.org/latest/Polygon/index.html
   *
   * The functions are for 2D triangulations in 2D space.
   * Projecting 3D points is possible with CGAL but not implemented.
   *
   * The generated boundary representation can be used to perform
   * regularized boolean set-operations.
   *
   * @param[in] cell The input deal.II cell iterator
   * @param[in] mapping The mapping used to map the vertices of the cell
   */
  template <typename KernelType>
  CGAL::Polygon_2<KernelType>
  dealii_cell_to_cgal_polygon(
    const typename Triangulation<2, 2>::cell_iterator &cell,
    const Mapping<2, 2>                               &mapping);



  /**
   * Convert a deal.II triangulation to a CGAL::Polygon_2.
   *
   * Triangulations that have holes are not supported.
   *
   * @param[in] tria The input deal.II triangulation
   * @param[in] mapping The mapping used to map the vertices of cells
   */
  template <typename KernelType>
  CGAL::Polygon_with_holes_2<KernelType>
  dealii_tria_to_cgal_polygon(const Triangulation<2, 2> &tria,
                              const Mapping<2, 2>       &mapping);



  /**
   * Constructs a Polygon_with_holes_2 from the input Polygon_2.
   * Further polygons for holes are optional.
   *
   * Polygon_with_holes_2 has function .outer_boundary() as well as
   * .holes() or .holes_begin() and .holes_end()
   *
   * @param[in] boundary_outside Polygon for the outer boundary
   * @param[in] boundary_holes Polygons for the holes
   */
  template <typename KernelType>
  CGAL::Polygon_with_holes_2<KernelType>
  polygon_to_polygon_with_holes(
    const CGAL::Polygon_2<KernelType>              &boundary_outside,
    const std::vector<CGAL::Polygon_2<KernelType>> &boundary_holes = {});



  /**
   * Perform a Regularized Boolean set-operation on two CGAL::Polygon_2.
   *
   * More information on this class is available at
   * https://doc.cgal.org/latest/Boolean_set_operations_2/index.html
   *
   * The output is generally a vector of CGAL::Polygon_2_with_holes.
   *
   * For the union operation the vector will always have length one.
   *
   * For the difference operation the second polygon is subtracted
   * from the first one.
   *
   * Corefinement is not supported as boolean operation.
   *
   * @param[in] polygon_1 The first input CGAL::Polygon_2
   * @param[in] polygon_2 The second input CGAL::Polygon_2
   * @param[in] boolean_operation The input BooleanOperation
   */
  template <typename KernelType>
  std::vector<CGAL::Polygon_with_holes_2<KernelType>>
  compute_boolean_operation(
    const CGAL::Polygon_with_holes_2<KernelType> &polygon_1,
    const CGAL::Polygon_with_holes_2<KernelType> &polygon_2,
    const BooleanOperation                       &boolean_operation);
} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
