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

#include <deal.II/base/config.h>

#include <deal.II/cgal/polygon.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h> //should be ok to delete is in grid tools

#ifdef DEAL_II_WITH_CGAL

#  include <CGAL/Boolean_set_operations_2.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
// Template implementations
namespace CGALWrappers
{
  template <typename KernelType>
  CGAL::Polygon_2<KernelType>
  dealii_cell_to_cgal_polygon(
    const typename Triangulation<2, 2>::cell_iterator &cell,
    const Mapping<2, 2>                               &mapping)
  {
    CGAL::Polygon_2<KernelType> polygon;
    const auto                 &vertices = mapping.get_vertices(cell);

    if (cell->reference_cell() == ReferenceCells::Triangle)
      {
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[0]));
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[1]));
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[2]));
      }
    else if (cell->reference_cell() == ReferenceCells::Quadrilateral)
      {
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[0]));
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[1]));
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[3]));
        polygon.push_back(
          CGALWrappers::dealii_point_to_cgal_point<CGAL::Point_2<KernelType>,
                                                   2>(vertices[2]));
      }
    else
      {
        DEAL_II_ASSERT_UNREACHABLE();
      }
    return polygon;
  }



  template <typename KernelType>
  CGAL::Polygon_with_holes_2<KernelType>
  dealii_tria_to_cgal_polygon(const Triangulation<2, 2> &tria,
                              const Mapping<2, 2>       &mapping)
  {
    // Outer boundary vertices counter clockwise
    // and hole vertices clockwise ordered
    auto boundaries =
      GridTools::extract_ordered_boundary_vertices(tria, mapping);

    CGAL::Polygon_2<KernelType>              outer_boundary;
    std::vector<CGAL::Polygon_2<KernelType>> holes;
    holes.reserve(boundaries.size() - 1);

    for (const auto &boundary : boundaries)
      {
        // Add vertices of current boundary to polygon
        CGAL::Polygon_2<KernelType> current_polygon;
        for (const auto &vertices : boundary)
          {
            current_polygon.push_back(
              dealii_point_to_cgal_point<CGAL::Point_2<KernelType>, 2>(
                vertices.second));
          }
        // Decide if outer boundary or hole
        if (current_polygon.is_counterclockwise_oriented())
          {
            outer_boundary = current_polygon;
          }
        else
          {
            holes.push_back(current_polygon);
          }
      }
    return CGAL::Polygon_with_holes_2<KernelType>(outer_boundary,
                                                  holes.begin(),
                                                  holes.end());
  }


  template <typename KernelType>
  CGAL::Polygon_with_holes_2<KernelType>
  polygon_to_polygon_with_holes(
    const CGAL::Polygon_2<KernelType>              &boundary_outside,
    const std::vector<CGAL::Polygon_2<KernelType>> &boundary_holes)
  {
    return CGAL::Polygon_with_holes_2<KernelType>(boundary_outside,
                                                  boundary_holes.begin(),
                                                  boundary_holes.end());
  }



  template <typename KernelType>
  std::vector<CGAL::Polygon_with_holes_2<KernelType>>
  compute_boolean_operation(
    const CGAL::Polygon_with_holes_2<KernelType> &polygon_1,
    const CGAL::Polygon_with_holes_2<KernelType> &polygon_2,
    const BooleanOperation                       &boolean_operation)
  {
    Assert(!(boolean_operation == BooleanOperation::compute_corefinement),
           ExcMessage("Corefinement has no usecase for 2D polygons"));

    std::vector<CGAL::Polygon_with_holes_2<KernelType>> polygon_out;

    if (boolean_operation == BooleanOperation::compute_intersection)
      {
        CGAL::intersection(polygon_1,
                           polygon_2,
                           std::back_inserter(polygon_out));
      }
    else if (boolean_operation == BooleanOperation::compute_difference)
      {
        CGAL::difference(polygon_1, polygon_2, std::back_inserter(polygon_out));
      }
    else if (boolean_operation == BooleanOperation::compute_union)
      {
        polygon_out.resize(1);
        CGAL::join(polygon_1, polygon_2, polygon_out[0]);
      }

    return polygon_out;
  }

// Explicit instantiations.
//
// We don't build the instantiations.inst file if deal.II isn't
// configured with CGAL, but doxygen doesn't know that and tries to
// find that file anyway for parsing -- which then of course it fails
// on. So exclude the following from doxygen consideration.
#    ifndef DOXYGEN
#      include "cgal/polygon.inst"
#    endif

} // namespace  CGALWrappers
#  endif


DEAL_II_NAMESPACE_CLOSE

#endif
