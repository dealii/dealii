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

#include <deal.II/grid/tria.h>

#ifdef DEAL_II_WITH_CGAL

#  include <CGAL/Boolean_set_operations_2.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
// Template implementations
namespace CGALWrappers
{
  template <typename KernelType>
  void
  dealii_cell_to_cgal_polygon(
    const typename Triangulation<2, 2>::cell_iterator &cell,
    const Mapping<2, 2>                               &mapping,
    CGAL::Polygon_2<KernelType>                       &polygon)
  {
    const auto &vertices = mapping.get_vertices(cell);
    polygon.clear();
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
  }



  template <typename KernelType>
  void
  dealii_tria_to_cgal_polygon(const Triangulation<2, 2>   &tria,
                              CGAL::Polygon_2<KernelType> &polygon)
  {
    std::map<unsigned int, unsigned int> face_vertex_indices;

    for (const auto &cell : tria.active_cell_iterators())
      {
        for (const unsigned int i : cell->face_indices())
          {
            const typename Triangulation<2, 2>::face_iterator &face =
              cell->face(i);
            if (face->at_boundary())
              {
                if (cell->reference_cell() == ReferenceCells::Quadrilateral &&
                    (i == 0 || i == 3))
                  {
                    face_vertex_indices[face->vertex_index(1)] =
                      face->vertex_index(0);
                  }
                else
                  {
                    face_vertex_indices[face->vertex_index(0)] =
                      face->vertex_index(1);
                  }
              }
          }
      }

    const auto &vertices = tria.get_vertices();
    polygon.clear();

    unsigned int current_index = face_vertex_indices.begin()->first;

    for (size_t i = face_vertex_indices.size(); i > 0; --i)
      {
        polygon.push_back(
          dealii_point_to_cgal_point<CGAL::Point_2<KernelType>, 2>(
            vertices[current_index]));
        auto it       = face_vertex_indices.find(current_index);
        current_index = it->second;
        face_vertex_indices.erase(it);
      }
  }



  template <typename KernelType>
  void
  compute_boolean_operation(
    const CGAL::Polygon_2<KernelType>                   &polygon_1,
    const CGAL::Polygon_2<KernelType>                   &polygon_2,
    const BooleanOperation                              &boolean_operation,
    std::vector<CGAL::Polygon_with_holes_2<KernelType>> &polygon_out)
  {
    Assert(!(boolean_operation == BooleanOperation::compute_corefinement),
           ExcMessage("Corefinement has no usecase for 2D polygons"));

    polygon_out.clear();

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
