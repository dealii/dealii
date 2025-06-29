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
  CGAL::Polygon_2<KernelType>
  dealii_tria_to_cgal_polygon(const Triangulation<2, 2> &tria,
                              const Mapping<2, 2>       &mapping)
  {
    // This map holds the two vertex indices of each face.
    // Counterclockwise first vertex index on first position,
    // counterclockwise second vertex index on second position.
    std::map<unsigned int, unsigned int> face_vertex_indices;
    std::map<unsigned int, Point<2>>     vertex_to_point;

    // Iterate over all active cells at the boundary
    for (const auto &cell : tria.active_cell_iterators())
      {
        for (const unsigned int f : cell->face_indices())
          {
            if (cell->face(f)->at_boundary())
              {
                // get mapped vertices of the cell
                const auto         v_mapped = mapping.get_vertices(cell);
                const unsigned int v0       = cell->face(f)->vertex_index(0);
                const unsigned int v1       = cell->face(f)->vertex_index(1);

                if (cell->reference_cell() == ReferenceCells::Triangle)
                  {
                    // add indices and first mapped vertex of the face
                    vertex_to_point[v0]     = v_mapped[f];
                    face_vertex_indices[v0] = v1;
                  }
                else if (cell->reference_cell() ==
                         ReferenceCells::Quadrilateral)
                  {
                    // Ensure that vertex indices of the face are in
                    // counterclockwise order inserted in the map.
                    if (f == 0 || f == 3)
                      {
                        // add indices and first mapped vertex of the face
                        vertex_to_point[v1] =
                          v_mapped[GeometryInfo<2>::face_to_cell_vertices(f,
                                                                          1)];
                        face_vertex_indices[v1] = v0;
                      }
                    else
                      {
                        // add indices and first mapped vertex of the face
                        vertex_to_point[v0] =
                          v_mapped[GeometryInfo<2>::face_to_cell_vertices(f,
                                                                          0)];
                        face_vertex_indices[v0] = v1;
                      }
                  }
                else
                  {
                    DEAL_II_ASSERT_UNREACHABLE();
                  }
              }
          }
      }

    CGAL::Polygon_2<KernelType> polygon;

    // Vertex to start counterclockwise insertion
    const unsigned int start_index   = face_vertex_indices.begin()->first;
    unsigned int       current_index = start_index;

    // As long as still entries in the map, use last vertex index to
    // find next vertex index in counterclockwise order
    while (face_vertex_indices.size() > 0)
      {
        const auto vertex_it = vertex_to_point.find(current_index);
        Assert(vertex_it != vertex_to_point.end(),
               ExcMessage("This should not occur, please report bug"));

        polygon.push_back(
          dealii_point_to_cgal_point<CGAL::Point_2<KernelType>, 2>(
            vertex_it->second));
        vertex_to_point.erase(vertex_it);

        const auto it = face_vertex_indices.find(current_index);
        // If the boundary is one closed loop, the next vertex index
        // must exist as key until the map is empty.
        Assert(it != face_vertex_indices.end(),
               ExcMessage("Triangulation might contain holes"));

        current_index = it->second;
        face_vertex_indices.erase(it);

        // Ensure that last vertex index is the start index.
        // This condition can be used to extend the code for
        // trianglations with holes
        Assert(
          !(face_vertex_indices.size() == 0 && current_index != start_index),
          ExcMessage(
            "This should not occur, reasons might be a non closed boundary or a bug in this function"));
      }

    return polygon;
  }



  template <typename KernelType>
  std::vector<CGAL::Polygon_with_holes_2<KernelType>>
  compute_boolean_operation(const CGAL::Polygon_2<KernelType> &polygon_1,
                            const CGAL::Polygon_2<KernelType> &polygon_2,
                            const BooleanOperation &boolean_operation)
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
