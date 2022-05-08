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

#ifndef dealii_cgal_triangulation_h
#define dealii_cgal_triangulation_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CGAL

#  include <deal.II/grid/tria.h>

#  include <boost/hana.hpp>

#  include <CGAL/Surface_mesh.h>
#  include <CGAL/Triangulation_3.h>
#  include <deal.II/cgal/utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace CGALWrappers
{
  /**
   * Build any CGAL triangulation from a list of deal.II points.
   *
   * The CGAL library offers several types of triangulations, including surface
   * triangulations, cell complexes, Delaunay triangulations in two and three
   * dimensions, and many others. For simple cases, all of these triangulations
   * allow you to incrementally build them from a list of points. This function
   * provides a convenient way to add a vector of deal.II points to any CGAL
   * triangulation that admits the insertion of new points via iterators.
   *
   * The resulting triangulation will be the convex-hull of the input points. If
   * you need to build more complicated triangulations, you can use the
   * functions that build a triangulation from a polyhedral surface domain.
   *
   * More information on the available CGAL triangulation classes is available
   * at https://doc.cgal.org/latest/Triangulation_2/index.html for two
   * dimensional triangulations, and at
   * https://doc.cgal.org/latest/Triangulation_3/index.html for three
   * dimensional triangulations.
   *
   * Notice that CGAL distinguishes between a triangulation and a polygonal or
   * polyhedral mesh. Generally speaking, a triangulation is made of simplices,
   * whereas a polygonal or polyhedral mesh is made of general polygons or
   * general polyhedrons. While CGAL implements the two concepts in a similar
   * fashion using half-edge data structures, some optimizations are performed
   * for triangulations, where only triangles or tetrahedra are used.
   *
   * @param[in] points The input points to build the triangulation from.
   * @param[out] triangulation The output triangulation.
   */
  template <int spacedim, typename CGALTriangulation>
  void
  add_points_to_cgal_triangulation(const std::vector<Point<spacedim>> &points,
                                   CGALTriangulation &triangulation);

  /**
   * Convert any compatible CGAL triangulation type to a deal.II triangulation.
   *
   * The conversion is done by looping over the finite vertices and finite cells
   * of the CGAL triangulation, and building with them a deal.II Triangulation
   * object.
   *
   * CGAL considers triangulations to be a partition of the entire
   * spacedim-dimensional space, where every face (including boundary faces) is
   * shared between neighboring cells. This is obtained by adding to the list of
   * vertices a special vertex (the vertex at infinity, representing a point on
   * the spacedim-dimensional sphere with infinite radius), and adding a
   * neighboring infinite cell to each boundary face.
   *
   * Valid conversions require that the capacity of the deal.II triangulation
   * matches that of the input CGAL triangulation, i.e., if the input CGAL
   * triangulation is a CGAL::Triangulation_2, then spacedim must be greater
   * or equal than two. If the input CGAL triangulation is a
   * CGAL::Triangulation_3, then spacedim must be greater or equal than three.
   *
   * Both CGAL::Triangulation_2 and CGAL::Triangulation_3 can store degenerate
   * two dimensional triangulations (i.e., a one-dimensional triangulation), or
   * a degenerate three dimensional triangulation (i.e., a two-dimensional
   * surface triangulation). The dimension of the triangulation, as returned by
   * `cgal_triangulation.dimension()`, must match the dimension `dim` of the
   * input Triangulation `dealii_triangulation`. If this is not the case, an
   * exception is thrown.
   *
   * @param[in] cgal_triangulation The input CGAL triangulation.
   * @param[out] dealii_triangulation The output deal.II triangulation.
   */
  template <typename CGALTriangulation, int dim, int spacedim>
  void
  cgal_triangulation_to_dealii_triangulation(
    const CGALTriangulation &     cgal_triangulation,
    Triangulation<dim, spacedim> &dealii_triangulation);



#  ifndef DOXYGEN
  // Template implementation
  template <int spacedim, typename CGALTriangulation>
  void
  add_points_to_cgal_triangulation(const std::vector<Point<spacedim>> &points,
                                   CGALTriangulation &triangulation)
  {
    Assert(triangulation.is_valid(),
           ExcMessage(
             "The triangulation you pass to this function should be a valid "
             "CGAL triangulation."));

    std::vector<typename CGALTriangulation::Point> cgal_points(points.size());
    std::transform(points.begin(),
                   points.end(),
                   cgal_points.begin(),
                   [](const auto &p) {
                     return CGALWrappers::dealii_point_to_cgal_point<
                       typename CGALTriangulation::Point>(p);
                   });

    triangulation.insert(cgal_points.begin(), cgal_points.end());
    Assert(triangulation.is_valid(),
           ExcMessage(
             "The Triangulation is no longer valid after inserting the points. "
             "Bailing out."));
  }



  template <typename CGALTriangulation, int dim, int spacedim>
  void
  cgal_triangulation_to_dealii_triangulation(
    const CGALTriangulation &     cgal_triangulation,
    Triangulation<dim, spacedim> &dealii_triangulation)
  {
    AssertThrow(cgal_triangulation.dimension() == dim,
                ExcMessage("The dimension of the input CGAL triangulation (" +
                           std::to_string(cgal_triangulation.dimension()) +
                           ") does not match the dimension of the output "
                           "deal.II triangulation (" +
                           std::to_string(dim) + ")."));

    Assert(dealii_triangulation.n_cells() == 0,
           ExcMessage("The output triangulation object needs to be empty."));

    // deal.II storage data structures
    std::vector<Point<spacedim>> vertices(
      cgal_triangulation.number_of_vertices());
    std::vector<CellData<dim>> cells;
    SubCellData                subcell_data;

    // CGAL storage data structures
    std::map<typename CGALTriangulation::Vertex_handle, unsigned int>
      vertex_map;
    {
      unsigned int i = 0;
      for (const auto &v : cgal_triangulation.finite_vertex_handles())
        {
          vertices[i] =
            CGALWrappers::cgal_point_to_dealii_point<spacedim>(v->point());
          vertex_map[v] = i++;
        }
    }

    const auto has_faces = boost::hana::is_valid(
      [](auto &&obj) -> decltype(obj.finite_face_handles()) {});

    const auto has_cells = boost::hana::is_valid(
      [](auto &&obj) -> decltype(obj.finite_cell_handles()) {});

    // Different loops for Triangulation_2 and Triangulation_3 types.
    if constexpr (decltype(has_faces(cgal_triangulation)){})
      {
        // This is a non-degenerate Triangulation_2
        if (cgal_triangulation.dimension() == 2)
          for (const auto &f : cgal_triangulation.finite_face_handles())
            {
              CellData<dim> cell(ReferenceCells::Triangle.n_vertices());
              for (unsigned int i = 0;
                   i < ReferenceCells::Triangle.n_vertices();
                   ++i)
                cell.vertices[i] = vertex_map[f->vertex(i)];
              cells.push_back(cell);
            }
        else if (cgal_triangulation.dimension() == 1)
          // This is a degenerate Triangulation_2, made of edges
          for (const auto &e : cgal_triangulation.finite_edges())
            {
              // An edge is idenfied by a face and a vertex index in the face
              const auto &  f = e.first;
              const auto &  i = e.second;
              CellData<dim> cell(ReferenceCells::Line.n_vertices());
              unsigned int  id = 0;
              // Since an edge is identified by a face (a triangle) and the
              // index of the opposite vertex to the edge, we can use this logic
              // to infer the indices of the vertices of the edge: loop over all
              // vertices, and keep only those that are not the opposite vertex
              // of the edge.
              for (unsigned int j = 0;
                   j < ReferenceCells::Triangle.n_vertices();
                   ++j)
                if (j != i)
                  cell.vertices[id++] = vertex_map[f->vertex(j)];
              cells.push_back(cell);
            }
        else
          {
            Assert(false, ExcInternalError());
          }
      }
    else if constexpr (decltype(has_cells(cgal_triangulation)){})
      {
        // This is a non-degenerate Triangulation_3
        if (cgal_triangulation.dimension() == 3)
          for (const auto &c : cgal_triangulation.finite_cell_handles())
            {
              CellData<dim> cell(ReferenceCells::Tetrahedron.n_vertices());
              for (unsigned int i = 0;
                   i < ReferenceCells::Tetrahedron.n_vertices();
                   ++i)
                cell.vertices[i] = vertex_map[c->vertex(i)];
              cells.push_back(cell);
            }
        else if (cgal_triangulation.dimension() == 2)
          // This is a degenerate Triangulation_3, made of triangles
          for (const auto &facet : cgal_triangulation.finite_facets())
            {
              // A facet is idenfied by a cell and the opposite vertex index in
              // the face
              const auto &  c = facet.first;
              const auto &  i = facet.second;
              CellData<dim> cell(ReferenceCells::Triangle.n_vertices());
              unsigned int  id = 0;
              // Since a face is identified by a cell (a tetrahedron) and the
              // index of the opposite vertex to the face, we can use this logic
              // to infer the indices of the vertices of the face: loop over all
              // vertices, and keep only those that are not the opposite vertex
              // of the face.
              for (unsigned int j = 0;
                   j < ReferenceCells::Tetrahedron.n_vertices();
                   ++j)
                if (j != i)
                  cell.vertices[id++] = vertex_map[c->vertex(j)];
              cells.push_back(cell);
            }
        else if (cgal_triangulation.dimension() == 1)
          // This is a degenerate Triangulation_3, made of edges
          for (const auto &edge : cgal_triangulation.finite_edges())
            {
              // An edge is idenfied by a cell and its two vertices
              const auto &[c, i, j] = edge;
              CellData<dim> cell(ReferenceCells::Line.n_vertices());
              cell.vertices[0] = vertex_map[c->vertex(i)];
              cell.vertices[1] = vertex_map[c->vertex(j)];
              cells.push_back(cell);
            }
        else
          {
            Assert(false, ExcInternalError());
          }
      }
    dealii_triangulation.create_triangulation(vertices, cells, subcell_data);
  }
#  endif
} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
