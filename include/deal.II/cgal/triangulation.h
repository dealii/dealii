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

#  include <CGAL/Polyhedron_3.h>
#  include <CGAL/Surface_mesh.h>
#  include <CGAL/Triangulation_2.h>
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

  /**
   * Construct a deal.II surface triangulation starting from a
   * CGAL::Surface_mesh or a CGAL::Polyhedron_3.
   *
   * These types can both represent a polyhedral surface made of general
   * polygons. deal.II only supports triangle or quadrilateral triangulations,
   * and this function will throw an exception if the input surface mesh
   * contains polygonal faces with more than 4 vertices per face.
   *
   * CGAL orients the faces of a polyhedral surface in a counter clock-wise
   * w.r.t. to the material side of the surface, i.e., the normal is directed
   * towards the outside of the surface if the surface is closed and represents
   * a bounded domain. The opposite orientation is also perfectly valid, and
   * would represent an unbounded domain complementary to the region
   * represented by the surface. This function does not change the orientation
   * of the surface mesh and will produce a triangulation with the same
   * orientation of the input mesh.
   *
   * @param[in] cgal_mesh The input CGAL surface mesh.
   * @param[out] triangulation The output deal.II triangulation.
   */
  template <typename CGAL_MeshType>
  void
  cgal_surface_mesh_to_dealii_triangulation(const CGAL_MeshType &cgal_mesh,
                                            Triangulation<2, 3> &triangulation);



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
              // An edge is idendified by a face and a vertex index in the face
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



  template <typename CGAL_MeshType>
  void
  cgal_surface_mesh_to_dealii_triangulation(const CGAL_MeshType &cgal_mesh,
                                            Triangulation<2, 3> &triangulation)
  {
    Assert(triangulation.n_cells() == 0,
           ExcMessage(
             "Triangulation must be empty upon calling this function."));

    const auto is_surface_mesh =
      boost::hana::is_valid([](auto &&obj) -> decltype(obj.faces()) {});

    const auto is_polyhedral =
      boost::hana::is_valid([](auto &&obj) -> decltype(obj.facets_begin()) {});

    // Collect Vertices and cells
    std::vector<dealii::Point<3>> vertices;
    std::vector<CellData<2>>      cells;
    SubCellData                   subcell_data;

    // Different loops for Polyhedron or Surface_mesh types
    if constexpr (decltype(is_surface_mesh(cgal_mesh)){})
      {
        AssertThrow(cgal_mesh.num_vertices() > 0,
                    ExcMessage("CGAL surface mesh is empty."));
        vertices.reserve(cgal_mesh.num_vertices());
        std::map<typename CGAL_MeshType::Vertex_index, unsigned int> vertex_map;
        {
          unsigned int i = 0;
          for (const auto &v : cgal_mesh.vertices())
            {
              vertices.emplace_back(CGALWrappers::cgal_point_to_dealii_point<3>(
                cgal_mesh.point(v)));
              vertex_map[v] = i++;
            }
        }

        // Collect CellData
        for (const auto &face : cgal_mesh.faces())
          {
            const auto face_vertices =
              CGAL::vertices_around_face(cgal_mesh.halfedge(face), cgal_mesh);

            AssertThrow(face_vertices.size() == 3 || face_vertices.size() == 4,
                        ExcMessage("Only triangle or quadrilateral surface "
                                   "meshes are supported in deal.II"));

            CellData<2> c(face_vertices.size());
            auto        it_vertex = c.vertices.begin();
            for (const auto &v : face_vertices)
              *(it_vertex++) = vertex_map[v];

            // Make sure the numberfing is consistent with the one in deal.II
            if (face_vertices.size() == 4)
              std::swap(c.vertices[3], c.vertices[2]);

            cells.emplace_back(c);
          }
      }
    else if constexpr (decltype(is_polyhedral(cgal_mesh)){})
      {
        AssertThrow(cgal_mesh.size_of_vertices() > 0,
                    ExcMessage("CGAL surface mesh is empty."));
        vertices.reserve(cgal_mesh.size_of_vertices());
        std::map<decltype(cgal_mesh.vertices_begin()), unsigned int> vertex_map;
        {
          unsigned int i = 0;
          for (auto it = cgal_mesh.vertices_begin();
               it != cgal_mesh.vertices_end();
               ++it)
            {
              vertices.emplace_back(
                CGALWrappers::cgal_point_to_dealii_point<3>(it->point()));
              vertex_map[it] = i++;
            }
        }

        // Loop over faces of Polyhedron, fill CellData
        for (auto face = cgal_mesh.facets_begin();
             face != cgal_mesh.facets_end();
             ++face)
          {
            auto               j                 = face->facet_begin();
            const unsigned int vertices_per_face = CGAL::circulator_size(j);
            AssertThrow(vertices_per_face == 3 || vertices_per_face == 4,
                        ExcMessage("Only triangle or quadrilateral surface "
                                   "meshes are supported in deal.II. You "
                                   "tried to read a mesh where a face has " +
                                   std::to_string(vertices_per_face) +
                                   " vertices per face."));

            CellData<2> c(vertices_per_face);
            auto        it = c.vertices.begin();
            for (unsigned int i = 0; i < vertices_per_face; ++i)
              {
                *(it++) = vertex_map[j->vertex()];
                ++j;
              }

            if (vertices_per_face == 4)
              std::swap(c.vertices[3], c.vertices[2]);

            cells.emplace_back(c);
          }
      }
    else
      {
        AssertThrow(false,
                    ExcInternalError(
                      "Unsupported CGAL surface triangulation type."));
      }
    triangulation.create_triangulation(vertices, cells, subcell_data);
  }
#  endif
} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
