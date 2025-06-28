// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
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

#include <deal.II/cgal/surface_mesh.h>

#ifdef DEAL_II_WITH_CGAL

#  include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#  include <CGAL/Simple_cartesian.h>
#  include <CGAL/Surface_mesh/Surface_mesh.h>


DEAL_II_NAMESPACE_OPEN

namespace
{
  template <typename dealiiFace, typename CGAL_Mesh>
  void
  add_facet(
    const dealiiFace                                               &face,
    const std::map<unsigned int, typename CGAL_Mesh::Vertex_index> &deal2cgal,
    CGAL_Mesh                                                      &mesh,
    const bool clockwise_ordering = true)
  {
    const auto reference_cell_type = face->reference_cell();
    std::vector<typename CGAL_Mesh::Vertex_index> indices;

    if (reference_cell_type == ReferenceCells::Line)
      {
        mesh.add_edge(deal2cgal.at(face->vertex_index(0)),
                      deal2cgal.at(face->vertex_index(1)));

        if (clockwise_ordering)
          std::reverse(indices.begin(), indices.end());
      }
    else if (reference_cell_type == ReferenceCells::Triangle)
      {
        indices = {deal2cgal.at(face->vertex_index(0)),
                   deal2cgal.at(face->vertex_index(1)),
                   deal2cgal.at(face->vertex_index(2))};

        if (clockwise_ordering)
          std::reverse(indices.begin(), indices.end());
      }

    else if (reference_cell_type == ReferenceCells::Quadrilateral)
      {
        if (clockwise_ordering)
          {
            indices = {deal2cgal.at(face->vertex_index(0)),
                       deal2cgal.at(face->vertex_index(2)),
                       deal2cgal.at(face->vertex_index(3)),
                       deal2cgal.at(face->vertex_index(1))};
          }
        else
          {
            indices = {deal2cgal.at(face->vertex_index(0)),
                       deal2cgal.at(face->vertex_index(1)),
                       deal2cgal.at(face->vertex_index(3)),
                       deal2cgal.at(face->vertex_index(2))};
          }
      }
    else
      DEAL_II_ASSERT_UNREACHABLE();

    [[maybe_unused]] const auto new_face = mesh.add_face(indices);
    Assert(new_face != mesh.null_face(),
           ExcInternalError("While trying to build a CGAL facet, "
                            "CGAL encountered a orientation problem that it "
                            "was not able to solve."));
  }



  template <typename dealiiFace, typename CGAL_Mesh>
  void
  map_vertices(
    const dealiiFace                                         &cell,
    std::map<unsigned int, typename CGAL_Mesh::Vertex_index> &deal2cgal,
    CGAL_Mesh                                                &mesh)
  {
    for (const auto i : cell->vertex_indices())
      {
        if (deal2cgal.find(cell->vertex_index(i)) == deal2cgal.end())
          {
            deal2cgal[cell->vertex_index(i)] =
              mesh.add_vertex(CGALWrappers::dealii_point_to_cgal_point<
                              typename CGAL_Mesh::Point>(cell->vertex(i)));
          }
      }
  }
} // namespace



#  ifndef DOXYGEN
// Template implementations
namespace CGALWrappers
{
  template <typename CGALPointType, int dim, int spacedim>
  void
  dealii_cell_to_cgal_surface_mesh(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Mapping<dim, spacedim>                               &mapping,
    CGAL::Surface_mesh<CGALPointType>                          &mesh)
  {
    Assert(dim > 1, ExcImpossibleInDim(dim));
    using Mesh           = CGAL::Surface_mesh<CGALPointType>;
    const auto &vertices = mapping.get_vertices(cell);
    std::map<unsigned int, typename Mesh::Vertex_index> deal2cgal;

    // Add all vertices to the mesh
    // Store CGAL ordering
    for (const auto &i : cell->vertex_indices())
      deal2cgal[cell->vertex_index(i)] = mesh.add_vertex(
        CGALWrappers::dealii_point_to_cgal_point<CGALPointType>(vertices[i]));

    // Add faces
    if (dim < 3)
      // simplices and quads are allowable faces for CGAL
      add_facet(cell, deal2cgal, mesh);
    else
      // in 3d, we build a surface mesh containing all the faces of the 3d cell.
      // Simplices, Tetrahedrons, and Pyramids have their faces numbered in the
      // same way as CGAL does (all faces of a bounding polyhedron are numbered
      // counter-clockwise, so that their normal points outwards). Hexahedrons
      // in deal.II have their faces numbered lexicographically, and one cannot
      // deduce the direction of the normals by just looking at the vertices.
      //
      // In order for CGAL to be able to produce the right orientation, we need
      // to reverse the order of the vertices for faces with even index.
      // However, in order to allow for all kinds of meshes in 3d, including
      // Moebius-loops, a deal.II face might even be rotated looking from one
      // cell, whereas it is according to the standard when looking at it from
      // the neighboring cell sharing that particular face. Therefore, when
      // building a cgal face we must also take into account the fact that a
      // face may have a non-standard orientation.
      for (const auto &f : cell->face_indices())
        {
          // Check for standard orientation of faces
          bool face_is_clockwise_oriented =
            cell->reference_cell() != ReferenceCells::Hexahedron ||
            (f % 2 == 0);
          // Make sure that we revert the orientation if required
          if (cell->face_orientation(f) == false)
            face_is_clockwise_oriented = !face_is_clockwise_oriented;
          add_facet(cell->face(f), deal2cgal, mesh, face_is_clockwise_oriented);
        }
  }



  template <typename CGALPointType, int dim, int spacedim>
  void
  dealii_tria_to_cgal_surface_mesh(
    const dealii::Triangulation<dim, spacedim> &tria,
    CGAL::Surface_mesh<CGALPointType>          &mesh)
  {
    Assert(tria.n_cells() > 0,
           ExcMessage(
             "Triangulation cannot be empty upon calling this function."));
    Assert(mesh.is_empty(),
           ExcMessage(
             "The surface mesh must be empty upon calling this function."));

    Assert(dim > 1, ExcImpossibleInDim(dim));
    using Mesh         = CGAL::Surface_mesh<CGALPointType>;
    using Vertex_index = typename Mesh::Vertex_index;

    std::map<unsigned int, Vertex_index> deal2cgal;
    if constexpr (dim == 2)
      {
        for (const auto &cell : tria.active_cell_iterators())
          {
            map_vertices(cell, deal2cgal, mesh);
            add_facet(cell, deal2cgal, mesh);
          }
      }
    else if constexpr (dim == 3 && spacedim == 3)
      {
        for (const auto &cell : tria.active_cell_iterators())
          {
            for (const auto &f : cell->face_indices())
              {
                if (cell->face(f)->at_boundary())
                  {
                    map_vertices(cell->face(f), deal2cgal, mesh);
                    add_facet(cell->face(f),
                              deal2cgal,
                              mesh,
                              (f % 2 == 0 || cell->n_vertices() != 8));
                  }
              }
          }
      }
    else
      {
        Assert(false, ExcImpossibleInDimSpacedim(dim, spacedim));
      }
  }

// Explicit instantiations.
//
// We don't build the instantiations.inst file if deal.II isn't
// configured with CGAL, but doxygen doesn't know that and tries to
// find that file anyway for parsing -- which then of course it fails
// on. So exclude the following from doxygen consideration.
#    ifndef DOXYGEN
#      include "cgal/surface_mesh.inst"
#    endif

} // namespace CGALWrappers
#  endif


DEAL_II_NAMESPACE_CLOSE

#endif
