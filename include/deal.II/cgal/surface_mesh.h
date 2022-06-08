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

#ifndef dealii_cgal_surface_mesh_h
#define dealii_cgal_surface_mesh_h

#include <deal.II/base/config.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>


#ifdef DEAL_II_WITH_CGAL
#  include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#  include <CGAL/Surface_mesh.h>
#  include <deal.II/cgal/point_conversion.h>


DEAL_II_NAMESPACE_OPEN

namespace CGALWrappers
{
  /**
   * Build a CGAL::Surface_mesh from a deal.II cell.
   *
   * The class Surface_mesh implements a halfedge data structure and can be used
   * to represent polyhedral surfaces. It is an edge-centered data structure
   * capable of maintaining incidence information of vertices, edges, and faces.
   * Each edge is represented by two halfedges with opposite orientation. The
   * orientation of a face is chosen so that the halfedges around a face are
   * oriented counterclockwise.
   *
   * More information on this class is available at
   * https://doc.cgal.org/latest/Surface_mesh/index.html
   *
   * The function will throw an exception in dimension one. In dimension two, it
   * generates a surface mesh of the quadrilateral cell or of the triangle cell,
   * while in dimension three it will generate the surface mesh of the cell,
   * i.e., a polyhedral mesh containing the faces of the input cell.
   *
   * The generated mesh is useful when performing geometric operations using
   * CGAL::Polygon_mesh_processing, i.e., to compute boolean operations on
   * cells, splitting, cutting, slicing, etc.
   *
   * For examples on how to use the resulting  CGAL::Surface_mesh see
   * https://doc.cgal.org/latest/Polygon_mesh_processing/
   *
   * @param[in] cell The input deal.II cell iterator
   * @param[in] mapping The mapping used to map the vertices of the cell
   * @param[out] mesh The output CGAL::Surface_mesh
   */
  template <typename CGALPointType, int dim, int spacedim>
  void
  dealii_cell_to_cgal_surface_mesh(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
    const dealii::Mapping<dim, spacedim> &                              mapping,
    CGAL::Surface_mesh<CGALPointType> &                                 mesh);

  /**
   * Convert a deal.II triangulation to a CGAL::Surface_mesh. The output depends
   * on the intrinsic dimension of the input deal.II triangulation.
   *
   * In 2D, i.e. with a
   * Triangulation<2> or a Triangulation<2,3>, the output is the
   * CGAL::Surface_mesh describing the whole triangulation.
   *
   * In 3D, the boundary the of the deal.II Triangulation is converted to
   * a CGAL::Surface_mesh by looping over all the boundary faces.
   *
   * @param[in] triangulation The input deal.II triangulation.
   * @param[out] mesh The output CGAL::Surface_mesh.
   */
  template <typename CGALPointType, int dim, int spacedim>
  void
  dealii_tria_to_cgal_surface_mesh(
    const dealii::Triangulation<dim, spacedim> &triangulation,
    CGAL::Surface_mesh<CGALPointType> &         mesh);
} // namespace CGALWrappers



DEAL_II_NAMESPACE_CLOSE

#endif
#endif
