// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cgal_additional_data_h
#define dealii_cgal_additional_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#ifdef DEAL_II_WITH_CGAL

#  include <CGAL/version.h>
#  if CGAL_VERSION_MAJOR >= 6
#    include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#  endif
#  include <CGAL/Mesh_facet_topology.h>

#  include <limits>

DEAL_II_NAMESPACE_OPEN
namespace CGALWrappers
{
  enum class FacetTopology
  {
    /**
     * Each vertex of the facet have to be on the surface, on a curve, or on a
     * corner.
     */
    facet_vertices_on_surface = CGAL::FACET_VERTICES_ON_SURFACE,

    /**
     * The three vertices of a facet belonging to a surface patch s have to be
     * on the same surface patch s, on a curve or on a corner.
     */
    facet_vertices_on_same_surface_patch =
      CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,

    /**
     * The three vertices of a facet belonging to a surface patch s have to be
     * on the same surface patch s, or on a curve incident to the surface patch
     * s or on a corner incident to the surface patch s.
     */
    facet_vertices_on_same_surface_patch_with_adjacency_check =
      CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK,
  };

  /**
   * Struct that must be used to pass additional
   * arguments to the CGAL::Mesh_criteria_3 class (see
   * https://doc.cgal.org/latest/Mesh_3/index.html for more information.)
   *
   * The arguments allow for fine control on the size, quality, and distribution
   * of the cells of the final triangulation. CGAL uses Boost named parameters
   * for these arguments in dimension three, i.e., they must be specified with
   * the syntax `CGAL::parameters::parameter_name=parameter_value`, irrespective
   * of their order. Accepted parameters are:
   *
   * - `CGAL::parameters::edge_size`: a constant
   * providing a uniform upper bound for the lengths of
   * curve edges. This parameter has to be set to a positive value when
   * 1-dimensional features protection is used.
   * - `CGAL::parameters::facet_angle`: a lower bound for the angles (in
   * degrees) of the surface mesh facets.
   * - `CGAL::parameters::facet_size`: a constant
   * describing a uniform upper bound for the radii
   * of the surface Delaunay balls.
   * - `CGAL::parameters::facet_distance`: a constant
   * describing a uniform upper bound for the distance
   * between the facet circumcenter and the center of its surface Delaunay ball.
   * - `CGAL::parameters::facet_topology`: the set of topological constraints
   * which have to be verified by each surface facet. The default value is
   * `CGAL::FACET_VERTICES_ON_SURFACE`. See the enum @p FacetToplogy CGAL::Mesh_facet_topology manual
   * page to get all possible values.
   * - `CGAL::parameters::cell_radius_edge_ratio`: an upper bound for the
   * radius-edge ratio of the mesh tetrahedra.
   * - `CGAL::parameters::cell_size`: a constant
   * describing a uniform upper bound for the
   * circumradii of the mesh tetrahedra.
   *
   * @note This struct must be instantiated with `dim=3`.
   */
  template <int dim>
  struct AdditionalData
  {
    /**
     * Uniform upper bound for the lengths of curve edges.
     * This parameter has to be set to a positive value when1-dimensional
     * features protection is used.
     */
    double edge_size;

    /**
     * Lower bound for the angles (in degrees) of the surface mesh facets.
     */
    double facet_angle;

    /**
     * Uniform upper bound for the radii of the surface Delaunay balls.
     */
    double facet_size;
    /**
     * Uniform upper bound for the distance between the facet circumcenter and
     * the center of its surface Delaunay ball.
     */
    double facet_distance;

    /**
     * Set of topological constraints which have to be verified by each surface
     * facet.
     */
    FacetTopology facet_topology;

    /**
     * upper bound for the radius-edge ratio of the mesh tetrahedra.
     */
    double cell_radius_edge_ratio;
    /**
     * Uniform upper bound for the circumradii of the mesh tetrahedra.
     */
    double cell_size;

    /**
     * Constructor.
     */
    AdditionalData(
      double        edge_s  = std::numeric_limits<double>::max(),
      double        facet_a = 0.,
      double        facet_s = 0.,
      double        facet_d = 0.,
      FacetTopology facet_t =
        dealii::CGALWrappers::FacetTopology::facet_vertices_on_surface,
      double cell_radius_edge_r = 0.,
      double cell_s             = 0.)
    {
      AssertThrow(
        dim == 3,
        ExcMessage(
          "These struct can be instantiated with 3d Triangulations only."));
      edge_size              = edge_s;
      facet_angle            = facet_a;
      facet_size             = facet_s;
      facet_distance         = facet_d;
      facet_topology         = facet_t;
      cell_radius_edge_ratio = cell_radius_edge_r;
      cell_size              = cell_s;
    }
  };

  /**
   * Specialization of the above struct when the object to be constructed is a
   * 2d triangulation embedded in the 3d space, i.e. a Triangulation<2,3>.
   * Only three parameters are accepted:
   * - `angular_bound` is a lower bound in degrees for the angles of mesh
   * facets.
   * - `radius_bound` is an upper bound on the radii of surface Delaunay balls.
   * A surface Delaunay ball is a ball circumscribing a mesh facet and centered
   * on the surface.
   * - `distance_bound` is an upper bound for the distance between the
   * circumcenter of a mesh facet and the center of a surface Delaunay ball of
   * this facet.
   */
  template <>
  struct AdditionalData<2>
  {
    /**
     * Lower bound in degrees for the angles of mesh facets.
     */
    double angular_bound;

    /**
     * Upper bound on the radii of surface Delaunay balls. A surface Delaunay
     * ball is a ball circumscribing a mesh facet and centered on the surface.
     */
    double radius_bound;

    /**
     * Upper bound for the distance between the circumcenter of a mesh facet
     * and the center of a surface Delaunay ball of this facet.
     */
    double distance_bound;

    /**
     * Constructor.
     */
    AdditionalData(double angular_b  = 0.,
                   double radius_b   = 0.,
                   double distance_b = 0.)
    {
      angular_bound  = angular_b;
      radius_bound   = radius_b;
      distance_bound = distance_b;
    }
  };

} // namespace CGALWrappers
DEAL_II_NAMESPACE_CLOSE

#else

DEAL_II_NAMESPACE_OPEN
namespace CGALWrappers
{
  /**
   * Empty structure for the case where CGAL is not available.
   */
  template <int dim>
  struct AdditionalData
  {};
} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE
#endif

#endif
