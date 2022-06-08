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

#ifndef dealii_cgal_additional_data_h
#define dealii_cgal_additional_data_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CGAL

#  include <CGAL/boost/graph/Named_function_parameters.h>

DEAL_II_NAMESPACE_OPEN
namespace CGALWrappers
{
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
   * `CGAL::FACET_VERTICES_ON_SURFACE`. See CGAL::Mesh_facet_topology manual
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
    double facet_angle;
    double facet_size;
    double facet_distance;
    double cell_radius_edge_ratio;
    double cell_size;

    AdditionalData(
      double facet_a            = CGAL::parameters::is_default_parameter(true),
      double facet_s            = CGAL::parameters::is_default_parameter(true),
      double facet_d            = CGAL::parameters::is_default_parameter(true),
      double cell_radius_edge_r = CGAL::parameters::is_default_parameter(true),
      double cell_s             = CGAL::parameters::is_default_parameter(true))
    {
      AssertThrow(
        dim == 3,
        ExcMessage(
          "These struct can be instantiated with 3D Triangulations only."));
      facet_angle            = facet_a;
      facet_size             = facet_s;
      facet_distance         = facet_d;
      cell_radius_edge_ratio = cell_radius_edge_r;
      cell_size              = cell_s;
    }
  };



  /**
   * Specialization of the above struct when the object to be constructed is a
   * 2D triangulation embedded in the 3D space, i.e. a Triangulation<2,3>.
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
    double angular_bound, radius_bound, distance_bound;
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
  template <int dim>
  struct AdditionalData
  {};
} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE
#endif

#endif
