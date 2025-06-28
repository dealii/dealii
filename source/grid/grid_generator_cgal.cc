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

#include <deal.II/base/ndarray.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_description.h>

#ifdef DEAL_II_WITH_CGAL
// Functions needed by the CGAL mesh generation utilities are inside
#  include <deal.II/cgal/triangulation.h>

#  include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#  include <CGAL/Labeled_mesh_domain_3.h>
#  include <CGAL/Mesh_triangulation_3.h>
#endif


DEAL_II_NAMESPACE_OPEN

// work around the problem that doxygen for some reason lists all template
// specializations in this file
#ifndef DOXYGEN

namespace GridGenerator
{
  template <int dim>
  void
  implicit_function(Triangulation<dim, 3> &tria,
                    const Function<3>     &dealii_implicit_function,
                    const CGALWrappers::AdditionalData<dim> &data,
                    const Point<3>                          &interior_point,
                    const double                            &outer_ball_radius)
  {
#  ifdef DEAL_II_WITH_CGAL
    Assert(dealii_implicit_function.n_components == 1,
           ExcMessage(
             "The implicit function must have exactly one component."));
    Assert(dealii_implicit_function.value(interior_point) < 0,
           ExcMessage(
             "The implicit function must be negative at the interior point."));
    Assert(outer_ball_radius > 0,
           ExcMessage("The outer ball radius must be positive."));
    Assert(tria.n_active_cells() == 0,
           ExcMessage("The triangulation must be empty."));

    if constexpr (dim == 3)
      {
        using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
        using NumberType = K::FT;
        using Point_3    = K::Point_3;
        using Sphere_3   = K::Sphere_3;

        using Mesh_domain = CGAL::Labeled_mesh_domain_3<K>;
        using Tr =
          CGAL::Mesh_triangulation_3<Mesh_domain,
                                     CGAL::Default,
                                     CGALWrappers::ConcurrencyTag>::type;
        using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr, int, int>;
        using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;


        auto cgal_implicit_function = [&](const Point_3 &p) {
          return NumberType(
            dealii_implicit_function.value(Point<3>(p.x(), p.y(), p.z())));
        };

        Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(
          cgal_implicit_function,
          K::Sphere_3(
            Point_3(interior_point[0], interior_point[1], interior_point[2]),
            outer_ball_radius * outer_ball_radius));

        Mesh_criteria criteria(CGAL::parameters::facet_size  = data.facet_size,
                               CGAL::parameters::facet_angle = data.facet_angle,
                               CGAL::parameters::facet_distance =
                                 data.facet_distance,
                               CGAL::parameters::cell_radius_edge_ratio =
                                 data.cell_radius_edge_ratio,
                               CGAL::parameters::cell_size = data.cell_size);

        auto cgal_triangulation = CGAL::make_mesh_3<C3t3>(domain, criteria);
        CGALWrappers::cgal_triangulation_to_dealii_triangulation(
          cgal_triangulation, tria);
      }
    else if constexpr (dim == 2)
      {
        // default triangulation for Surface_mesher
        using Tr           = CGAL::Surface_mesh_default_triangulation_3;
        using C2t3         = CGAL::Complex_2_in_triangulation_3<Tr>;
        using GT           = Tr::Geom_traits;
        using Sphere_3     = GT::Sphere_3;
        using Point_3      = GT::Point_3;
        using FT           = GT::FT;
        using Function     = FT (*)(Point_3);
        using Surface_3    = CGAL::Implicit_surface_3<GT, Function>;
        using Surface_mesh = CGAL::Surface_mesh<Point_3>;


        auto cgal_implicit_function = [&](const Point_3 &p) {
          return FT(
            dealii_implicit_function.value(Point<3>(p.x(), p.y(), p.z())));
        };

        Surface_3 surface(cgal_implicit_function,
                          Sphere_3(Point_3(interior_point[0],
                                           interior_point[1],
                                           interior_point[2]),
                                   outer_ball_radius * outer_ball_radius));

        Tr           tr;
        C2t3         c2t3(tr);
        Surface_mesh mesh;

        CGAL::Surface_mesh_default_criteria_3<Tr> criteria(data.angular_bound,
                                                           data.radius_bound,
                                                           data.distance_bound);
        CGAL::make_surface_mesh(c2t3,
                                surface,
                                criteria,
                                CGAL::Non_manifold_tag());
        CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, mesh);
        CGALWrappers::cgal_surface_mesh_to_dealii_triangulation(mesh, tria);
      }
    else
      {
        Assert(false, ExcImpossibleInDim(dim));
      }

#  else

    (void)tria;
    (void)dealii_implicit_function;
    (void)data;
    (void)interior_point;
    (void)outer_ball_radius;
    AssertThrow(false, ExcMessage("This function needs CGAL to be installed."));

#  endif
  }



  void
  surface_mesh_to_volumetric_mesh(const Triangulation<2, 3> &surface_tria,
                                  Triangulation<3>          &vol_tria,
                                  const CGALWrappers::AdditionalData<3> &data)
  {
#  ifdef DEAL_II_WITH_CGAL
    Assert(
      surface_tria.n_cells() > 0,
      ExcMessage(
        "The input triangulation cannot be empty when calling this function."));
    Assert(
      vol_tria.n_cells() == 0,
      ExcMessage(
        "The output triangulation must be empty when calling this function."));
    using K       = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point_3 = K::Point_3;

    using Mesh_domain =
      CGAL::Polyhedral_mesh_domain_with_features_3<K,
                                                   CGAL::Surface_mesh<Point_3>>;
    using Tr            = CGAL::Mesh_triangulation_3<Mesh_domain,
                                          CGAL::Default,
                                          CGALWrappers::ConcurrencyTag>::type;
    using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
    using C3t3 =
      CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                              Mesh_domain::Corner_index,
                                              Mesh_domain::Curve_index>;

    CGAL::Surface_mesh<Point_3> mesh;
    // This function "fills" the missing arrow of the following diagram.
    //  Tria<2,3>                           Tria<3>
    //      |                                 ^
    //      |                                 |
    //      |                                 |
    //      |                                 |
    //      V                                 |
    // CGAL::Surface_mesh -----------> CGAL::C3t3
    CGALWrappers::dealii_tria_to_cgal_surface_mesh(surface_tria, mesh);
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);
    Mesh_domain domain(mesh);
    domain.detect_features();
    Mesh_criteria criteria(CGAL::parameters::facet_size  = data.facet_size,
                           CGAL::parameters::facet_angle = data.facet_angle,
                           CGAL::parameters::facet_distance =
                             data.facet_distance,
                           CGAL::parameters::cell_radius_edge_ratio =
                             data.cell_radius_edge_ratio,
                           CGAL::parameters::cell_size = data.cell_size);
    const auto cgal_triangulation = CGAL::make_mesh_3<C3t3>(domain, criteria);
    CGALWrappers::cgal_triangulation_to_dealii_triangulation(cgal_triangulation,
                                                             vol_tria);

#  else

    (void)surface_tria;
    (void)vol_tria;
    (void)data;
    AssertThrow(false, ExcMessage("This function needs CGAL to be installed."));

#  endif
  }
} // namespace GridGenerator

// explicit instantiations
#  include "grid/grid_generator_cgal.inst"

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE
