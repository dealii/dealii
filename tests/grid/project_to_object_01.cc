// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <array>
#include <cmath>
#include <numeric>

#include "../tests.h"

int
main()
{
  initlog();

  // test the internal cross derivative function: the stencil should be order 2.
  {
    using namespace dealii::GridTools::internal::ProjectToObject;

    constexpr int structdim = 2;
    auto          objective =
      [](const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &weights) {
        return std::sin(2.0 * weights[0]) * std::cos(3.0 * weights[1]) +
               std::cos(weights[2]);
      };
    Tensor<1, GeometryInfo<structdim>::vertices_per_cell> c0;
    for (const unsigned int row_n : GeometryInfo<structdim>::vertex_indices())
      c0[row_n] = 1.0;

    const std::array<double, 3> exact{
      {2.0 * std::cos(3.0) * std::cos(2.0) +
         3.0 * std::sin(3.0) * std::sin(2.0),
       -2.0 * std::cos(3.0) * std::cos(2.0) -
         3.0 * std::sin(3.0) * std::sin(2.0),
       -3.0 * std::sin(3.0) * std::sin(2.0) + std::sin(1.0)}};
    const std::array<CrossDerivative, 3> cross_derivatives{
      {{0, 1}, {1, 0}, {1, 2}}};

    for (unsigned int cross_derivative_n = 0;
         cross_derivative_n < cross_derivatives.size();
         ++cross_derivative_n)
      {
        deallog << "testing cross derivative " << cross_derivative_n
                << std::endl;
        std::vector<double> steps;
        std::vector<double> errors;
        for (unsigned int step_n = 5; step_n < 10; ++step_n)
          {
            const double current_step = std::pow(0.5, double(step_n));
            steps.push_back(std::log(current_step));
            errors.push_back(std::log(std::abs(
              cross_stencil<structdim>(cross_derivatives[cross_derivative_n],
                                       c0,
                                       current_step,
                                       objective) -
              exact[cross_derivative_n])));
          }

        const double mean_step =
          std::accumulate(steps.begin(), steps.end(), 0.0) / steps.size();
        const double mean_error =
          std::accumulate(errors.begin(), errors.end(), 0.0) / errors.size();

        double numerator   = 0.0;
        double denominator = 0.0;
        for (std::size_t i = 0; i < steps.size(); ++i)
          {
            numerator += (steps[i] - mean_step) * (errors[i] - mean_error);
            denominator += Utilities::fixed_power<2>(steps[i] - mean_step);
          }
        deallog << "slope is nearly 3: "
                << (std::abs(numerator / denominator - 3.0) < 0.05)
                << std::endl;
      }
  }

  // Test project_to_object in 2D with an annulus
  {
    PolarManifold<2>             polar_manifold;
    Triangulation<2>             triangulation;
    constexpr types::manifold_id polar_id = 42;
    GridGenerator::hyper_shell(triangulation, Point<2>(), 1.0, 2.0);
    triangulation.set_manifold(polar_id, polar_manifold);
    triangulation.set_all_manifold_ids(42);

    triangulation.refine_global(2);

    auto cell = triangulation.begin_active();
    while (!cell->at_boundary())
      ++cell;
    unsigned int face_n = 0;
    while (!cell->face(face_n)->at_boundary())
      ++face_n;
    const auto face = cell->face(face_n);

    // easy test: project the arithmetic mean of the vertices onto the
    // geodesic
    {
      const Point<2> trial_point = 0.9 * face->center();
      const Point<2> projected_point =
        GridTools::project_to_object(face, trial_point);

      const std::array<Point<2>, 2> vertices{
        {face->vertex(0), face->vertex(1)}};
      const std::array<double, 2> weights{{0.5, 0.5}};
      const auto                  vertices_view =
        make_array_view(vertices.begin(), vertices.end());
      const auto weights_view = make_array_view(weights.begin(), weights.end());

      deallog << std::endl
              << "Project the arithmetic mean of two vertices on a circular"
              << std::endl
              << "arc. The result should be the geodesic midpoint:" << std::endl
              << std::endl
              << "vertex 0 distance from origin: "
              << face->vertex(0).distance(Point<2>()) << std::endl
              << "vertex 1 distance from origin: "
              << face->vertex(1).distance(Point<2>()) << std::endl
              << "trial point distance from origin: "
              << trial_point.distance(Point<2>()) << std::endl
              << "projected point distance from origin: "
              << projected_point.distance(Point<2>()) << std::endl
              << "projected point:        " << projected_point << std::endl
              << "true geodesic midpoint: "
              << face->get_manifold().get_new_point(vertices_view, weights_view)
              << std::endl;
    }

    // test that we can project a convex combination point
    {
      const std::array<Point<2>, 2> vertices{
        {face->vertex(0), face->vertex(1)}};
      const std::array<double, 2> weights{{0.125, 0.875}};
      const auto                  vertices_view =
        make_array_view(vertices.begin(), vertices.end());
      const auto weights_view = make_array_view(weights.begin(), weights.end());

      const Point<2> trial_point =
        face->get_manifold().get_new_point(vertices_view, weights_view);
      const Point<2> projected_point =
        GridTools::project_to_object(face, trial_point);
      deallog << std::endl
              << "Project the convex combination of two vertices:" << std::endl
              << "projected point:          " << projected_point << std::endl
              << "Convex combination point: " << trial_point << std::endl
              << "Distance:                 "
              << trial_point.distance(projected_point) << std::endl;
    }

    // easy test: project with spacedim == structdim (should be the identity)
    {
      const Point<2> trial_point{0.417941, 4242424242.4242};
      const Point<2> projected_point =
        GridTools::project_to_object(cell, trial_point);
      deallog << std::endl
              << "Check that projecting with spacedim == structdim" << std::endl
              << "is the identity map:" << std::endl
              << std::endl
              << "trial point:     " << trial_point << std::endl
              << "projected point: " << projected_point << std::endl;
    }

    // harder test: project a vertex onto the geodesic
    {
      const Point<2> trial_point = face->vertex(0);
      const Point<2> projected_point =
        GridTools::project_to_object(cell->face(face_n), trial_point);

      deallog << std::endl
              << "Project a vertex onto the geodesic:" << std::endl
              << std::endl
              << "vertex distance from origin:          "
              << face->vertex(1).distance(Point<2>()) << std::endl
              << "trial point distance from origin:     "
              << trial_point.distance(Point<2>()) << std::endl
              << "projected point distance from origin: "
              << projected_point.distance(Point<2>()) << std::endl
              << "projected point: " << projected_point << std::endl
              << "actual vertex:   " << face->vertex(0) << std::endl;
    }
  }

  // Test project_to_object with a 2D surface in 3D
  {
    TorusManifold<2>    torus_manifold(2.0, 1.0);
    Triangulation<2, 3> triangulation;
    GridGenerator::torus(triangulation, 2.0, 1.0);

    constexpr types::manifold_id torus_id = 42;
    triangulation.set_manifold(torus_id, torus_manifold);
    triangulation.set_all_manifold_ids(torus_id);

    // an ill-posed problem: project a point along the axis of symmetry onto a
    // cell face. Make sure that we end up with something that is on the
    // manifold and is near the equator (for the bottom cells) or on the top
    // (for the top cells).
    const Point<3> trial_point(0.0, 100.0, 0.0);
    deallog << "Test for robustness by projecting points with nonunique"
            << std::endl
            << "minimizers. The output here has been eyeballed as decent."
            << std::endl;

    MappingQGeneric<2, 3> mapping(6);
    for (auto &cell : triangulation.active_cell_iterators())
      {
        const Point<3> projected_point =
          GridTools::project_to_object(cell, trial_point);

        // Ensure that the point we found is both on the manifold and (up to
        // error in the polynomial approximation of the mapping) on the cell.
        Assert((torus_manifold.push_forward(
                  torus_manifold.pull_back(projected_point)) -
                projected_point)
                   .norm() < 1e-14,
               ExcInternalError());
        const Point<2> unit_point =
          mapping.transform_real_to_unit_cell(cell, projected_point);
        Assert(GeometryInfo<2>::is_inside_unit_cell(unit_point, 1.0e-5),
               ExcInternalError());
      }
  }

  // refine in 3D a few times so that we can observe that the projection error
  // drops proportionally as the grid is refined
  for (unsigned int n_refinements = 4; n_refinements < 7; ++n_refinements)
    {
      deallog << "====================================================="
              << std::endl
              << "Number of global refinements: " << n_refinements << std::endl
              << "====================================================="
              << std::endl;

      SphericalManifold<3>         spherical_manifold;
      Triangulation<3>             triangulation;
      constexpr types::manifold_id spherical_id = 42;
      GridGenerator::hyper_shell(triangulation, Point<3>(), 1.0, 2.0);
      triangulation.set_manifold(spherical_id, spherical_manifold);
      triangulation.set_all_manifold_ids(42);

      triangulation.refine_global(n_refinements);

      auto cell = triangulation.begin_active();
      while (!cell->at_boundary())
        ++cell;
      unsigned int face_n = 0;
      while (!cell->face(face_n)->at_boundary())
        ++face_n;
      const auto face = cell->face(face_n);

      // easy test: project the arithmetic mean of the vertices onto the face
      {
        const Point<3> trial_point = 1.1 * face->center();
        const Point<3> projected_point =
          GridTools::project_to_object(face, trial_point);

        deallog << std::endl
                << "Project a reweighed center onto a face in 3D:" << std::endl
                << std::endl
                << "vertex 0 distance from origin: "
                << face->vertex(0).distance(Point<3>()) << std::endl
                << "vertex 1 distance from origin: "
                << face->vertex(1).distance(Point<3>()) << std::endl
                << "trial point distance from origin: "
                << trial_point.distance(Point<3>()) << std::endl
                << "projected point distance from origin: "
                << projected_point.distance(Point<3>()) << std::endl;

        const std::array<Point<3>, 4> vertices{
          {face->vertex(0), face->vertex(1), face->vertex(2), face->vertex(3)}};
        const std::array<double, 4> weights{{0.25, 0.25, 0.25, 0.25}};
        const auto                  vertices_view =
          make_array_view(vertices.begin(), vertices.end());
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());

        const Point<3> true_midpoint =
          face->get_manifold().get_new_point(vertices_view, weights_view);
        deallog << "projected point:    " << projected_point << std::endl
                << "true line midpoint: " << true_midpoint << std::endl
                << "distance less than 5.0e-6: "
                << (projected_point.distance(true_midpoint) < 5.0e-6)
                << std::endl;
      }

      // easy test: project with spacedim == structdim
      {
        const Point<3> trial_point{0.417941, 4242424242.4242, 1};
        const Point<3> projected_point =
          GridTools::project_to_object(cell, trial_point);
        deallog << std::endl
                << "Check that projecting with spacedim == structdim "
                << "is the identity map:" << std::endl
                << std::endl
                << "trial point:     " << trial_point << std::endl
                << "projected point: " << projected_point << std::endl;
      }

      // project a vertex onto the surface
      {
        const Point<3> trial_point = face->vertex(0);
        const Point<3> projected_point =
          GridTools::project_to_object(face, trial_point);

        deallog << std::endl
                << "Project a vertex:" << std::endl
                << std::endl
                << "vertex 0 distance from origin:        "
                << face->vertex(0).distance(Point<3>()) << std::endl
                << "trial point distance from origin:     "
                << trial_point.distance(Point<3>()) << std::endl
                << "projected point distance from origin: "
                << projected_point.distance(Point<3>()) << std::endl
                << "projected point: " << projected_point << std::endl
                << "actual vertex:   " << trial_point << std::endl
                << "distance:        " << trial_point.distance(projected_point)
                << std::endl;
      }

      // project (nearly) a vertex onto the surface
      {
        const std::array<Point<3>, 2> vertices{
          {face->vertex(0), face->vertex(1)}};
        const std::array<double, 2> weights{{1.0 - 1.0 / 2048.0, 1.0 / 2048.0}};
        const auto                  vertices_view =
          make_array_view(vertices.begin(), vertices.end());
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());

        const Point<3> trial_point =
          face->get_manifold().get_new_point(vertices_view, weights_view);
        const Point<3> projected_point =
          GridTools::project_to_object(face, trial_point);

        deallog << std::endl
                << "Project a point near a vertex:" << std::endl
                << std::endl
                << "vertex 0 distance from origin:        "
                << face->vertex(0).distance(Point<3>()) << std::endl
                << "trial point distance from origin:     "
                << trial_point.distance(Point<3>()) << std::endl
                << "projected point distance from origin: "
                << projected_point.distance(Point<3>()) << std::endl
                << "projected point: " << projected_point << std::endl
                << "trial point:     " << trial_point << std::endl
                << "distance:        " << projected_point.distance(trial_point)
                << std::endl;
      }

      // test that we can recover a point that is on the face
      {
        const std::array<Point<3>, 4> vertices{
          {face->vertex(0), face->vertex(1), face->vertex(2), face->vertex(3)}};
        const std::array<double, 4> weights{{0.0625, 0.5, 0.25, 0.1875}};
        const auto                  vertices_view =
          make_array_view(vertices.begin(), vertices.end());
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());

        const Point<3> trial_point =
          face->get_manifold().get_new_point(vertices_view, weights_view);
        const Point<3> projected_point =
          GridTools::project_to_object(face, trial_point);

        deallog << std::endl
                << "Project a weighed face point onto the same face in 3D:"
                << std::endl
                << std::endl
                << "trial point distance from origin:     "
                << trial_point.distance(Point<3>()) << std::endl
                << "projected point distance from origin: "
                << projected_point.distance(Point<3>()) << std::endl
                << "Error is less than 1.0e-3:           "
                << (projected_point.distance(trial_point) < 1.0e-3)
                << std::endl;
      }

      // test that we can recover a point that is moved off the surface along a
      // normal vector
      {
        const std::array<Point<3>, 4> vertices{
          {face->vertex(0), face->vertex(1), face->vertex(2), face->vertex(3)}};
        const std::array<double, 4> weights{{0.0625, 0.5, 0.25, 0.1875}};
        const auto                  vertices_view =
          make_array_view(vertices.begin(), vertices.end());
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());

        Point<3> trial_point =
          face->get_manifold().get_new_point(vertices_view, weights_view);
        Tensor<1, 3> normal_vector = 0.1 * (trial_point - Point<3>());
        trial_point += normal_vector;
        const Point<3> projected_point =
          GridTools::project_to_object(face, trial_point);

        deallog << std::endl
                << "Project a point offset from the face onto the face in 3D:"
                << std::endl
                << std::endl
                << "trial point distance from origin:     "
                << trial_point.distance(Point<3>()) << std::endl
                << "projected point distance from origin: "
                << projected_point.distance(Point<3>()) << std::endl
                << "Error (vs. going along the normal) less than 0.01: "
                << (projected_point.distance(trial_point - normal_vector) <
                    0.01)
                << std::endl;
      }

      // test that we can recover a point that is on a line in 3D
      {
        const auto                    line = face->line(0);
        const std::array<Point<3>, 2> vertices{
          {line->vertex(0), line->vertex(1)}};
        const std::array<double, 2> weights{{0.125, 1.0 - 0.125}};
        const auto                  vertices_view =
          make_array_view(vertices.begin(), vertices.end());
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());

        const Point<3> trial_point =
          face->get_manifold().get_new_point(vertices_view, weights_view);
        const Point<3> projected_point =
          GridTools::project_to_object(face->line(0), trial_point);

        deallog << std::endl
                << "Project a weighed line point onto the same line in 3D:"
                << std::endl
                << std::endl
                << "trial point:                          " << trial_point
                << std::endl
                << "projected point:                      " << projected_point
                << std::endl
                << "trial point distance from origin:     "
                << trial_point.distance(Point<3>()) << std::endl
                << "projected point distance from origin: "
                << projected_point.distance(Point<3>()) << std::endl
                << "Error less than 1.0e-7:               "
                << (projected_point.distance(trial_point) < 1.0e-7)
                << std::endl;
      }
    }
  deallog << "OK" << std::endl;
}
