// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This test verifies that the transfinite interpolation works on a torus

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <memory>

#include "../tests.h"


template <int dim>
class GradingManifold : public ChartManifold<dim>
{
public:
  GradingManifold(const Point<dim> center,
                  const double     grading,
                  const double     epsilon)
    : center(center)
    , grading(grading)
    , epsilon(epsilon)
    , polar_manifold(center)
  {}

  virtual Point<dim>
  pull_back(const Point<dim> &space_point) const final override
  {
    auto point = polar_manifold.pull_back(space_point);
    Assert(point[0] >= 0., ExcInternalError());
    point[0] = std::pow(point[0] + epsilon, 1. / grading) -
               std::pow(epsilon, 1. / grading) + 1.e-14;
    const auto chart_point = polar_manifold.push_forward(point);
    return chart_point;
  }

  virtual Point<dim>
  push_forward(const Point<dim> &chart_point) const final override
  {
    auto point = polar_manifold.pull_back(chart_point);
    point[0]   = std::pow(point[0] + std::pow(epsilon, 1. / grading), grading) -
               epsilon + 1.e-14;
    Assert(point[0] >= 0., ExcInternalError());
    return polar_manifold.push_forward(point);
  }

  std::unique_ptr<Manifold<dim, dim>>
  clone() const final override
  {
    return std::make_unique<GradingManifold<dim>>(center, grading, epsilon);
  }

private:
  const Point<dim> center;
  const double     grading;
  const double     epsilon;

  PolarManifold<dim> polar_manifold;
};


int
main()
{
  initlog();
  deallog << std::setprecision(9);

  const std::vector<Point<2>> vertices{
    {0.5 * 2.5, -std::sqrt(3.) / 2. * 2.5}, // 0
    {-0.9 + 0.7, -0.025},                   // 1
    {-0.9 + 0.7, 0.025},                    // 2
    {0.5 * 2.5, std::sqrt(3.) / 2. * 2.5},  // 3
    {2.5, -0.5 * 2.5},                      // 4
    {2.5, -0.416667},                       // 5
    {2.5, 0.416667},                        // 6
    {2.5, 0.5 * 2.5},                       // 7
  };

  std::vector<CellData<2>> cells(3);
  cells[0].vertices = {0, 4, 1, 5};
  cells[1].vertices = {1, 5, 2, 6};
  cells[2].vertices = {2, 6, 3, 7};

  Triangulation<2> triangulation;
  triangulation.create_triangulation(vertices, cells, SubCellData());

  {
    triangulation.set_all_manifold_ids(3);
    const auto center_cell  = std::next(triangulation.begin_active());
    const auto lower_radial = center_cell->face(2);
    const auto upper_radial = center_cell->face(3);
    lower_radial->set_manifold_id(1);
    upper_radial->set_manifold_id(2);
  }

  GradingManifold<2> lower_radial{Point<2>{-0.9 + 0.7, -0.025}, 3.0, 0.01};
  triangulation.set_manifold(1, lower_radial);

  GradingManifold<2> upper_radial{Point<2>{-0.9 + 0.7, 0.025}, 3.0, 0.01};
  triangulation.set_manifold(2, upper_radial);

  TransfiniteInterpolationManifold<2> transfinite;
  transfinite.initialize(triangulation);
  triangulation.set_manifold(3, transfinite);

  /*
   * triangulation.refine_global(8) fails with:
   *
   *   An error occurred in line [...] in function
   *       [...]compute_chart_points(const ArrayView<const Point<spacedim>> &,
   *   ArrayView<Point<dim>>) The violated condition was: false Additional
   *   information: Looking at cell 1_0: with vertices: -0.2 -0.025    2.5
   *   -0.416667    -0.2 0.025    2.5 0.416667 Transformation to chart
   *   coordinates: -0.1271 -0.0177875 -> 0.1875 0.25 -0.1271 -0.015564 -> 20 20
   *   Looking at cell 0_0: with vertices:
   *   1.25 -2.16506    2.5 -1.25    -0.2 -0.025    2.5 -0.416667
   *   Transformation to chart coordinates:
   *   -0.1271 -0.0177875 -> 0.211251 1.01047
   *   -0.1271 -0.015564 -> 0.214091 1.01181
   *   Looking at cell 2_0: with vertices:
   *   -0.2 0.025    2.5 0.416667    1.25 2.16506    2.5 1.25
   *   Transformation to chart coordinates:
   *   -0.1271 -0.0177875 -> 0.253794 -0.0323121
   *   -0.1271 -0.015564 -> 0.251301 -0.0309146
   *
   * The following is a simplified version that fails with a slightly
   * different error message.
   */

  std::vector<Point<2>> points{Point<2>{-0.1271, -0.0177875},
                               Point<2>{-0.1271, -0.015564}};
  ArrayView<Point<2>>   points_view{points.data(), 2};

  std::vector<double> weights{0.5, 0.5};
  ArrayView<double>   weights_view{weights.data(), 2};

  deallog << triangulation.get_manifold(3).get_new_point(points_view,
                                                         weights_view)
          << std::endl;

  return 0;
}
