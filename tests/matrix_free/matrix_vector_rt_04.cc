// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test is the same as matrix_vector_rt_01.cc but with non-affine cells in
// standard orientation.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

#include "matrix_vector_rt_common.h"

// This class is taken from
// https://github.com/exadg/exadg/blob/master/include/exadg/grid/deformed_cube_manifold.h
template <int dim>
class DeformedCubeManifold : public dealii::ChartManifold<dim, dim, dim>
{
public:
  DeformedCubeManifold(const double       left,
                       const double       right,
                       const double       deformation,
                       const unsigned int frequency = 1)
    : left(left)
    , right(right)
    , deformation(deformation)
    , frequency(frequency)
  {}

  dealii::Point<dim>
  push_forward(const dealii::Point<dim> &chart_point) const override
  {
    double sinval = deformation;
    for (unsigned int d = 0; d < dim; ++d)
      sinval *= std::sin(frequency * dealii::numbers::PI *
                         (chart_point[d] - left) / (right - left));
    dealii::Point<dim> space_point;
    for (unsigned int d = 0; d < dim; ++d)
      space_point[d] = chart_point[d] + sinval;
    return space_point;
  }

  dealii::Point<dim>
  pull_back(const dealii::Point<dim> &space_point) const override
  {
    dealii::Point<dim> x = space_point;
    dealii::Point<dim> one;
    for (unsigned int d = 0; d < dim; ++d)
      one[d] = 1.;

    // Newton iteration to solve the nonlinear equation given by the point
    dealii::Tensor<1, dim> sinvals;
    for (unsigned int d = 0; d < dim; ++d)
      sinvals[d] = std::sin(frequency * dealii::numbers::PI * (x[d] - left) /
                            (right - left));

    double sinval = deformation;
    for (unsigned int d = 0; d < dim; ++d)
      sinval *= sinvals[d];
    dealii::Tensor<1, dim> residual = space_point - x - sinval * one;
    unsigned int           its      = 0;
    while (residual.norm() > 1e-12 && its < 100)
      {
        dealii::Tensor<2, dim> jacobian;
        for (unsigned int d = 0; d < dim; ++d)
          jacobian[d][d] = 1.;
        for (unsigned int d = 0; d < dim; ++d)
          {
            double sinval_der = deformation * frequency / (right - left) *
                                dealii::numbers::PI *
                                std::cos(frequency * dealii::numbers::PI *
                                         (x[d] - left) / (right - left));
            for (unsigned int e = 0; e < dim; ++e)
              if (e != d)
                sinval_der *= sinvals[e];
            for (unsigned int e = 0; e < dim; ++e)
              jacobian[e][d] += sinval_der;
          }

        x += invert(jacobian) * residual;

        for (unsigned int d = 0; d < dim; ++d)
          sinvals[d] = std::sin(frequency * dealii::numbers::PI *
                                (x[d] - left) / (right - left));

        sinval = deformation;
        for (unsigned int d = 0; d < dim; ++d)
          sinval *= sinvals[d];
        residual = space_point - x - sinval * one;
        ++its;
      }
    AssertThrow(residual.norm() < 1e-12,
                dealii::ExcMessage("Newton for point did not converge."));
    return x;
  }

  std::unique_ptr<dealii::Manifold<dim>>
  clone() const override
  {
    return std::make_unique<DeformedCubeManifold<dim>>(left,
                                                       right,
                                                       deformation,
                                                       frequency);
  }

private:
  const double       left;
  const double       right;
  const double       deformation;
  const unsigned int frequency;
};

template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  const unsigned int        frequency   = 2;
  const double              deformation = 0.05;
  DeformedCubeManifold<dim> manifold(0.0, 1.0, deformation, frequency);
  tria.set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  std::vector<bool> vertex_touched(tria.n_vertices(), false);

  for (auto cell : tria.cell_iterators())
    {
      for (const auto &v : cell->vertex_indices())
        {
          if (vertex_touched[cell->vertex_index(v)] == false)
            {
              Point<dim> &vertex    = cell->vertex(v);
              Point<dim>  new_point = manifold.push_forward(vertex);
              vertex                = new_point;
              vertex_touched[cell->vertex_index(v)] = true;
            }
        }
    }


  FE_RaviartThomasNodal<dim> fe(fe_degree - 1);
  DoFHandler<dim>            dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  deallog << "Using " << dof.get_fe().get_name() << std::endl;
  deallog << "Number of cells: " << dof.get_triangulation().n_active_cells()
          << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl
          << std::endl;
  do_test<dim, fe_degree, double>(dof, constraints, TestType::values);
  do_test<dim, fe_degree, double>(dof, constraints, TestType::gradients);
  do_test<dim, fe_degree, double>(dof, constraints, TestType::divergence);
}
