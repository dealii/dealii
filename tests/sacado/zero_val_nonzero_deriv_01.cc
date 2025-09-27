// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests the computation of the first derivatives of a function with zero
// input can have a non-zero first derivative.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <Sacado.hpp>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim = 2;
  Triangulation<dim> triangulation;
  FE_Q<dim>          fe(1);
  QGauss<dim>        quadrature_formula(2);
  DoFHandler<dim>    dof_handler(triangulation);
  Vector<double>     zero_solution;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values);

  GridGenerator::hyper_cube(triangulation, -1, 1);
  dof_handler.distribute_dofs(fe);
  zero_solution.reinit(dof_handler.n_dofs());

  const FEValuesExtractors::Scalar extractor_sclr(0);

  using ad_type                              = Sacado::Fad::DFad<double>;
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
  DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      // Extract (zero'd) DoF values
      Vector<double> local_dof_values(cell->get_fe().dofs_per_cell);
      cell->get_dof_values(zero_solution, local_dof_values);

      // Configure AD-equivalent DoF values. These
      // are the independent variables "x".
      std::vector<ad_type> local_ad_dof_values(fe.n_dofs_per_cell());
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        local_ad_dof_values[i] =
          ad_type(fe.n_dofs_per_cell(), i, local_dof_values[i]);

      // Get solution values while tracking the influence
      // that perturbing the DoF values may have on the
      // derivative of the dependent variable "f(x)".
      std::vector<ad_type> local_qp_soln_values_ad(fe.n_dofs_per_cell());
      fe_values[extractor_sclr].get_function_values_from_local_dof_values(
        local_ad_dof_values, local_qp_soln_values_ad);

      // Compute the dependent variable "f(x)".
      ad_type f_x = 0.0;
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        f_x += local_qp_soln_values_ad[i];

      deallog << "f(x): " << f_x.val() << std::endl;
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        deallog << "df/dx_" << i << ": " << f_x.dx(i) << std::endl;
    }

  deallog << "OK" << std::endl;
}
