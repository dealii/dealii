// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// there used to be a bug in MappingCartesian, where the inverse jacobian
// was not computed correctly. Check that all works as intended.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check(const Triangulation<dim> &tria)
{
  MappingCartesian<dim> mapping;
  FE_Q<dim>             fe(1);
  DoFHandler<dim>       dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<dim> quadrature(1);

  UpdateFlags update_flags =
    update_quadrature_points | update_JxW_values | update_jacobians |
    update_jacobian_grads | update_jacobian_pushed_forward_grads |
    update_jacobian_2nd_derivatives |
    update_jacobian_pushed_forward_2nd_derivatives |
    update_jacobian_3rd_derivatives |
    update_jacobian_pushed_forward_3rd_derivatives | update_inverse_jacobians |
    // Transformation dependence
    update_covariant_transformation | update_contravariant_transformation |
    update_transformation_values | update_transformation_gradients |
    // Volume data
    update_volume_elements;

  FEValues<dim> fe_values(mapping, fe, quadrature, update_flags);

  fe_values.reinit(dof_handler.begin_active());

  for (unsigned int i = 0; i < dim; ++i)
    deallog << fe_values.inverse_jacobian(0)[i] << std::endl;

  deallog << std::endl;

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_cube(coarse_grid);
    check(coarse_grid);
  }
}
