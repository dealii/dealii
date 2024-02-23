// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// the Hessian of the RT element was a tensor of NaN's. This doesn't make much
// sense

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe,
                          quadrature,
                          update_covariant_transformation | update_hessians);

  fe_values.reinit(dof.begin_active());

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << fe_values.shape_hessian_component(0, 0, 0)[i][j] << std::endl;

  // compare the hessian with
  // itself. this fails if the
  // values are NaN's which the
  // Hessian consists of at the
  // time this test is written
  Assert(fe_values.shape_hessian_component(0, 0, 0) ==
           fe_values.shape_hessian_component(0, 0, 0),
         ExcInternalError());
}



template <int dim>
void
test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  test(tr, FE_RaviartThomas<dim>(1));
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
