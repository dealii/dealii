// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// make sure FEValuesExtractors::SymmetricTensor can be default-constructed
// but that it produces an invalid, unusable object

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe, quadrature, update_values);
  fe_values.reinit(dof.begin_active());

  const FEValuesExtractors::SymmetricTensor<2> extr; // invalid object
  try
    {
      fe_values[extr]; // invalid access
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
      goto ok;
    }

  Assert(false, ExcMessage("No exception!?"));

ok:;
}



template <int dim>
void
test()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe(FE_Q<dim>(1), dim * (dim + 1) / 2);
  test(tr, fe);
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog << std::setprecision(2);

  test<1>();
  test<2>();
  test<3>();
}
