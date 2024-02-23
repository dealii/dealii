// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// an extract of _09 that failed at the time of writing the test
// with an internal error


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_values.h>

#include "../tests.h"



template <int dim>
void
test()
{
  FESystem<dim> fe(FE_Nothing<dim>(), 2);
  FEValues<dim> fe_values(fe, QGauss<dim>(2), update_values);

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
