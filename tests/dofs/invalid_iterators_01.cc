// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that assigning one invalid iterator to another works. this
// test is for ::DoFHandler

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
check()
{
  typename DoFHandler<dim>::active_cell_iterator invalid_1;

  // try copy constructor
  typename DoFHandler<dim>::active_cell_iterator invalid_2 = invalid_1;

  // now also try copy operator
  invalid_1 = invalid_2;

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
