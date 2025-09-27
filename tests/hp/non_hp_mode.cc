// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// warn users that they should have nonzero active FE indices in non-hp mode


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  FE_Q<dim> fe(1);

  DoFHandler<dim> dof_handler(triangulation);

  // Choose index out of range here.
  dof_handler.begin_active()->set_active_fe_index(1);

  try
    {
      dof_handler.distribute_dofs(fe);
    }
  catch (ExceptionBase &e)
    {
      deallog << "Successfully failed to distribute DoFs with a single "
              << "nonzero active FE index" << std::endl;
    }
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  test<2>();

  deallog << "OK" << std::endl;
}
