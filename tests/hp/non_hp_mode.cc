// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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
