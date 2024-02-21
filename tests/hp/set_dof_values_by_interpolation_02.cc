// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// cell->set_dof_values_by_interpolation can not work properly in the hp-
// context when called on non-active cells because these do not have a
// finite element associated with them
//
// this test verifies that if we call the function on active cells with no
// explicitly given fe_index that we get the same result as from
// cell->set_dof_values, and that if we call it with an fe_index for a Q1
// element that we simply get the vertex dof values set on the space of that
// cell

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 2



template <int dim>
void
test()
{
  // create a hp::DoFHandler with different finite elements on the
  // cells. note that we skip setting active_fe_indices on inactive
  // elements
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  tr.refine_global(2);

  hp::FECollection<dim> fe;
  for (unsigned int i = 1; i < 5; ++i)
    fe.push_back(FE_Q<dim>(i));

  DoFHandler<dim> dof_handler(tr);
  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    if (cell->has_children() == false)
      cell->set_active_fe_index(cell->index() % fe.size());

  dof_handler.distribute_dofs(fe);

  // create a mostly arbitrary FE field
  Vector<double> solution1(dof_handler.n_dofs());
  Vector<double> solution2(dof_handler.n_dofs());

  // do the test
  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      solution1 = 0;
      solution2 = 0;

      // set values without specifying an explicit fe_index
      Vector<double> local(cell->get_fe().dofs_per_cell);
      for (unsigned int i = 0; i < local.size(); ++i)
        local(i) = i;

      cell->set_dof_values_by_interpolation(local, solution1);

      // then do the same with the "correct", local fe_index
      cell->set_dof_values_by_interpolation(local,
                                            solution2,
                                            cell->active_fe_index());

      // now verify correctness
      AssertThrow(solution1 == solution2, ExcInternalError());
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(PRECISION);
  deallog.get_file_stream().setf(std::ios::fixed);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
