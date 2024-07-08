// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the case where the two children on a face with a hanging node have
// different FE_Nothing structure. while all cells in the _12 test had degrees
// of freedom, here, cells with active_fe_index==1 are simply FE_Nothing:

// active_fe_index
// *-----*
// |     |
// |  0  |
// |     |
// *--*--*
// | 1| 0|
// *--*--*
// | 1| 0|
// *--*--*


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim>        triangulation;
  std::vector<unsigned int> sub(dim, 1);
  sub[dim - 1] = 2;

  Point<2> p1;
  Point<2> p2(1.0, 2.0);
  GridGenerator::subdivided_hyper_rectangle(triangulation, sub, p1, p2);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Nothing<dim>());

  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.begin_active()->set_active_fe_index(1);
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(0),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    if (cell->index() % 2 == 0)
      cell->set_active_fe_index(1);
    else
      cell->set_active_fe_index(0);

  dof_handler.distribute_dofs(fe_collection);

  AffineConstraints<double> constraints;

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  constraints.close();

  deallog << "   Number of constraints:        " << constraints.n_constraints()
          << std::endl;
  {
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();

    for (; cell != endc; ++cell)
      {
        deallog << cell << ' ' << cell->active_fe_index() << std::endl << "   ";
        std::vector<types::global_dof_index> local_dof_indices(
          cell->get_fe().dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < cell->get_fe().dofs_per_cell; ++i)
          deallog << local_dof_indices[i]
                  << (constraints.is_constrained(local_dof_indices[i]) ? "*" :
                                                                         "")
                  << ' ';
        deallog << std::endl;
      }
  }
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<2>();

  deallog << "OK" << std::endl;
}
