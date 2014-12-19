// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check the case where the two children on a face with a hanging node have
// different FENothing structure. compared to _12 and _13, the situation is
// even more wicked since in the picture below 0=FE_Q(1), 1=FE_Nothing,
// 2=FE_Q(2), i.e. 2 dominates 1

// active_fe_index
// *-----*
// |     |
// |  0  |
// |     |
// *--*--*
// | 1| 2|
// *--*--*
// | 1| 2|
// *--*--*

// DoF index
// 2--7--3
// |     |
// 4  8  5
// |     |
// 0--6--1
// |  |  |
// *--11-12
// |  |  |
// *--9--10


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/vector_tools.h>


#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  std::vector<unsigned int> sub(dim, 1);
  sub[dim-1] = 2;

  GridGenerator::subdivided_hyper_rectangle (triangulation, sub, Point<dim>(), Point<dim>());
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  hp::FECollection<dim>    fe_collection;
  fe_collection.push_back (FE_Q<dim>(1));
  fe_collection.push_back (FE_Nothing<dim>());
  fe_collection.push_back (FE_Q<dim>(2));

  hp::DoFHandler<dim>      dof_handler (triangulation);

  dof_handler.begin_active()->set_active_fe_index(1);
  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(0),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->index() % 2 == 0)
      cell->set_active_fe_index (1);
    else
      cell->set_active_fe_index (2);

  dof_handler.distribute_dofs (fe_collection);

  ConstraintMatrix constraints;

  DoFTools::make_hanging_node_constraints (dof_handler, constraints);

  constraints.close();

  deallog << "   Number of constraints:        "
          << constraints.n_constraints()
          << std::endl;
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell != endc; cell++)
      {
        deallog << cell << ' ' << cell->active_fe_index() << std::endl
                << "   ";
        std::vector<types::global_dof_index> local_dof_indices (cell->get_fe().dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
          deallog << local_dof_indices[i]
                  << (constraints.is_constrained(local_dof_indices[i]) ?
                      "*" : "")
                  << ' ';
        deallog << std::endl;
      }
  }

  constraints.print (deallog.get_file_stream());
}



int main ()
{
  std::ofstream logfile("output");
  logfile.precision(3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();

  deallog << "OK" << std::endl;
}
