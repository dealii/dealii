// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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


// cell->get_interpolated_dof_values can not work properly in the hp
// context when called on non-active cells because these do not have a
// finite element associated with them
//
// this test verifies that if we call the function on active cells with no
// explicitly given fe_index that we get the same result as from
// cell->get_dof_values, and that if we call it with an fe_index for a Q1
// element that we simply get the vertex dof values

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_q.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2




template <int dim>
void test ()
{
  // create a hp::DoFHandler with different finite elements on the
  // cells. note that we skip setting active_fe_indices on inactive
  // elements
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  tr.refine_global (2);

  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<5; ++i)
    fe.push_back (FE_Q<dim>(i));

  hp::DoFHandler<dim> dof_handler(tr);
  for (typename hp::DoFHandler<dim>::cell_iterator cell=dof_handler.begin();
       cell!=dof_handler.end(); ++cell)
    if (cell->has_children() == false)
      cell->set_active_fe_index (cell->index() % fe.size());

  dof_handler.distribute_dofs (fe);

  // create a mostly arbitrary FE field
  Vector<double> solution(dof_handler.n_dofs());
  for (unsigned int i=0; i<solution.size(); ++i)
    solution(i) = i;

  // do the test
  for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      // get values without specifying an explicit fe_index
      Vector<double> local1 (cell->get_fe().dofs_per_cell);
      cell->get_interpolated_dof_values (solution, local1);

      // then do the same with the "correct", local fe_index
      Vector<double> local2 (cell->get_fe().dofs_per_cell);
      cell->get_interpolated_dof_values (solution, local2,
					 cell->active_fe_index());

      // and do it a third time with the fe_index for a Q1 element
      Vector<double> local3 (fe[0].dofs_per_cell);
      cell->get_interpolated_dof_values (solution, local3,
					 0);

      // now verify correctness
      Assert (local1 == local2, ExcInternalError());

      // also for the second test. note that vertex dofs always come first in
      // local1, so we can easily compare
      for (unsigned int i=0; i<fe[0].dofs_per_cell; ++i)
	Assert (local1[i] == local3[i],
		ExcInternalError());
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}



