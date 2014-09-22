// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// we used to be able to call DoFCellAccessor::get_dof_indices also for
// inactive cells, check that this is now forbidden.

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>
#include <vector>





template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);

  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (1);

  dof_handler.distribute_dofs (fe);

  // loop over all cells, active or
  // not
  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin();
       cell != dof_handler.end(); ++cell)
    {
      try
        {
          cell->get_dof_indices (local_dof_indices);
        }
      catch (...)
        {
          deallog << "Assertion: cell not active." << std::endl;
          continue;
        }


      deallog << "Cell = " << cell
              << ", DoFs=";
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        {
          Assert (local_dof_indices[i] != DoFHandler<dim>::invalid_dof_index,
                  ExcInternalError());
          deallog << local_dof_indices[i] << ' ';
        }

      deallog << std::endl;
    }
}



int main ()
{
  deal_II_exceptions::disable_abort_on_exception();
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}

