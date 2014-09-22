// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// an extract of the _04 test that produced different (integer)
// results on different machines, inexplicably


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> v (dof.get_fe().dofs_per_cell);

  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      cell->get_dof_indices (v);
      deallog << "cell=" << cell << std::endl;
      for (unsigned int i=0; i<v.size(); ++i)
        deallog << v[i] << std::endl;
    }
}




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1., 1.);

  FESystem<dim> fe(FE_Q<dim>(2),2);
  DoFHandler<dim> dof(tr);

  dof.distribute_dofs(fe);
  std::vector<types::global_dof_index> new_dofs (dof.n_dofs());
  DoFRenumbering::boost::compute_Cuthill_McKee(new_dofs, dof);

  for (unsigned int i=0; i<new_dofs.size(); ++i)
    deallog << new_dofs[i] << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  check<2> ();
}
