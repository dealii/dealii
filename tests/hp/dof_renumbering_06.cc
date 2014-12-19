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


// Check DoFRenumbering::boost::minimum_degree


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
print_dofs (const hp::DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> v;
  for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      v.resize (cell->get_fe().dofs_per_cell);
      deallog << "Cell " << cell << " -- ";
      cell->get_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
        deallog << v[i] << ' ';
      deallog << std::endl;
    }
}



template <int dim>
void
check_renumbering(hp::DoFHandler<dim> &dof)
{
  for (unsigned int i=0; i<dof.get_fe().size(); ++i)
    {
      const FiniteElement<dim> &element = dof.get_fe()[i];
      deallog << element.get_name() << std::endl;
    }

  DoFRenumbering::boost::minimum_degree(dof);
  print_dofs (dof);
}


template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  hp::DoFHandler<dim> dof(tr);
  {
    bool coin = false;
    for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
         cell != dof.end(); ++cell)
      {
        cell->set_active_fe_index (coin ? 0 : 1);
        coin = !coin;
      }
  }

  FESystem<dim> e1 (FE_Q<dim>(2), 2, FE_DGQ<dim>(1), 2);
  FESystem<dim> e2 (FE_Q<dim>(1), 2, FE_DGQ<dim>(2), 2);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (e1);
  fe_collection.push_back (e2);

  dof.distribute_dofs(fe_collection);
  check_renumbering(dof);
  dof.clear();
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
