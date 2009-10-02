//----------------------------  renumber_component_wise.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  renumber_component_wise.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include "../tests.h"
#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <hp/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_system.h>
#include <hp/fe_collection.h>

#include <fstream>



template <int dim>
void
print_dofs (const hp::DoFHandler<dim> &dof)
{
  std::vector<unsigned int> v;
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
check_renumbering(hp::DoFHandler<dim>& dof)
{
				   // Prepare a reordering of
				   // components for later use
  std::vector<unsigned int> order(dof.get_fe().n_components());
  for (unsigned int i=0; i<order.size(); ++i)
    order[i] = order.size()-i-1;

  DoFRenumbering::component_wise (dof, order);
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
  std::ofstream logfile ("renumber_component_wise/output");
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
