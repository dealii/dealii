//----------------------------  dof_renumbering.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_renumbering.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <numerics/dof_renumbering.h>

#include <fstream>



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  std::vector<unsigned int> v (dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      deallog << "Cell " << cell << std::endl;
      cell->get_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
	deallog << v[i] << ' ';
      deallog << std::endl;
    }
}



template <int dim>
void
check ()
{
  CosineFunction<dim> cosine;
  
  Triangulation<dim> tr;  
  if (dim==2)
    GridGenerator::hyper_ball(tr);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);
  
  FESystem<dim> element (FE_Q<dim>(3), 2, FE_DGQ<dim>(1), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  print_dofs (dof);
  
  DoFRenumbering::Cuthill_McKee (dof);              print_dofs (dof);
  DoFRenumbering::Cuthill_McKee (dof, true);        print_dofs (dof);
  DoFRenumbering::Cuthill_McKee (dof, true, true);  print_dofs (dof);

  std::vector<unsigned int> order(element.n_components());
  for (unsigned int i=0; i<order.size(); ++i) order[i] = order.size()-i-1;
  DoFRenumbering::component_wise (dof, order);      print_dofs (dof);

  std::vector<bool> selected_dofs (dof.n_dofs(), false);
  for (unsigned int i=0; i<dof.n_dofs(); ++i) if (i%2==0) selected_dofs[i] = true;
  DoFRenumbering::sort_selected_dofs_back (dof, selected_dofs);
  print_dofs (dof);
}


int main ()
{
  ofstream logfile ("dof_renumbering.output");
  logfile.precision (2);
  logfile.setf(ios::fixed);  
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
