//----------------------------  dof_renumbering.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003 by the deal.II authors
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
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>

#include <fstream>



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  std::vector<unsigned int> v (dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
	deallog << v[i] << ' ';
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs (const MGDoFHandler<dim> &dof, unsigned int level)
{
  std::vector<unsigned int> v (dof.get_fe().dofs_per_cell);
  for (typename MGDoFHandler<dim>::cell_iterator cell=dof.begin(level);
       cell != dof.end(level); ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_mg_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
	deallog << v[i] << ' ';
      deallog << std::endl;
    }
}



template <int dim>
void
check ()
{
  Functions::CosineFunction<dim> cosine;
  
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
  
  FESystem<dim> element (FE_Q<dim>(2), 2, FE_DGQ<dim>(1), 1);
  MGDoFHandler<dim> mgdof(tr);
  DoFHandler<dim>& dof = mgdof;
  dof.distribute_dofs(element);

				   // Prepare a reordering of
				   // components for later use
  std::vector<unsigned int> order(element.n_components());
  for (unsigned int i=0; i<order.size(); ++i) order[i] = order.size()-i-1;

				   // Check global ordering
  print_dofs (dof);
  
  DoFRenumbering::Cuthill_McKee (dof);
  print_dofs (dof);
  DoFRenumbering::Cuthill_McKee (dof, true);
  print_dofs (dof);
  DoFRenumbering::Cuthill_McKee (dof, true, true);
  print_dofs (dof);

  DoFRenumbering::component_wise (dof, order);
  print_dofs (dof);

  std::vector<bool> selected_dofs (dof.n_dofs(), false);
  for (unsigned int i=0; i<dof.n_dofs(); ++i) if (i%2==0) selected_dofs[i] = true;
  DoFRenumbering::sort_selected_dofs_back (dof, selected_dofs);
  print_dofs (dof);

				   // Check level ordering
  for (unsigned int level=0;level<tr.n_levels();++level)
    {
      print_dofs (mgdof, level);

// Reinsert after fixing
//        DoFRenumbering::Cuthill_McKee (mgdof, level);
//        print_dofs (mgdof, level);
//        DoFRenumbering::Cuthill_McKee (mgdof, level, true);
//        print_dofs (mgdof, level);

      DoFRenumbering::component_wise (mgdof, order);
      print_dofs (mgdof, level);
    }
}


int main ()
{
  std::ofstream logfile ("dof_renumbering.output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
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
