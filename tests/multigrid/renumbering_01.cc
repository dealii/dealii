//----------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2000 - 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// check the new DoFRenumbering::component_wise function that handles
// MGDoFHandlers and renumbers all MG and non-MG dofs

#include "../tests.h"
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_renumbering.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgq.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void check()
{
  FESystem<dim> fe(FE_DGP<dim>(1), 1, FE_DGQ<dim>(2), 2, FE_Q<dim>(3), 1);
  deallog << fe.get_name() << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement ();
  
  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  DoFRenumbering::component_wise (mgdof);

  typename MGDoFHandler<dim>::cell_iterator
    cell = mgdof.begin(),
    endc = mgdof.end();
  std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
  std::vector<unsigned int> mg_dof_indices (fe.dofs_per_cell);
  for (; cell!=endc; ++cell)
    {
      if (!cell->has_children())
	{
	  deallog << "Global numbering: ";
	  cell->get_dof_indices (local_dof_indices);
	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    deallog << local_dof_indices[i] << ' ';
	  deallog << std::endl;
	}
      
      deallog << "MG levelwise numbering: ";
      cell->get_mg_dof_indices (mg_dof_indices);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	deallog << mg_dof_indices[i] << ' ';
      deallog << std::endl;

				       // assert at locally on each
				       // cell that dofs with lower
				       // component also have lower
				       // dof index
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int j=i+1; j<fe.dofs_per_cell; ++j)
	  if (fe.system_to_component_index(i).first <
	      fe.system_to_component_index(j).first)
	    Assert (mg_dof_indices[i] < mg_dof_indices[j],
		    ExcInternalError())
	  else
	    if (fe.system_to_component_index(i).first >
		fe.system_to_component_index(j).first)
	      Assert (mg_dof_indices[i] > mg_dof_indices[j],
		      ExcInternalError());
    }
}


int main()
{
  std::ofstream logfile("renumbering_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}
