//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>

using namespace dealii;


template <int dim>
void check ()
{
  deallog << "Dimension " << dim << std::endl;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(18/dim);

  deallog << "Cells " << tr.n_cells()
	  << " active " << tr.n_active_cells()
	  << " memory " << tr.memory_consumption()
	  << std::endl;

  FE_Q<dim> q1(1);
  FE_Q<dim> q3(3);
  FESystem<dim> sys1(q3, 1);
  FESystem<dim> sys2(q3, 10);
  DoFHandler<dim> dof(tr);
  
  dof.distribute_dofs(q1);
  deallog << "Dofs Q1 " << dof.n_dofs()
	  << " memory " << dof.memory_consumption()
	  << std::endl;
  
  dof.distribute_dofs(q3);
  deallog << "Dofs Q3 " << dof.n_dofs()
	  << " memory " << dof.memory_consumption()
	  << std::endl;
  
  dof.distribute_dofs(sys1);
  deallog << "Dofs Sys1 " << dof.n_dofs()
	  << " memory " << dof.memory_consumption()
	  << std::endl;  
  
  dof.distribute_dofs(sys2);
  deallog << "Dofs Sys2 " << dof.n_dofs()
	  << " memory " << dof.memory_consumption()
	  << std::endl;  
}

int main()
{
  deallog.log_execution_time(true);
  check<2>();
  check<3>();
}

