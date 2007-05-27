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

#include <iomanip>

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

  deallog << "Cells " << std::setw(12)  << tr.n_cells()
	  << " active " << std::setw(12)  << tr.n_active_cells()
	  << " memory " << std::setw(12)  << tr.memory_consumption()
	  << " quotient " << (1./tr.n_cells()*tr.memory_consumption())
	  << std::endl;

  FE_Q<dim> q1(1);
  FE_Q<dim> q3(3);
  FESystem<dim> sys1(q3, 1);
  FESystem<dim> sys2(q3, 10);
  DoFHandler<dim> dof(tr);
  
  dof.distribute_dofs(q1);
  deallog << "Dofs Q1   " << std::setw(12) << dof.n_dofs()
	  << " memory " << std::setw(12) << dof.memory_consumption()
	  << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
	  << std::endl;
  
  dof.distribute_dofs(q3);
  deallog << "Dofs Q3   " << std::setw(12) << dof.n_dofs()
	  << " memory " << std::setw(12) << dof.memory_consumption()
	  << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
	  << std::endl;
  
  dof.distribute_dofs(sys1);
  deallog << "Dofs Sys1 " << std::setw(12) << dof.n_dofs()
	  << " memory " << std::setw(12) << dof.memory_consumption()
	  << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
	  << std::endl;  
  
  dof.distribute_dofs(sys2);
  deallog << "Dofs Sys2 " << std::setw(12) << dof.n_dofs()
	  << " memory " << std::setw(12) << dof.memory_consumption()
	  << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
	  << std::endl;  
}

int main()
{
  deallog.log_execution_time(true);
  check<2>();
  check<3>();
}

