//----------------------------  get_dof_indices_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  get_dof_indices_01.cc  ---------------------------


// make sure we can call DoFCellAccessor::get_dof_indices also for
// inactive cells

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
  std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin();
       cell != dof_handler.end(); ++cell)
    {
      cell->get_dof_indices (local_dof_indices);

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
  std::ofstream logfile("get_dof_indices_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
  
  return 0;
}

