//----------------------------  crash_16.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  crash_16.cc  ---------------------------


// an extract of hp_constraints_q_system_x_01. something still goes wrong,
// akin to what happens in crash_15. (it turned out to be bogus index
// computations.)

char logname[] = "crash_16/output";


#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>

#include <grid/tria.h>
#include <hp/dof_handler.h>
#include <lac/constraint_matrix.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/vectors.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>

#include <fstream>
#include <vector>



template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
				  FE_DGQ<dim>(j), 1));

  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  hp::DoFHandler<dim>        dof_handler(triangulation);

				   // distribute fe_indices randomly
  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (rand() % fe.size());
  dof_handler.distribute_dofs (fe);

				   // loop over all lines and make sure that
				   // all the DoF indices on these lines are
				   // identical
  std::vector<unsigned int> indices_1;
  std::vector<unsigned int> indices_2;
  
  for (typename hp::DoFHandler<dim>::active_line_iterator
	 line = dof_handler.begin_active_line();
       line != dof_handler.end_line(); ++line)
    {
      deallog << "line=" << line << std::endl;
      
      for (unsigned int f=0; f<line->n_active_fe_indices(); ++f)
	{
	  indices_1.resize (fe[line->nth_active_fe_index(f)].dofs_per_line +
			    2 * fe[line->nth_active_fe_index(f)].dofs_per_vertex);
	  line->get_dof_indices (indices_1,
				 line->nth_active_fe_index(f));
	  
	  deallog << "  fe index=" << line->nth_active_fe_index(f)
		  << ", indices=";
	  for (unsigned int i=0; i<indices_1.size(); ++i)
	    deallog << indices_1[i] << ' ';
	  
	  deallog << std::endl;
	}

      for (unsigned int f=0; f<line->n_active_fe_indices(); ++f)
	{
	  indices_1.resize (fe[line->nth_active_fe_index(f)].dofs_per_line +
			    2 * fe[line->nth_active_fe_index(f)].dofs_per_vertex);
	  line->get_dof_indices (indices_1,
				 line->nth_active_fe_index(f));
	  for (unsigned int g=f+1; g<line->n_active_fe_indices(); ++g)
	    if (fe[line->nth_active_fe_index(f)].dofs_per_line
		==
		fe[line->nth_active_fe_index(g)].dofs_per_line)
	      {
		indices_2.resize (fe[line->nth_active_fe_index(g)].dofs_per_line +
				  2 * fe[line->nth_active_fe_index(g)].dofs_per_vertex);
		line->get_dof_indices (indices_2,
				       line->nth_active_fe_index(g));
		Assert (indices_1 == indices_2,
			ExcInternalError());
	      }
	}
    }
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<3> ();
}

