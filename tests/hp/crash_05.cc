//----------------------------  crash_05.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  crash_05.cc  ---------------------------


// a crash when using different fe's on different cells


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (1);

				   // use the same element, but with
				   // different indices
  hp::FECollection<dim> fe_collection;
  for (unsigned int i=0; i<tria.n_active_cells(); ++i)
    fe_collection.push_back(FE_Q<dim> (1));

  hp::DoFHandler<dim> dof_handler(tria);

  unsigned int fe_index = 0;
  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell, ++fe_index)
    {
      deallog << "Setting fe_index=" << fe_index << " on cell " << cell
	      << std::endl;
      cell->set_active_fe_index (fe_index);
    }
  
  dof_handler.distribute_dofs(fe_collection);

  std::vector<unsigned int> local_dof_indices;
  for (typename hp::DoFHandler<dim>::active_cell_iterator
         cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      local_dof_indices.resize (cell->get_fe().dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      deallog << cell << std::endl;
      for (unsigned int i=0; i<local_dof_indices.size(); ++i)
        deallog << local_dof_indices[i] << ' ';
      deallog << std::endl;
    }
}


int main ()
{
  std::ofstream logfile("crash_05/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
