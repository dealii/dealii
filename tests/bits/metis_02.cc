//----------------------------  metis_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  metis_02.cc  ---------------------------


// check GridTools::partition_triangulation. generate some output

#include "../tests.h"
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_q.h>
#include <numerics/data_out.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;
  
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (4-dim);
  for (unsigned int i=0; i<11-2*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
      for (unsigned int index=0; cell != triangulation.end(); ++cell, ++index)
        if (index % (3*dim) == 0)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement ();
    }
  

                                   // subdivide into 5 subdomains
  GridTools::partition_triangulation (5, triangulation);

                                   // generate a field where the value equals
                                   // the subdomain number, and output it
  FE_Q<dim> fe(2);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  
  Vector<double> partitions (triangulation.n_active_cells());
  {
    typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
    for (unsigned int index=0; cell != triangulation.end(); ++cell, ++index)
      partitions(index) = cell->subdomain_id();
  }
  Vector<double> dof_part (dof_handler.n_dofs());
  std::vector<unsigned int> assoc (dof_handler.n_dofs());
  DoFTools::get_subdomain_association (dof_handler, assoc);

  deallog << "Cell association:" << std::endl;
  for (unsigned int i=0; i<partitions.size(); ++i)
    deallog << i << ' ' << partitions(i) << std::endl;
  
  deallog << "DoF association:" << std::endl;
  for (unsigned int i=0; i<assoc.size(); ++i)
    deallog << i << ' ' << assoc[i] << std::endl;
  
}



int main ()
{
  std::ofstream logfile("metis_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test<1> ();
      test<2> ();
      test<3> ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
