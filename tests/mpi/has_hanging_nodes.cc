// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Test has_hanging_nodes() in sequential and parallel triangulations

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>            
#include <deal.II/base/conditional_ostream.h>  
#include <deal.II/base/index_set.h>            
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>          
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include <deal.II/distributed/tria.h>            
#include <deal.II/distributed/grid_refinement.h> 

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/lac/sparsity_tools.h>


#include <fstream>
#include <iostream>

using namespace dealii;

//#define DISTRIBUTED	

template<int dim>
void test_int(const unsigned int n_global = 0,
              const unsigned int n_local = 0)
{
   
  MPI_Comm mpi_communicator (MPI_COMM_WORLD);
  const unsigned int n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator));
   
  ConditionalOStream pcout (std::cout,
			    (Utilities::MPI::this_mpi_process(mpi_communicator)
			     == 0));


  parallel::distributed::Triangulation<dim> tria_distrib(mpi_communicator);

  Triangulation<dim> tria_sequential;


  {
    Triangulation<dim>   triangulation1;
    Triangulation<dim>   triangulation2;
    GridGenerator::hyper_cube (triangulation1, -1,0); //create a square [-1,0]^d domain
    GridGenerator::hyper_cube (triangulation2, -1,0); //create a square [-1,0]^d domain
    Point<dim> shift_vector;
    shift_vector[0] = 1.0;
    GridTools::shift(shift_vector,triangulation2);
    GridGenerator::merge_triangulations (triangulation1, triangulation2, tria_distrib);
    GridGenerator::merge_triangulations (triangulation1, triangulation2, tria_sequential);
  }
  
  tria_distrib.refine_global (n_global);
  tria_sequential.refine_global(n_global);
  
  for (unsigned int i = 0; i < n_local; i++)
    {
      // refine both the distributed and the sequential triangulation. the
      // following relies on the fact that the first cell of the "global"
      // triangulation happens to be on processor 0
      if (this_mpi_process == 0)
	{
	  for (typename Triangulation<dim>::active_cell_iterator
		 cell = tria_distrib.begin_active();
	       cell != tria_distrib.end(); ++cell)
	    if (cell->is_locally_owned())
	      {
		cell->set_refine_flag ();
		break;
	      } 
	}
  		
  				         
      tria_sequential.begin_active()->set_refine_flag ();
		             
  		
      tria_distrib.prepare_coarsening_and_refinement ();
      tria_distrib.execute_coarsening_and_refinement ();
  			
      tria_sequential.prepare_coarsening_and_refinement ();
      tria_sequential.execute_coarsening_and_refinement ();
    }
  
  deallog << n_global << " " << n_local
	  << " distrib=" << tria_distrib.has_hanging_nodes()
	  << "; shared=" << tria_sequential.has_hanging_nodes()
	  << std::endl;
}


int main (int argc, char *argv[])
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
   
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
  test_int<2>(0,0);
  test_int<2>(2,0);
  test_int<2>(3,1);
   
  test_int<3>(0,0);
  test_int<3>(2,0);
  test_int<3>(3,1);   
}
