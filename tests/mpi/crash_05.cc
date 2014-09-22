// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// test hits an exception in hanging_node constraints when using RT elements in parallel
// reported by francesco.cattoglio on the mailing list on 2013/11/22
/*
905: --------------------------------------------------------
905: An error occurred in line <1898> of file </scratch/deal-trunk/deal.II/include/deal.II/lac/constraint_matrix.h> in function
905:     void dealii::ConstraintMatrix::add_line(dealii::ConstraintMatrix::size_type)
905: The violated condition was: 
905:     line != numbers::invalid_size_type
905: The name and call sequence of the exception was:
905:     ExcInternalError()
905: Additional Information: 
905: (none)
905: --------------------------------------------------------
*/

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <cstdlib>
#include <numeric>


template<int dim>
void test()
{
    parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
    GridGenerator::hyper_ball(triangulation);
    
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
	typename Triangulation<dim>::active_cell_iterator
	  cell = triangulation.begin_active();
	//cell->set_refine_flag ();
	(++cell)->set_refine_flag ();
      }

    triangulation.execute_coarsening_and_refinement();
    deallog << "n_cells: " << triangulation.n_global_active_cells() << std::endl;
    FE_RaviartThomas<dim> fe(0); // crashing
    //FESystem<dim> fe(FE_Q<dim>(2), 1); // working
    DoFHandler<dim> dof_handler(triangulation);

    dof_handler.distribute_dofs(fe);
	
    ConstraintMatrix constraints;
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    IndexSet relevant_set;
    DoFTools::extract_locally_relevant_dofs (dof_handler, relevant_set);
    deallog << "relevant set:" << std::endl;
    relevant_set.print(deallog.get_file_stream());
    deallog << "constraints:" << std::endl;
    constraints.print(deallog.get_file_stream());
    
    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{  
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  test<2>();
  
}
