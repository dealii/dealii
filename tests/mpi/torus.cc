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

/*
 * Author: Guido Kanschat, Texas A&M University, 2009
 */

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>


typedef parallel::distributed::Triangulation<2,3>::cell_iterator cell_iterator;
DeclException1(ExcMissingCell,
	       cell_iterator,
	       << "Trying to find cell " << arg1 << " but it doesn't appear to be in the list");

int main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll log;
  
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      static std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);
    }
      
  parallel::distributed::Triangulation<2,3> triangulation(MPI_COMM_WORLD,
							  typename Triangulation<2,3>::MeshSmoothing
							  (Triangulation<2,3>::smoothing_on_refinement |
							   Triangulation<2,3 >::smoothing_on_coarsening));
  GridGenerator::torus(triangulation, 1, 0.2);

  // create a set of all cells, and insert all cells into it
  std::set<parallel::distributed::Triangulation<2,3>::cell_iterator> cells;
  for (parallel::distributed::Triangulation<2,3>::cell_iterator cell= triangulation.begin(0);
       cell!=triangulation.end(0);
       ++cell)
    {
      cells.insert (cell);
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "Adding cell " << cell << std::endl;
    }
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "List contains " << cells.size() << " items" << std::endl;

  // verify that every cell is in there
  for(parallel::distributed::Triangulation<2,3>::cell_iterator cell= triangulation.begin(0);
      cell!=triangulation.end(0);
      ++cell)
    Assert (cells.find(cell)!=cells.end(),
	    ExcMissingCell(cell));

  // refine triangulation and verify that every coarse mesh cell is in there
  triangulation.refine_global(2);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "List contains " << cells.size() << " items" << std::endl;
  for(parallel::distributed::Triangulation<2,3>::cell_iterator cell= triangulation.begin(0);
      cell!=triangulation.end(0);
      ++cell)
    Assert (cells.find(cell)!=cells.end(),
	    ExcMissingCell(cell));
}
