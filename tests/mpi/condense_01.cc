// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// test ConstraintMatrix::condense(in, out)

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_tools.h>

#include <fstream>


template<int dim>
void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube (triangulation, -0.5, 0.5);

  triangulation.refine_global(1);

  /*
  for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active (); cell != triangulation.end (); ++cell)
    if (cell->is_locally_owned ())
    {
      if(cell->center().square() < 0.3*0.3)
        cell->set_refine_flag ();
    }
  */
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active ();
  if (cell->is_locally_owned ())
    cell->set_refine_flag ();
  
  triangulation.execute_coarsening_and_refinement ();
  if (myid == 0)
    deallog << "#cells = " << triangulation.n_global_active_cells() << std::endl;

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close ();

  IndexSet locally_owned_dofs, locally_relevant_dofs;
  DoFTools::extract_locally_owned_dofs (dof_handler, locally_owned_dofs);
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

  PETScWrappers::MPI::Vector vec (locally_owned_dofs, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector vec_ghosted (locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  vec = 1;
  vec_ghosted = vec;
  vec = -1;
  
  constraints.condense (vec_ghosted, vec);

  vec.print(deallog.get_file_stream());

  if (myid==0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  test<2>();
}
