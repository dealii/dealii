// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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



// LA::MPI::SparseMatrix, make sure it has the same
// entries in Trilinos and PETSc (also when using
// a TrilinosWrappers::SparsityPattern instead)

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparsity_tools.h>

#include "gla.h"

template <class LA, int dim>
void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  ConstraintMatrix cm;
  cm.close();

  parallel::distributed::Triangulation<dim>
  triangulation (MPI_COMM_WORLD,
                 typename Triangulation<dim>::MeshSmoothing
                 (Triangulation<dim>::smoothing_on_refinement |
                  Triangulation<dim>::smoothing_on_coarsening));
  const double R0      = 6371000.-2890000;
  const double R1      = 6371000.-  35000.;
  GridGenerator::hyper_shell (triangulation,
                              Point<dim>(),
                              R0,
                              R1,
                              (dim==3) ? 96 : 12,
                              true);

  FE_Q<dim>                                 temperature_fe(1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs (temperature_fe);

  IndexSet owned = dof_handler.locally_owned_dofs();
  IndexSet relevant;
  DoFTools::extract_locally_relevant_dofs (dof_handler, relevant);

  CompressedSimpleSparsityPattern sp (relevant);
  typename LA::MPI::SparseMatrix matrix;
  DoFTools::make_sparsity_pattern (dof_handler, sp,
                                   cm, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  SparsityTools::distribute_sparsity_pattern(sp,
                                             dof_handler.n_locally_owned_dofs_per_processor(),
                                             MPI_COMM_WORLD, relevant);
  sp.compress();
  matrix.reinit (owned, owned, sp, MPI_COMM_WORLD);


  matrix.print(deallog.get_file_stream());


  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



template <int dim>
void test_trilinos_alternative ()
{
  typedef LA_Trilinos LA;

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  ConstraintMatrix cm;
  cm.close();

  parallel::distributed::Triangulation<dim>
  triangulation (MPI_COMM_WORLD,
                 typename Triangulation<dim>::MeshSmoothing
                 (Triangulation<dim>::smoothing_on_refinement |
                  Triangulation<dim>::smoothing_on_coarsening));
  const double R0      = 6371000.-2890000;
  const double R1      = 6371000.-  35000.;
  GridGenerator::hyper_shell (triangulation,
                              Point<dim>(),
                              R0,
                              R1,
                              (dim==3) ? 96 : 12,
                              true);

  FE_Q<dim>                                 temperature_fe(1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs (temperature_fe);

  IndexSet owned = dof_handler.locally_owned_dofs();
  IndexSet relevant;
  DoFTools::extract_locally_relevant_dofs (dof_handler, relevant);

  TrilinosWrappers::SparsityPattern sp (owned, MPI_COMM_WORLD);
  typename LA::MPI::SparseMatrix matrix;
  DoFTools::make_sparsity_pattern (dof_handler, sp,
                                   cm, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();
  matrix.reinit (sp);

  matrix.print(deallog.get_file_stream());

  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}

int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  {
    deallog.push("PETSc");
    test<LA_PETSc,2>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos,2>();
    deallog.pop();
    deallog.push("Trilinos_alt");
    test_trilinos_alternative<2>();
    deallog.pop();
  }

  // compile, don't run
  //if (myid==9999)
  //  test<LA_Dummy>();


}
