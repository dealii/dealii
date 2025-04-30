// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check n_nonzero_elements() for an empty matrix. We generate a few dummy
// non-zero entries.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<2>       fe(1);
  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Table<2, DoFTools::Coupling> coupling(1, 1);
  coupling.fill(DoFTools::none);

  const IndexSet &system_partitioning = dof_handler.locally_owned_dofs();
  const IndexSet  system_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dof_handler);


  // create an empty sparsity pattern
  LinearAlgebra::TpetraWrappers::SparsityPattern sparsity;
  sparsity.reinit(system_partitioning,
                  system_partitioning,
                  system_relevant_set,
                  MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  coupling,
                                  sparsity,
                                  AffineConstraints<double>(),
                                  false,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));
  sparsity.compress();

  // attach a sparse matrix to it
  LinearAlgebra::TpetraWrappers::SparseMatrix<double> A;
  A.reinit(sparsity);

  // See how many nonzero elements it reports. Note that we have a
  // nonlocal Trilinos graph for the off-processor entries and we need
  // to enter a few dummy elements:
  deallog << A.n_nonzero_elements() << std::endl;
  A.print(deallog.get_file_stream());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
