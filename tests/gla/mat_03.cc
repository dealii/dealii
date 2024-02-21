// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// document bug in assembling a LA::MPI::SparseMatrix
// this was due to handing a wrong IndexSet to the AffineConstraints<double>

#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparsity_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"

#include "gla.h"

template <class LA, int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  AffineConstraints<double> cm;
  cm.close();

  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));
  const double R0 = 6371000. - 2890000;
  const double R1 = 6371000. - 35000.;
  GridGenerator::hyper_shell(
    triangulation, Point<dim>(), R0, R1, (dim == 3) ? 96 : 12, true);

  FE_Q<dim> temperature_fe(2);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(temperature_fe);

  const IndexSet &owned = dof_handler.locally_owned_dofs();
  const IndexSet  relevant =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  // this causes a crash in PETSc, but is ignored in Trilinos:
  // DynamicSparsityPattern sp (owned);
  DynamicSparsityPattern         sp(relevant);
  typename LA::MPI::SparseMatrix matrix;
  DoFTools::make_sparsity_pattern(dof_handler,
                                  sp,
                                  cm,
                                  false,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));
  SparsityTools::distribute_sparsity_pattern(sp,
                                             owned,
                                             MPI_COMM_WORLD,
                                             relevant);
  sp.compress();
  matrix.reinit(owned, owned, sp, MPI_COMM_WORLD);

  if (myid != 0)
    {
      types::global_dof_index                    row      = 21;
      unsigned int                               n_values = 4;
      types::global_dof_index                    cols[]   = {21, 22, 23, 39};
      typename LA::MPI::SparseMatrix::value_type vals[]   = {1, 201, 401, 101};
      matrix.add(row, n_values, cols, vals, false, true);
    }

  matrix.compress(VectorOperation::add);

  matrix.print(deallog.get_file_stream());


  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  {
    deallog.push("PETSc");
    test<LA_PETSc, 2>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos, 2>();
    deallog.pop();
  }

  // compile, don't run
  // if (myid==9999)
  //  test<LA_Dummy>();
}
