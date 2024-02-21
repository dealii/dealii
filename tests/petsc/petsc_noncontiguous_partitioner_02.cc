// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PETScWrappers::CommunicationPattern for repartitioning of a
// degrees of freedom (Morton order -> layer partitioning)
// Copy-pasted from tests/base/mpi_noncontiguous_partitioner_02.cc

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/petsc_communication_pattern.h>

#include "../tests.h"


template <int dim>
void
test(const MPI_Comm comm, const bool do_revert, const unsigned int dir)
{
  const unsigned int degree           = 2;
  const unsigned int n_refinements    = 2;
  const unsigned int n_cells_1D       = Utilities::pow(2, n_refinements);
  const unsigned int n_points_1D      = (degree + 1) * n_cells_1D;
  const double       delta            = 1.0 / n_cells_1D;
  const unsigned int n_points_cell_1D = degree + 1;
  const unsigned int n_points_cell    = Utilities::pow(degree + 1, dim);
  const unsigned int n_points_face    = Utilities::pow(n_points_1D, dim - 1);

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(comm);
  const unsigned int my_rank =
    do_revert ? (n_procs - 1 - Utilities::MPI::this_mpi_process(comm)) :
                Utilities::MPI::this_mpi_process(comm);

  parallel::distributed::Triangulation<dim> tria(comm);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  FE_DGQ<dim> fe(degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const unsigned n_dofs = dof_handler.n_dofs();


  std::vector<types::global_dof_index> indices_has, indices_want;

  auto norm_point_to_lex = [&](const Point<dim> &c) {
    // convert normalized point [0, 1] to lex
    if (dim == 2)
      return std::floor(c[0]) + n_cells_1D * std::floor(c[1]);
    else
      return std::floor(c[0]) + n_cells_1D * std::floor(c[1]) +
             n_cells_1D * n_cells_1D * std::floor(c[2]);
  };

  // ... has (symm)
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_active() && cell->is_locally_owned())
      {
        auto c = cell->center();
        for (unsigned int i = 0; i < dim; ++i)
          c[i] = c[i] / delta;

        const auto lid = static_cast<unsigned int>(norm_point_to_lex(c));

        for (unsigned int i = lid * n_points_cell;
             i < (lid + 1) * n_points_cell;
             i++)
          indices_has.push_back(i);
      }


  const unsigned int div = n_points_1D / n_procs;
  const unsigned int rem = n_points_1D % n_procs;

  const unsigned int start = div * (my_rank + 0) + std::min(my_rank + 0, rem);
  const unsigned int end   = div * (my_rank + 1) + std::min(my_rank + 1, rem);

  if (dim == 2 && dir == 0)
    {
      for (unsigned int j = 0, c = start * n_points_face; j < end; ++j)
        for (unsigned int i = start; i < n_points_1D; ++i)
          indices_want.push_back(c++);
    }
  else if (dim == 2 && dir == 1)
    {
      for (unsigned int j = 0; j < n_points_1D; ++j)
        for (unsigned int i = start; i < end; ++i)
          indices_want.push_back(j * n_points_face + i);
    }
  else
    Assert(false, StandardExceptions::ExcNotImplemented());

  if (do_revert)
    std::reverse(indices_want.begin(), indices_want.end());

  PETScWrappers::CommunicationPattern vector;
  vector.reinit(indices_has, indices_want, comm);

  AlignedVector<double> src(indices_has.size());
  for (unsigned int i = 0; i < indices_has.size(); ++i)
    src[i] = indices_has[i];


  AlignedVector<double> dst(indices_want.size());

  vector.export_to_ghosted_array(ArrayView<const double>(src.data(),
                                                         src.size()),
                                 ArrayView<double>(dst.data(), dst.size()));

  for (size_t i = 0; i < src.size(); ++i)
    deallog << static_cast<int>(src[i]) << ' ';
  deallog << std::endl;
  for (size_t i = 0; i < dst.size(); ++i)
    deallog << static_cast<int>(dst[i]) << ' ';
  deallog << std::endl << std::endl;


  for (size_t i = 0; i < dst.size(); ++i)
    AssertDimension(dst[i], indices_want[i]);
}

template <int dim>
void
test_dim(const MPI_Comm comm, const bool do_revert)
{
  for (int dir = 0; dir < dim; ++dir)
    test<dim>(comm, do_revert, dir);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  test_dim<2>(comm, /*do_revert=*/false);
  test_dim<2>(comm, /*do_revert=*/true);
}
