// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test verifies that DataOut correctly handles both DoF data and
// cell data when locally owned DoFs are intentionally made non-contiguous.

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>

#include <algorithm>

#include "../tests.h"



template <int dim>
void
force_noncontiguous_locally_owned_dofs(DoFHandler<dim> &dof_handler,
                                       const MPI_Comm   comm)
{
  const IndexSet locally_owned = dof_handler.locally_owned_dofs();
  const auto     local_n       = locally_owned.n_elements();

  const unsigned int my_rank = Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(comm);

  // How many DoFs does each rank own?
  const std::vector<types::global_dof_index> counts =
    Utilities::MPI::all_gather(comm, local_n);

  // number of ranks that have at least (k+1) locally owned DoFs
  auto active_ranks_at =
    [&](const types::global_dof_index k) -> types::global_dof_index {
    return static_cast<types::global_dof_index>(std::count_if(
      counts.begin(), counts.end(), [k](const auto c) { return c > k; }));
  };

  // Enumerate locally owned DoFs in ascending order, assign each
  // a local ordinal k = 0..local_n-1, and map that to a new global index by
  // interleaving ranks:
  std::vector<types::global_dof_index> new_numbers(local_n);
  types::global_dof_index              k = 0;
  for (const auto i : locally_owned)
    {
      // How many items were emitted in all previous k-levels?
      types::global_dof_index base = 0;
      for (types::global_dof_index kk = 0; kk < k; ++kk)
        base += active_ranks_at(kk);

      // Within this k-level, how many ranks < my_rank participate?
      const types::global_dof_index offset =
        static_cast<types::global_dof_index>(
          std::count_if(counts.begin(),
                        counts.begin() + my_rank,
                        [k](const auto c) { return c > k; }));

      const types::global_dof_index new_i = base + offset;

      new_numbers[locally_owned.index_within_set(i)] = new_i;
      ++k;
    }

  dof_handler.renumber_dofs(new_numbers);

  // locally owned IndexSet should be non-contiguous
  if (n_ranks > 1 && local_n > 1)
    AssertThrow(dof_handler.locally_owned_dofs().n_intervals() > 1,
                ExcMessage(
                  "Test setup failed: locally owned DoFs are contiguous."));
}



template <int dim>
void
run(const MPI_Comm comm)
{
  parallel::distributed::Triangulation<dim> tria(comm);
  GridGenerator::subdivided_hyper_cube(tria, 2, 0.0, 1.0);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  force_noncontiguous_locally_owned_dofs(dof_handler, comm);

  const IndexSet locally_owned = dof_handler.locally_owned_dofs();

  TrilinosWrappers::MPI::Vector solution(locally_owned, comm);
  for (const auto i : locally_owned)
    solution[i] = static_cast<double>(i);
  solution.compress(VectorOperation::insert);


  // -------------------------
  // DoF-data path test
  // -------------------------
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "u");
  data_out.build_patches();

  const auto &patches = data_out.get_patches();
  AssertThrow(!patches.empty(), ExcInternalError());

  double min_v = std::numeric_limits<double>::infinity();
  double max_v = -std::numeric_limits<double>::infinity();
  for (const auto &p : patches)
    {
      for (double v : p.data)
        {
          min_v = std::min(min_v, v);
          max_v = std::max(max_v, v);
        }

      deallog << "min_v=" << min_v << ", max_v=" << max_v << std::endl;
    }
  deallog << "patches=" << patches.size() << std::endl;



  // -------------------------
  // Cell-data path test
  // -------------------------
  const unsigned int n_global_active = tria.n_global_active_cells();

  IndexSet locally_owned_cells(n_global_active);
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      locally_owned_cells.add_index(cell->global_active_cell_index());
  locally_owned_cells.compress();

  TrilinosWrappers::MPI::Vector cell_data(locally_owned_cells, comm);

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      cell_data[cell->global_active_cell_index()] =
        static_cast<double>(cell->global_active_cell_index());

  cell_data.compress(VectorOperation::insert);

  DataOut<dim> cell_out;
  cell_out.attach_dof_handler(dof_handler);
  cell_out.add_data_vector(cell_data, "cell", DataOut<dim>::type_cell_data);
  cell_out.build_patches();

  const auto &cell_patches = cell_out.get_patches();

  double min_c = std::numeric_limits<double>::infinity();
  double max_c = -std::numeric_limits<double>::infinity();
  for (const auto &p : cell_patches)
    {
      for (double v : p.data)
        {
          min_c = std::min(min_c, v);
          max_c = std::max(max_c, v);
        }

      deallog << "min_c=" << min_c << ", max_c=" << max_c << std::endl;
    }
  deallog << "cell_patches=" << cell_patches.size() << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  run<2>(MPI_COMM_WORLD);
}
