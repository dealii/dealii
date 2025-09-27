// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/process_grid.h>

#ifdef DEAL_II_WITH_SCALAPACK

#  include <deal.II/lac/scalapack.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Internal function to determine dimension of process grid based on the total
   * number of cores, matrix dimensions and matrix block sizes
   *
   * Amesos heuristics:
   * https://github.com/trilinos/Trilinos/blob/master/packages/amesos/src/Amesos_Scalapack.cpp#L142-L166
   *
   * Elemental default grid: El::Grid::Grid(mpi_comm,...)
   * https://github.com/elemental/Elemental/blob/master/src/core/Grid.cpp#L67-L91
   */
  inline std::pair<int, int>
  compute_processor_grid_sizes(const MPI_Comm     mpi_comm,
                               const unsigned int m,
                               const unsigned int n,
                               const unsigned int block_size_m,
                               const unsigned int block_size_n)
  {
    // Few notes from the ScaLAPACK user guide:
    // It is possible to predict the best grid shape given the number of
    // processes available: Pr x Pc <= P This, however, depends on the task to
    // be done. LU , QR and QL factorizations perform better for “flat” process
    // grids (Pr < Pc ) For large N, Pc = 2*Pr is a good choice, whereas for
    // small N, one should choose small Pr Square or near square grids are more
    // optimal for Cholesky factorization. LQ and RQ factorizations take
    // advantage of “tall” grids (Pr > Pc )

    // Below we always try to create 2d processor grids:

    const int n_processes = Utilities::MPI::n_mpi_processes(mpi_comm);

    // Get the total number of cores we can occupy in a rectangular dense matrix
    // with rectangular blocks when every core owns only a single block:
    const int n_processes_heuristic = int(std::ceil((1. * m) / block_size_m)) *
                                      int(std::ceil((1. * n) / block_size_n));
    const int Np = std::min(n_processes_heuristic, n_processes);

    // Now we need to split Np into  Pr x Pc. Assume we know the shape/ratio
    // Pc =: ratio * Pr
    // therefore
    // Np = Pc * Pc / ratio
    // for quadratic matrices the ratio equals 1
    const double ratio = double(n) / m;
    int          Pc    = static_cast<int>(std::sqrt(ratio * Np));

    // one could rounds up Pc to the number which has zero remainder from the
    // division of Np while ( Np % Pc != 0 )
    //  ++Pc;
    // but this affects the grid shape dramatically, i.e. 10 cores 3x3 becomes
    // 2x5.
    // limit our estimate to be in [2, Np]
    int n_process_columns = std::min(Np, std::max(2, Pc));
    // finally, get the rows:
    int n_process_rows = Np / n_process_columns;

    Assert(n_process_columns >= 1 && n_process_rows >= 1 &&
             n_processes >= n_process_rows * n_process_columns,
           ExcMessage(
             "error in process grid: " + std::to_string(n_process_rows) + "x" +
             std::to_string(n_process_columns) + "=" +
             std::to_string(n_process_rows * n_process_columns) + " out of " +
             std::to_string(n_processes)));

    return std::make_pair(n_process_rows, n_process_columns);

    // For example,
    // 320x320 with 32x32 blocks and 16 cores:
    // Pc = 1.0 * Pr  => 4x4 grid
    // Pc = 0.5 * Pr  => 8x2 grid
    // Pc = 2.0 * Pr  => 3x5 grid
  }
} // namespace

namespace Utilities
{
  namespace MPI
  {
    ProcessGrid::ProcessGrid(
      const MPI_Comm                               mpi_comm,
      const std::pair<unsigned int, unsigned int> &grid_dimensions)
      : mpi_communicator(mpi_comm)
      , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
      , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
      , n_process_rows(grid_dimensions.first)
      , n_process_columns(grid_dimensions.second)
    {
      Assert(grid_dimensions.first > 0,
             ExcMessage("Number of process grid rows has to be positive."));
      Assert(grid_dimensions.second > 0,
             ExcMessage("Number of process grid columns has to be positive."));

      Assert(
        grid_dimensions.first * grid_dimensions.second <= n_mpi_processes,
        ExcMessage(
          "Size of process grid is larger than number of available MPI processes."));

      // processor grid order.
      const bool column_major = false;

      // Initialize Cblas context from the provided communicator
      blacs_context     = Csys2blacs_handle(mpi_communicator);
      const char *order = (column_major ? "Col" : "Row");
      // Note that blacs_context can be modified below. Thus Cblacs2sys_handle
      // may not return the same MPI communicator.
      Cblacs_gridinit(&blacs_context, order, n_process_rows, n_process_columns);

      // Blacs may modify the grid size on processes which are not used
      // in the grid. So provide copies below:
      int procrows_ = n_process_rows;
      int proccols_ = n_process_columns;
      Cblacs_gridinfo(blacs_context,
                      &procrows_,
                      &proccols_,
                      &this_process_row,
                      &this_process_column);

      // If this MPI core is not on the grid, flag it as inactive and
      // skip all jobs
      // Note that a different condition is used in FORTRAN code here
      // https://stackoverflow.com/questions/18516915/calling-blacs-with-more-processes-than-used
      if (this_process_row < 0 || this_process_column < 0)
        mpi_process_is_active = false;
      else
        mpi_process_is_active = true;

      // Create an auxiliary communicator which has root and all inactive cores.
      // Assume that inactive cores start with
      // id=n_process_rows*n_process_columns
      const unsigned int n_active_mpi_processes =
        n_process_rows * n_process_columns;
      Assert(mpi_process_is_active ||
               this_mpi_process >= n_active_mpi_processes,
             ExcInternalError());

      std::vector<int> inactive_with_root_ranks;
      inactive_with_root_ranks.push_back(0);
      for (unsigned int i = n_active_mpi_processes; i < n_mpi_processes; ++i)
        inactive_with_root_ranks.push_back(i);

      // Get the group of processes in mpi_communicator
      int       ierr = 0;
      MPI_Group all_group;
      ierr = MPI_Comm_group(mpi_communicator, &all_group);
      AssertThrowMPI(ierr);

      // Construct the group containing all ranks we need:
      MPI_Group inactive_with_root_group;
      const int n = inactive_with_root_ranks.size();
      ierr        = MPI_Group_incl(all_group,
                            n,
                            inactive_with_root_ranks.data(),
                            &inactive_with_root_group);
      AssertThrowMPI(ierr);

      // Create the communicator based on inactive_with_root_group.
      // Note that on all the active MPI processes (except for the one with
      // rank 0) the resulting MPI_Comm mpi_communicator_inactive_with_root
      // will be MPI_COMM_NULL.
      const int mpi_tag =
        Utilities::MPI::internal::Tags::process_grid_constructor;

      ierr = MPI_Comm_create_group(mpi_communicator,
                                   inactive_with_root_group,
                                   mpi_tag,
                                   &mpi_communicator_inactive_with_root);
      AssertThrowMPI(ierr);

      ierr = MPI_Group_free(&all_group);
      AssertThrowMPI(ierr);
      ierr = MPI_Group_free(&inactive_with_root_group);
      AssertThrowMPI(ierr);

      // Double check that the process with rank 0 in subgroup is active:
      if constexpr (running_in_debug_mode())
        {
          if (mpi_communicator_inactive_with_root != MPI_COMM_NULL &&
              Utilities::MPI::this_mpi_process(
                mpi_communicator_inactive_with_root) == 0)
            Assert(mpi_process_is_active, ExcInternalError());
        }
    }



    ProcessGrid::ProcessGrid(const MPI_Comm     mpi_comm,
                             const unsigned int n_rows_matrix,
                             const unsigned int n_columns_matrix,
                             const unsigned int row_block_size,
                             const unsigned int column_block_size)
      : ProcessGrid(mpi_comm,
                    compute_processor_grid_sizes(mpi_comm,
                                                 n_rows_matrix,
                                                 n_columns_matrix,
                                                 row_block_size,
                                                 column_block_size))
    {}



    ProcessGrid::ProcessGrid(const MPI_Comm     mpi_comm,
                             const unsigned int n_rows,
                             const unsigned int n_columns)
      : ProcessGrid(mpi_comm, std::make_pair(n_rows, n_columns))
    {}



    ProcessGrid::~ProcessGrid()
    {
      if (mpi_process_is_active)
        Cblacs_gridexit(blacs_context);

      if (mpi_communicator_inactive_with_root != MPI_COMM_NULL)
        Utilities::MPI::free_communicator(mpi_communicator_inactive_with_root);
    }



    template <typename NumberType>
    void
    ProcessGrid::send_to_inactive(NumberType *value, const int count) const
    {
      Assert(count > 0, ExcInternalError());
      if (mpi_communicator_inactive_with_root != MPI_COMM_NULL)
        {
          const int ierr =
            MPI_Bcast(value,
                      count,
                      Utilities::MPI::mpi_type_id_for_type<decltype(*value)>,
                      0 /*from root*/,
                      mpi_communicator_inactive_with_root);
          AssertThrowMPI(ierr);
        }
    }

  } // namespace MPI
} // namespace Utilities

// instantiations

template void
Utilities::MPI::ProcessGrid::send_to_inactive<double>(double *,
                                                      const int) const;
template void
Utilities::MPI::ProcessGrid::send_to_inactive<float>(float *, const int) const;
template void
Utilities::MPI::ProcessGrid::send_to_inactive<int>(int *, const int) const;

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK
