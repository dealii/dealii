// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_process_grid_h
#define dealii_process_grid_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>

#ifdef DEAL_II_WITH_SCALAPACK

DEAL_II_NAMESPACE_OPEN


// Forward declaration of class ScaLAPACKMatrix for ProcessGrid
#  ifndef DOXYGEN
template <typename NumberType>
class ScaLAPACKMatrix;
#  endif

namespace Utilities
{
  namespace MPI
  {
    /**
     * A class taking care of setting up a two-dimensional processor grid.
     * For example an MPI communicator with 5 processes can be arranged into a
     * 2x2 grid with the 5-th processor being inactive:
     * @code
     *      |   0     |   1
     * -----| ------- |-----
     * 0    |   P0    |  P1
     *      |         |
     * -----| ------- |-----
     * 1    |   P2    |  P3
     * @endcode
     *
     * A shared pointer to this class is provided to ScaLAPACKMatrix matrices to
     * perform block-cyclic distribution.
     *
     * Note that this class allows to set up a process grid which has fewer
     * MPI cores than the total number of cores in the communicator.
     *
     * Currently the only place where one would use a ProcessGrid object is
     * in connection with a ScaLAPACKMatrix object.
     */
    class ProcessGrid
    {
    public:
      // Declare class ScaLAPACK as friend to provide access to private members.
      template <typename NumberType>
      friend class dealii::ScaLAPACKMatrix;

      /**
       * Constructor for a process grid with @p n_rows and @p n_columns for a given @p mpi_communicator.
       * The product of rows and columns should be less or equal to the total
       * number of cores
       * in the @p mpi_communicator.
       */
      ProcessGrid(const MPI_Comm     mpi_communicator,
                  const unsigned int n_rows,
                  const unsigned int n_columns);

      /**
       * Constructor for a process grid for a given @p mpi_communicator.
       * In this case the process grid is heuristically chosen based on the
       * dimensions and block-cyclic distribution of a target matrix provided
       * in @p n_rows_matrix, @p n_columns_matrix, @p row_block_size and @p column_block_size.
       *
       * The maximum number of MPI cores one can utilize is
       * $\min\{\frac{M}{MB}\frac{N}{NB}, Np\}$, where $M,N$ are the matrix
       * dimension and $MB,NB$ are the block sizes and $Np$ is the number of
       * processes in the @p mpi_communicator. This function then creates a 2d processor grid
       * assuming the ratio between number of process row $p$ and columns $q$ to
       * be equal the ratio between matrix dimensions $M$ and $N$.
       *
       * For example, a square matrix $640x640$ with the block size $32$
       * and the @p mpi_communicator with 11 cores will result in the $3x3$
       * process grid.
       */
      ProcessGrid(const MPI_Comm     mpi_communicator,
                  const unsigned int n_rows_matrix,
                  const unsigned int n_columns_matrix,
                  const unsigned int row_block_size,
                  const unsigned int column_block_size);

      /**
       * Destructor.
       */
      ~ProcessGrid();

      /**
       * Return the number of rows in the processes grid.
       */
      unsigned int
      get_process_grid_rows() const;

      /**
       * Return the number of columns in the processes grid.
       */
      unsigned int
      get_process_grid_columns() const;

      /**
       * Return row of this process in the process grid.
       *
       * It's negative for inactive processes.
       */
      int
      get_this_process_row() const;

      /**
       * Return column of this process in the process grid.
       *
       * It's negative for inactive processes.
       */
      int
      get_this_process_column() const;

      /**
       * Send @p count values stored consequently starting at @p value from
       * the process with rank zero to processes which
       * are not in the process grid.
       */
      template <typename NumberType>
      void
      send_to_inactive(NumberType *value, const int count = 1) const;

      /**
       * Return <code>true</code> if the process is active within the grid.
       */
      bool
      is_process_active() const;

    private:
      /**
       * A private constructor which takes grid dimensions as an
       * <code>std::pair</code>.
       */
      ProcessGrid(const MPI_Comm                               mpi_communicator,
                  const std::pair<unsigned int, unsigned int> &grid_dimensions);

      /**
       * An MPI communicator with all processes (active and inactive).
       */
      MPI_Comm mpi_communicator;

      /**
       * An MPI communicator with inactive processes and the process with rank
       * zero.
       */
      MPI_Comm mpi_communicator_inactive_with_root;

      /**
       * BLACS context. This is equivalent to MPI communicators and is
       * used by ScaLAPACK.
       */
      int blacs_context;

      /**
       * Rank of this MPI process.
       */
      const unsigned int this_mpi_process;

      /**
       * Total number of MPI processes.
       */
      const unsigned int n_mpi_processes;

      /**
       * Number of rows in the processes grid.
       */
      int n_process_rows;

      /**
       * Number of columns in the processes grid.
       */
      int n_process_columns;

      /**
       * Row of this process in the grid.
       *
       * It's negative for in-active processes.
       */
      int this_process_row;

      /**
       * Column of this process in the grid.
       *
       * It's negative for in-active processes.
       */
      int this_process_column;

      /**
       * A flag which is true for processes within the 2d process grid.
       */
      bool mpi_process_is_active;
    };

    /*--------------------- Inline functions --------------------------------*/

#  ifndef DOXYGEN

    inline unsigned int
    ProcessGrid::get_process_grid_rows() const
    {
      return n_process_rows;
    }



    inline unsigned int
    ProcessGrid::get_process_grid_columns() const
    {
      return n_process_columns;
    }



    inline int
    ProcessGrid::get_this_process_row() const
    {
      return this_process_row;
    }



    inline int
    ProcessGrid::get_this_process_column() const
    {
      return this_process_column;
    }



    inline bool
    ProcessGrid::is_process_active() const
    {
      return mpi_process_is_active;
    }


#  endif // ifndef DOXYGEN

  } // end of namespace MPI

} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK

#endif
