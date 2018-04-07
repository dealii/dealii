// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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


#ifndef dealii_matrix_free_helper_functions_h
#define dealii_matrix_free_helper_functions_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

DEAL_II_NAMESPACE_OPEN



namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * A struct that collects all information related to parallelization with
     * threads: The work is subdivided into tasks that can be done
     * independently.
     */
    struct TaskInfo
    {
      /**
       * Constructor.
       */
      TaskInfo ();

      /**
       * Clear all the data fields and resets them to zero.
       */
      void clear ();

      /**
       * Return the memory consumption of the class.
       */
      std::size_t memory_consumption () const;

      unsigned int block_size;
      unsigned int n_blocks;
      unsigned int block_size_last;
      unsigned int position_short_block;
      bool use_multithreading;
      bool use_partition_partition;
      bool use_coloring_only;

      std::vector<unsigned int> partition_color_blocks_row_index;
      std::vector<unsigned int> partition_color_blocks_data;
      unsigned int evens;
      unsigned int odds;
      unsigned int n_blocked_workers;
      unsigned int n_workers;

      std::vector<unsigned int> partition_evens;
      std::vector<unsigned int> partition_odds;
      std::vector<unsigned int> partition_n_blocked_workers;
      std::vector<unsigned int> partition_n_workers;
    };



    /**
     * A struct that collects all information related to the size of the
     * problem and MPI parallelization.
     */
    struct SizeInfo
    {
      /**
       * Constructor.
       */
      SizeInfo ();

      /**
       * Clear all data fields and resets the sizes to zero.
       */
      void clear();

      /**
       * Prints minimum, average, and maximal memory consumption over the MPI
       * processes.
       */
      template <typename StreamType>
      void print_memory_statistics (StreamType &out,
                                    std::size_t data_length) const;

      /**
       * Determines the position of cells with ghosts for distributed-memory
       * calculations.
       */
      void make_layout (const unsigned int n_active_cells_in,
                        const unsigned int vectorization_length_in,
                        std::vector<unsigned int> &boundary_cells,
                        std::vector<unsigned int> &irregular_cells);

      unsigned int n_active_cells;
      unsigned int n_macro_cells;
      unsigned int boundary_cells_start;
      unsigned int boundary_cells_end;
      unsigned int vectorization_length;

      /**
       * index sets to describe the layout of cells: locally owned cells and
       * locally active cells
       */
      IndexSet locally_owned_cells;
      IndexSet ghost_cells;

      /**
       * MPI communicator
       */
      MPI_Comm communicator;
      unsigned int my_pid;
      unsigned int n_procs;
    };

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
