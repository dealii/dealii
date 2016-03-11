// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
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


#ifndef dealii__matrix_free_helper_functions_h
#define dealii__matrix_free_helper_functions_h


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
       * Clears all the data fields and resets them to zero.
       */
      void clear ();

      /**
       * Returns the memory consumption of the class.
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
       * Clears all data fields and resets the sizes to zero.
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

    /**
     * Data type to identify cell type.
     */
    enum CellType {cartesian=0, affine=1, general=2, undefined=3};


    /**
     * A class that is used to compare floating point arrays (e.g.
     * std::vectors, Tensor<1,dim>, etc.). The idea of this class is to
     * consider two arrays as equal if they are the same within a given
     * tolerance. We use this comparator class within an std::map<> of the
     * given arrays. Note that this comparison operator does not satisfy all
     * the mathematical properties one usually wants to have (consider e.g.
     * the numbers a=0, b=0.1, c=0.2 with tolerance 0.15; the operator gives
     * a<c, but neither of a<b? or b<c? is satisfied). This is not a problem
     * in the use cases for this class, but be careful when using it in other
     * contexts.
     */
    template<typename Number>
    struct FPArrayComparator
    {
      FPArrayComparator (const Number scaling);

      bool operator() (const std::vector<Number> &v1,
                       const std::vector<Number> &v2) const;

      template <int dim>
      bool operator ()(const Tensor<1,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t1,
                       const Tensor<1,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t2) const;

      template <int dim>
      bool operator ()(const Tensor<2,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t1,
                       const Tensor<2,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t2) const;

      Number tolerance;
    };

    // Note: Implementation in matrix_free.templates.h

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
