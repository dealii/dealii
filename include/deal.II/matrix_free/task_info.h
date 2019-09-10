// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
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


#ifndef dealii_matrix_free_task_info_h
#define dealii_matrix_free_task_info_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  /**
   * An interface for the worker object that runs the various operations we
   * want to perform during the matrix-free loop.
   *
   * @author Katharina Kormann, Martin Kronbichler, 2018
   */
  struct MFWorkerInterface
  {
  public:
    virtual ~MFWorkerInterface() = default;

    /// Starts the communication for the update ghost values operation
    virtual void
    vector_update_ghosts_start() = 0;

    /// Finishes the communication for the update ghost values operation
    virtual void
    vector_update_ghosts_finish() = 0;

    /// Starts the communication for the vector compress operation
    virtual void
    vector_compress_start() = 0;

    /// Finishes the communication for the vector compress operation
    virtual void
    vector_compress_finish() = 0;

    /// Zeros part of the vector accroding to a given range as stored in
    /// DoFInfo
    virtual void
    zero_dst_vector_range(const unsigned int range_index) = 0;

    virtual void
    cell_loop_pre_range(const unsigned int range_index) = 0;

    virtual void
    cell_loop_post_range(const unsigned int range_index) = 0;

    /// Runs the cell work specified by MatrixFree::loop or
    /// MatrixFree::cell_loop
    virtual void
    cell(const std::pair<unsigned int, unsigned int> &cell_range) = 0;

    /// Runs the body of the work on interior faces specified by
    /// MatrixFree::loop
    virtual void
    face(const std::pair<unsigned int, unsigned int> &face_range) = 0;

    /// Runs the body of the work on boundary faces specified by
    /// MatrixFree::loop
    virtual void
    boundary(const std::pair<unsigned int, unsigned int> &face_range) = 0;
  };



  namespace MatrixFreeFunctions
  {
    // forward declaration of internal data structure
    template <typename Number>
    struct ConstraintValues;

    /**
     * A struct that collects all information related to parallelization with
     * threads: The work is subdivided into tasks that can be done
     * independently.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2011, 2018
     */
    struct TaskInfo
    {
      // enum for choice of how to build the task graph. Odd add versions with
      // preblocking and even versions with postblocking. partition_partition
      // and partition_color are deprecated but kept for backward
      // compatibility.
      enum TasksParallelScheme
      {
        none,
        partition_partition,
        partition_color,
        color
      };

      /**
       * Constructor.
       */
      TaskInfo();

      /**
       * Clears all the data fields and resets them
       * to zero.
       */
      void
      clear();

      /**
       * Runs the matrix-free loop.
       */
      void
      loop(MFWorkerInterface &worker) const;

      /**
       * Determines the position of cells with ghosts for distributed-memory
       * calculations.
       */
      void
      collect_boundary_cells(const unsigned int n_active_cells,
                             const unsigned int n_active_and_ghost_cells,
                             const unsigned int vectorization_length,
                             std::vector<unsigned int> &boundary_cells);

      /**
       * Sets up the blocks for running the cell loop based on the options
       * controlled by the input arguments.
       *
       * @param boundary_cells A list of cells that need to exchange data prior
       * to performing computations. These will be given a certain id in the
       * partitioning.
       *
       * @param cells_close_to_boundary A list of cells that interact with
       * boundaries between different subdomains and are involved in sending
       * data to neighboring processes.
       *
       * @param dofs_per_cell Gives an expected value for the number of degrees
       * of freedom on a cell, which is used to determine the block size for
       * interleaving cell and face integrals.
       *
       * @param cell_vectorization_categories This set of categories defines
       * the cells that should be grouped together inside the lanes of a
       * vectorized array. This can be the polynomial degree in an hp-element
       * or a user-provided grouping.
       *
       * @param cell_vectorization_categories_strict Defines whether the
       * categories defined by the previous variables should be separated
       * strictly or whether it is allowed to insert lower categories into the
       * next high one(s).
       *
       * @param renumbering When leaving this function, the vector contains a
       * new numbering of the cells that aligns with the grouping stored in
       * this class.
       *
       * @param incompletely_filled_vectorization Given the vectorized layout
       * of this class, some cell batches might have components in the
       * vectorized array (SIMD lanes) that are not used and do not carray
       * valid data. This array indicates the cell batches where this occurs
       * according to the renumbering returned by this function.
       */
      void
      create_blocks_serial(
        const std::vector<unsigned int> &boundary_cells,
        const std::vector<unsigned int> &cells_close_to_boundary,
        const unsigned int               dofs_per_cell,
        const std::vector<unsigned int> &cell_vectorization_categories,
        const bool                       cell_vectorization_categories_strict,
        std::vector<unsigned int> &      renumbering,
        std::vector<unsigned char> &     incompletely_filled_vectorization);

      /**
       * First step in the block creation for the task-parallel blocking setup.
       *
       * @param boundary_cells A list of cells that need to exchange data prior
       * to performing computations. These will be given a certain id in the
       * partitioning.
       *
       * @param renumbering When leaving this function, the vector contains a
       * new numbering of the cells that aligns with the grouping stored in
       * this class (before actually creating the tasks).
       *
       * @param incompletely_filled_vectorization Given the vectorized layout
       * of this class, some cell batches might have components in the
       * vectorized array (SIMD lanes) that are not used and do not carray
       * valid data. This array indicates the cell batches where this occurs
       * according to the renumbering returned by this function.
       */
      void
      initial_setup_blocks_tasks(
        const std::vector<unsigned int> &boundary_cells,
        std::vector<unsigned int> &      renumbering,
        std::vector<unsigned char> &     incompletely_filled_vectorization);

      /**
       * This helper function determines a block size if the user decided not
       * to force a block size through MatrixFree::AdditionalData. This is
       * computed based on the number of hardware threads on the system and
       * the number of macro cells that we should work on.
       */
      void
      guess_block_size(const unsigned int dofs_per_cell);

      /**
       * This method goes through all cells that have been filled into @p
       * dof_indices and finds out which cells can be worked on independently
       * and which ones are neighboring and need to be done at different times
       * when used in parallel.
       *
       * The strategy is based on a two-level approach. The outer level is
       * subdivided into partitions similar to the type of neighbors in
       * Cuthill-McKee, and the inner level is subdivided via colors (for
       * chunks within the same color, can work independently). One task is
       * represented by a chunk of cells. The cell chunks are formed before
       * subdivision into partitions and colors.
       *
       * @param connectivity (in/out) Determines whether cells `i` and `j` are
       * conflicting, expressed by an entry in position (i,j).
       *
       * @param renumbering (in/out) At output, the element j of this variable
       * gives the original number of the cell that is reordered to place j by
       * the ordering due to the thread graph.
       *
       * @param irregular_cells (in/out) Informs the current function whether
       * some SIMD lanes in VectorizedArray would not be filled for a given
       * cell batch index.
       *
       * @param hp_bool Defines whether we are in hp mode or not
       */
      void
      make_thread_graph_partition_color(
        DynamicSparsityPattern &    connectivity,
        std::vector<unsigned int> & renumbering,
        std::vector<unsigned char> &irregular_cells,
        const bool                  hp_bool);

      /**
       * This function goes through all cells that have been filled into @p
       * dof_indices and finds out which cells can be worked on independently
       * and which ones are neighboring and need to be done at different times
       * when used in parallel.
       *
       * The strategy is based on a two-level approach. The outer level is
       * subdivided into partitions similar to the type of neighbors in
       * Cuthill-McKee, and the inner level is again subdivided into Cuthill-
       * McKee-like partitions (partitions whose level differs by more than 2
       * can be worked on independently). One task is represented by a chunk
       * of cells. The cell chunks are formed after subdivision into the two
       * levels of partitions.
       *
       * @param cell_active_fe_index The active FE index corresponding to the
       * individual indices in the list of all cell indices, in order to be
       * able to not place cells with different indices into the same cell
       * batch with vectorization.
       *
       * @param connectivity (in/out) Determines whether cells `i` and `j` are
       * conflicting, expressed by an entry in position (i,j).
       *
       * @param renumbering (in/out) At output, the element j of this variable
       * gives the original number of the cell that is reordered to place j by
       * the ordering due to the thread graph.
       *
       * @param irregular_cells (in/out) Informs the current function whether
       * some SIMD lanes in VectorizedArray would not be filled for a given
       * cell batch index.
       *
       * @param hp_bool Defines whether we are in hp mode or not
       */
      void
      make_thread_graph_partition_partition(
        const std::vector<unsigned int> &cell_active_fe_index,
        DynamicSparsityPattern &         connectivity,
        std::vector<unsigned int> &      renumbering,
        std::vector<unsigned char> &     irregular_cells,
        const bool                       hp_bool);

      /**
       * Either calls make_thread_graph_partition_color() or
       * make_thread_graph_partition_partition() accessible from the outside,
       * depending on the setting in the data structure.
       *
       * @param cell_active_fe_index The active FE index corresponding to the
       * individual indices in the list of all cell indices, in order to be
       * able to not place cells with different indices into the same cell
       * batch with vectorization.
       *
       * @param connectivity (in/out) Determines whether cells `i` and `j` are
       * conflicting, expressed by an entry in position (i,j).
       *
       * @param renumbering (in/out) At output, the element j of this variable
       * gives the original number of the cell that is reordered to place j by
       * the ordering due to the thread graph.
       *
       * @param irregular_cells (in/out) Informs the current function whether
       * some SIMD lanes in VectorizedArray would not be filled for a given
       * cell batch index.
       *
       * @param hp_bool Defines whether we are in hp mode or not
       */
      void
      make_thread_graph(const std::vector<unsigned int> &cell_active_fe_index,
                        DynamicSparsityPattern &         connectivity,
                        std::vector<unsigned int> &      renumbering,
                        std::vector<unsigned char> &     irregular_cells,
                        const bool                       hp_bool);

      /**
       * This function computes the connectivity between blocks of cells from
       * the connectivity between the individual cells.
       */
      void
      make_connectivity_cells_to_blocks(
        const std::vector<unsigned char> &irregular_cells,
        const DynamicSparsityPattern &    connectivity_cells,
        DynamicSparsityPattern &          connectivity_blocks) const;

      /**
       * Function to create coloring on the second layer within each
       * partition.
       */
      void
      make_coloring_within_partitions_pre_blocked(
        const DynamicSparsityPattern &   connectivity,
        const unsigned int               partition,
        const std::vector<unsigned int> &cell_partition,
        const std::vector<unsigned int> &partition_list,
        const std::vector<unsigned int> &partition_size,
        std::vector<unsigned int> &      partition_color_list);

      /**
       * Function to create partitioning on the second layer within each
       * partition.
       */
      void
      make_partitioning_within_partitions_post_blocked(
        const DynamicSparsityPattern &   connectivity,
        const std::vector<unsigned int> &cell_active_fe_index,
        const unsigned int               partition,
        const unsigned int               cluster_size,
        const bool                       hp_bool,
        const std::vector<unsigned int> &cell_partition,
        const std::vector<unsigned int> &partition_list,
        const std::vector<unsigned int> &partition_size,
        std::vector<unsigned int> &      partition_partition_list,
        std::vector<unsigned char> &     irregular_cells);

      /**
       * This function creates partitions according to the provided connectivity
       * graph.
       *
       * @param connectivity Connectivity between (blocks of cells)
       *
       * @param cluster_size The number of cells in each partition should be a
       * multiple of cluster_size (for blocking later on)
       *
       * @param cell_partition Saves of each (block of cells) to which
       * partition the block belongs
       *
       * @param partition_list partition_list[j] gives the old number of the
       * block that should be renumbered to j due to the partitioning
       *
       * @param partition_size Vector pointing to start of each partition (on
       * output)
       *
       * @param partition number of partitions created
       */
      void
      make_partitioning(const DynamicSparsityPattern &connectivity,
                        const unsigned int            cluster_size,
                        std::vector<unsigned int> &   cell_partition,
                        std::vector<unsigned int> &   partition_list,
                        std::vector<unsigned int> &   partition_size,
                        unsigned int &                partition) const;

      /**
       * Update fields of task info for task graph set up in
       * make_thread_graph.
       */
      void
      update_task_info(const unsigned int partition);

      /**
       * Creates a task graph from a connectivity structure.
       */
      void
      create_flow_graph();

      /**
       * Returns the memory consumption of the class.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Prints minimum, average, and maximal memory consumption over the MPI
       * processes.
       */
      template <typename StreamType>
      void
      print_memory_statistics(StreamType &out, std::size_t data_length) const;

      /**
       * Number of physical cells in the mesh, not cell batches after
       * vectorization
       */
      unsigned int n_active_cells;

      /**
       * Number of physical ghost cells in the mesh which are subject to
       * special treatment and should not be included in loops
       */
      unsigned int n_ghost_cells;

      /**
       * Number of lanes in the SIMD array that are used for vectorization
       */
      unsigned int vectorization_length;

      /**
       * Block size information for multithreading
       */
      unsigned int block_size;

      /**
       * Number of blocks for multithreading
       */
      unsigned int n_blocks;

      /**
       * Parallel scheme applied by multithreading
       */
      TasksParallelScheme scheme;

      /**
       * The blocks are organized by a vector-of-vector concept, and this data
       * field @p partition_row_index stores the distance from one 'vector' to
       * the next within the linear storage of all data to the two-level
       * partitioning.
       */
      std::vector<unsigned int> partition_row_index;

      /**
       * This is a linear storage of all partitions, building a range of
       * indices of the form cell_partition_data[idx] to
       * cell_partition_data[idx+1] within the integer list of all cells in
       * MatrixFree, subdivided into chunks by @p partition_row_index.
       */
      std::vector<unsigned int> cell_partition_data;

      /**
       * This is a linear storage of all partitions of inner faces, building a
       * range of indices of the form face_partition_data[idx] to
       * face_partition_data[idx+1] within the integer list of all interior
       * faces in MatrixFree, subdivided into chunks by @p
       * partition_row_index.
       */
      std::vector<unsigned int> face_partition_data;

      /**
       * This is a linear storage of all partitions of boundary faces,
       * building a range of indices of the form boundary_partition_data[idx]
       * to boundary_partition_data[idx+1] within the integer list of all
       * boundary faces in MatrixFree, subdivided into chunks by @p
       * partition_row_index.
       */
      std::vector<unsigned int> boundary_partition_data;

      /**
       * This is a linear storage of all partitions of interior faces on
       * boundaries to other processors that are not locally used, building a
       * range of indices of the form ghost_face_partition_data[idx] to
       * ghost_face_partition_data[idx+1] within the integer list of all such
       * faces in MatrixFree, subdivided into chunks by @p
       * partition_row_index.
       */
      std::vector<unsigned int> ghost_face_partition_data;

      /**
       * This is a linear storage of all partitions of faces for multigrid
       * levels that have a coarser neighbor and are only included in certain
       * residual computations but not in smoothing, building a range of
       * indices of the form refinement_edge_face_partition_data[idx] to
       * refinement_edge_face_partition_data[idx+1] within the integer list of
       * all such faces in MatrixFree, subdivided into chunks by @p
       * partition_row_index.
       */
      std::vector<unsigned int> refinement_edge_face_partition_data;

      /**
       * Thread information (which chunk to start 'even' partitions from) to
       * be handed to the dynamic task scheduler
       */
      std::vector<unsigned int> partition_evens;

      /**
       * Thread information (which chunk to start 'odd' partitions from) to be
       * handed to the dynamic task scheduler
       */
      std::vector<unsigned int> partition_odds;

      /**
       * Thread information regarding the dependencies for partitions handed
       * to the dynamic task scheduler
       */
      std::vector<unsigned int> partition_n_blocked_workers;

      /**
       * Thread information regarding the dependencies for partitions handed
       * to the dynamic task scheduler
       */
      std::vector<unsigned int> partition_n_workers;

      /**
       * Number of even partitions accumulated over the field @p
       * partitions_even
       */
      unsigned int evens;

      /**
       * Number of odd partitions accumulated over the field @p
       * partitions_odd
       */
      unsigned int odds;

      /**
       * Number of blocked workers accumulated over the field @p
       * partition_n_blocked_workers
       */
      unsigned int n_blocked_workers;

      /**
       * Number of workers accumulated over the field @p partition_n_workers
       */
      unsigned int n_workers;

      /**
       * Stores whether a particular task is at an MPI boundary and needs data
       * exchange
       */
      std::vector<unsigned char> task_at_mpi_boundary;

      /**
       * MPI communicator
       */
      MPI_Comm communicator;

      /**
       * Rank of MPI process
       */
      unsigned int my_pid;

      /**
       * Number of MPI rank for the current communicator
       */
      unsigned int n_procs;
    };

    /**
     * Typedef to deprecated name.
     */
    using SizeInfo DEAL_II_DEPRECATED = TaskInfo;

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
