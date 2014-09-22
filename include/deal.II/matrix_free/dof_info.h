// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


#ifndef __deal2__matrix_free_dof_info_h
#define __deal2__matrix_free_dof_info_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/matrix_free/helper_functions.h>

#include <deal.II/base/std_cxx11/array.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * The class that stores the indices of the degrees of freedom for all the
     * cells. Essentially, this is a smart number cache in the style of a
     * DoFHandler that also embeds the description of constraints directly on
     * the cell level without the need to refer to the external
     * ConstraintMatrix.
     *
     * This class only stores index relations. The weights for hanging node
     * constraints are stored in a different field. This is because a
     * different field allows for the same compressed weight data on different
     * DoFHandlers for vector-valued problems. There, the indices might be
     * constrained differently on different components (e.g. Dirichlet
     * conditions only on selected components), whereas the weights from
     * hanging nodes are the same and need to be stored only once. The
     * combination will be handled in the MatrixFree class.
     *
     * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
     */
    struct DoFInfo
    {
      /**
       * Default empty constructor.
       */
      DoFInfo ();

      /**
       * Copy constructor.
       */
      DoFInfo (const DoFInfo &dof_info);

      /**
       * Clears all data fields in this class.
       */
      void clear ();


      /**
       * Returns a pointer to the first index in the DoF row @p row.
       */
      const unsigned int *begin_indices (const unsigned int row) const;

      /**
       * Returns a pointer to the one past the last DoF index in the row @p
       * row.
       */
      const unsigned int *end_indices (const unsigned int row) const;

      /**
       * Returns the number of entries in the indices field for the given row.
       */
      unsigned int row_length_indices (const unsigned int row) const;

      /**
       * Returns a pointer to the first constraint indicator in the row @p
       * row.
       */
      const std::pair<unsigned short,unsigned short> *
      begin_indicators (const unsigned int row) const;

      /**
       * Returns a pointer to the one past the last constraint indicator in
       * the row @p row.
       */
      const std::pair<unsigned short,unsigned short> *
      end_indicators (const unsigned int row) const;

      /**
       * Returns the number of entries in the constraint indicator field for
       * the given row.
       */
      unsigned int row_length_indicators (const unsigned int row) const;

      /**
       * Returns a pointer to the first index in the DoF row @p row for plain
       * indices (i.e., the entries where constraints are not embedded).
       */
      const unsigned int *begin_indices_plain (const unsigned int row) const;

      /**
       * Returns a pointer to the one past the last DoF index in the row @p
       * row (i.e., the entries where constraints are not embedded).
       */
      const unsigned int *end_indices_plain (const unsigned int row) const;

      /**
       * Returns the FE index for a given finite element degree. If not in hp
       * mode, this function always returns index 0. If an index is not found
       * in hp mode, it returns max_fe_degree, i.e., one index past the last
       * valid one.
       */
      unsigned int fe_index_from_degree (const unsigned int fe_degree) const;


      /**
       * Returns the FE index for a given finite element degree. If not in hp
       * mode or if the index is not found, this function always returns index
       * 0. Hence, this function does not check whether the given degree is
       * actually present.
       */
      unsigned int
      fe_index_from_dofs_per_cell (const unsigned int dofs_per_cell) const;

      /**
       * This internal method takes the local indices on a cell and fills them
       * into this class. It resolves the constraints and distributes the
       * results. Ghost indices, i.e., indices that are located on another
       * processor, get a temporary number by this function, and will later be
       * assigned the correct index after all the ghost indices have been
       * collected by the call to @p assign_ghosts.
       */
      void read_dof_indices (const std::vector<types::global_dof_index> &local_indices,
                             const std::vector<unsigned int> &lexicographic_inv,
                             const ConstraintMatrix          &constraints,
                             const unsigned int               cell_number,
                             ConstraintValues<double> &constraint_values,
                             bool                            &cell_at_boundary);

      /**
       * This method assigns the correct indices to ghost indices from the
       * temporary numbering employed by the @p read_dof_indices function. The
       * numbers are localized with respect to the MPI process, and ghosts
       * start at the end of the locally owned range. This way, we get direct
       * access to all vector entries.
       */
      void assign_ghosts(const std::vector<unsigned int> &boundary_cells);

      /**
       * Reorganizes cells for serial (non-thread-parallelized) such that
       * boundary cells are places in the middle. This way, computations and
       * communication can be overlapped. Should only be called by one DoFInfo
       * object when used on a system of several DoFHandlers.
       */
      void compute_renumber_serial (const std::vector<unsigned int> &boundary_cells,
                                    const SizeInfo                  &size_info,
                                    std::vector<unsigned int>       &renumbering);

      /**
       * Reorganizes cells in the hp case without parallelism such that all
       * cells with the same FE index are placed consecutively. Should only be
       * called by one DoFInfo object when used on a system of several
       * DoFHandlers.
       */
      void compute_renumber_hp_serial (SizeInfo                  &size_info,
                                       std::vector<unsigned int> &renumbering,
                                       std::vector<unsigned int> &irregular_cells);

      /**
       * Computes the initial renumbering of cells such that all cells with
       * ghosts are put first. This is the first step before building the
       * thread graph and used to overlap computations and communication.
       */
      void compute_renumber_parallel (const std::vector<unsigned int> &boundary_cells,
                                      SizeInfo                        &size_info,
                                      std::vector<unsigned int>       &renumbering);

      /**
       * This method reorders the way cells are gone through based on a given
       * renumbering of the cells. It also takes @p vectorization_length cells
       * together and interprets them as one cell only, as is needed for
       * vectorization.
       */
      void reorder_cells (const SizeInfo                  &size_info,
                          const std::vector<unsigned int> &renumbering,
                          const std::vector<unsigned int> &constraint_pool_row_index,
                          const std::vector<unsigned int> &irregular_cells,
                          const unsigned int               vectorization_length);

      /**
       * This helper function determines a block size if the user decided not
       * to force a block size through MatrixFree::AdditionalData. This is
       * computed based on the number of hardware threads on the system
       *  and the number of macro cells that we should work on.
       */
      void guess_block_size (const SizeInfo &size_info,
                             TaskInfo       &task_info);

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
       */
      void
      make_thread_graph_partition_color (SizeInfo                  &size_info,
                                         TaskInfo                  &task_info,
                                         std::vector<unsigned int> &renumbering,
                                         std::vector<unsigned int> &irregular_cells,
                                         const bool                 hp_bool);

      /**
       * This function goes through all cells that have been filled into @p
       * dof_indices and finds out which cells can be worked on independently
       * and which ones are neighboring and need to be done at different times
       * when used in parallel.
       *
       * The strategy is based on a two-level approach. The outer level is
       * subdivided into partitions similar to the type of neighbors in
       * Cuthill-McKee, and the inner level is again subdivided into
       * Cuthill-McKee-like partitions (partitions whose level differs by more
       * than 2 can be worked on independently). One task is represented by a
       * chunk of cells. The cell chunks are formed after subdivision into the
       * two levels of partitions.
       */
      void
      make_thread_graph_partition_partition (SizeInfo                  &size_info,
                                             TaskInfo                  &task_info,
                                             std::vector<unsigned int> &renumbering,
                                             std::vector<unsigned int> &irregular_cells,
                                             const bool                 hp_bool);

      /**
       * This function computes the connectivity of the currently stored
       * indices and fills the structure into a sparsity pattern. The
       * parameter block_size can be used to specify whether several cells
       * should be treated as one.
       */
      void
      make_connectivity_graph (const SizeInfo                  &size_info,
                               const TaskInfo                  &task_info,
                               const std::vector<unsigned int> &renumbering,
                               const std::vector<unsigned int> &irregular_cells,
                               const bool                       do_blocking,
                               CompressedSimpleSparsityPattern &connectivity) const;

      /**
       * Renumbers the degrees of freedom to give good access for this class.
       */
      void renumber_dofs (std::vector<types::global_dof_index> &renumbering);

      /**
       * Returns the memory consumption in bytes of this class.
       */
      std::size_t memory_consumption() const;

      /**
       * Prints a detailed summary of memory consumption in the different
       * structures of this class to the given output stream.
       */
      template <typename STREAM>
      void print_memory_consumption(STREAM         &out,
                                    const SizeInfo &size_info) const;

      /**
       * Prints a representation of the indices in the class to the given
       * output stream.
       */
      template <typename Number>
      void print (const std::vector<Number>       &constraint_pool_data,
                  const std::vector<unsigned int> &constraint_pool_row_index,
                  std::ostream                    &out) const;

      /**
       * Stores the rowstart indices of the compressed row storage in the @p
       * dof_indices and @p constraint_indicator fields. These two fields are
       * always accessed together, so it is simpler to keep just one variable
       * for them. This also obviates keeping two rowstart vectors in synch.
       *
       * In addition, the third field stores whether a particular cell has a
       * certain structure in the indices, like indices for vector-valued
       * problems or for cells where not all vector components are filled.
       */
      std::vector<std_cxx11::array<unsigned int, 3> > row_starts;

      /**
       * Stores the indices of the degrees of freedom for each cell. These
       * indices are computed in MPI-local index space, i.e., each processor
       * stores the locally owned indices as numbers between <tt>0</tt> and
       * <tt>n_locally_owned_dofs-1</tt> and ghost indices in the range
       * <tt>n_locally_owned_dofs</tt> to
       * <tt>n_locally_owned_dofs+n_ghost_dofs</tt>. The translation between
       * this MPI-local index space and the global numbering of degrees of
       * freedom is stored in the @p vector_partitioner data structure.

       * This array also includes the indirect contributions from constraints,
       * which are described by the @p constraint_indicator field. Because of
       * variable lengths of rows, this would be a vector of a
       * vector. However, we use one contiguous memory region and store the
       * rowstart in the variable @p row_starts.
       */
      std::vector<unsigned int> dof_indices;

      /**
       * This variable describes the position of constraints in terms of the
       * local numbering of degrees of freedom on a cell. The first number
       * stores the distance from one constrained degree of freedom to the
       * next. This allows to identify the position of constrained DoFs as we
       * loop through the local degrees of freedom of the cell when reading
       * from or writing to a vector. The second number stores the index of
       * the constraint weights, stored in the variable constraint_pool_data.
       */
      std::vector<std::pair<unsigned short,unsigned short> > constraint_indicator;

      /**
       * This stores the parallel partitioning that can be used to set up
       * vectors. The partitioner includes the description of the local range
       * in the vector, and also includes how the ghosts look like. This
       * enables initialization of vectors based on the DoFInfo field.
       */
      std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> vector_partitioner;

      /**
       * This stores a (sorted) list of all locally owned degrees of freedom
       * that are constrained.
       */
      std::vector<unsigned int> constrained_dofs;

      /**
       * Stores the rowstart indices of the compressed row storage in the @p
       * plain_dof_indices fields.
       */
      std::vector<unsigned int> row_starts_plain_indices;

      /**
       * Stores the indices of the degrees of freedom for each cell. This
       * array does not include the indirect contributions from constraints,
       * which are included in @p dof_indices. Because of variable lengths of
       * rows, this would be a vector of a vector. However, we use one
       * contiguous memory region and store the rowstart in the variable @p
       * row_starts_plain_indices.
       */
      std::vector<unsigned int> plain_dof_indices;

      /**
       * Stores the dimension of the underlying DoFHandler. Since the indices
       * are not templated, this is the variable that makes the dimension
       * accessible in the (rare) cases it is needed inside this class.
       */
      unsigned int dimension;

      /**
       * Stores the number of components in the DoFHandler where the indices
       * have been read from.
       */
      unsigned int n_components;

      /**
       * Stores the number of degrees of freedom per cell.
       */
      std::vector<unsigned int> dofs_per_cell;

      /**
       * Stores the number of degrees of freedom per face.
       */
      std::vector<unsigned int> dofs_per_face;

      /**
       * Informs on whether plain indices are cached.
       */
      bool store_plain_indices;

      /**
       * Stores the index of the active finite element in the hp case.
       */
      std::vector<unsigned int> cell_active_fe_index;

      /**
       * Stores the maximum degree of different finite elements for the hp
       * case.
       */
      unsigned int max_fe_index;

      /**
       * This variable stores the dofs per cell and the finite element degree
       * associated for all fe indices in the underlying element for easier
       * access to data in the hp case.
       */
      std::vector<std::pair<unsigned int,unsigned int> > fe_index_conversion;

      /**
       * Temporarily stores the numbers of ghosts during setup. Cleared when
       * calling @p assign_ghosts. Then, all information is collected by the
       * partitioner.
       */
      std::vector<types::global_dof_index> ghost_dofs;
    };


    /*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

    inline
    const unsigned int *
    DoFInfo::begin_indices (const unsigned int row) const
    {
      AssertIndexRange (row, row_starts.size()-1);
      const unsigned int index = row_starts[row][0];
      AssertIndexRange(index, dof_indices.size()+1);
      return dof_indices.empty() ? 0 : &dof_indices[0] + index;
    }



    inline
    const unsigned int *
    DoFInfo::end_indices (const unsigned int row) const
    {
      AssertIndexRange (row, row_starts.size()-1);
      const unsigned int index = row_starts[row+1][0];
      AssertIndexRange(index, dof_indices.size()+1);
      return dof_indices.empty() ? 0 : &dof_indices[0] + index;
    }



    inline
    unsigned int
    DoFInfo::row_length_indices (const unsigned int row) const
    {
      AssertIndexRange (row, row_starts.size()-1);
      return (row_starts[row+1][0] - row_starts[row][0]);
    }



    inline
    const std::pair<unsigned short,unsigned short> *
    DoFInfo::begin_indicators (const unsigned int row) const
    {
      AssertIndexRange (row, row_starts.size()-1);
      const unsigned int index = row_starts[row][1];
      AssertIndexRange (index, constraint_indicator.size()+1);
      return constraint_indicator.empty() ? 0 : &constraint_indicator[0] + index;
    }



    inline
    const std::pair<unsigned short,unsigned short> *
    DoFInfo::end_indicators (const unsigned int row) const
    {
      AssertIndexRange (row, row_starts.size()-1);
      const unsigned int index = row_starts[row+1][1];
      AssertIndexRange (index, constraint_indicator.size()+1);
      return constraint_indicator.empty() ? 0 : &constraint_indicator[0] + index;
    }



    inline
    unsigned int
    DoFInfo::row_length_indicators (const unsigned int row) const
    {
      AssertIndexRange (row, row_starts.size()-1);
      return (row_starts[row+1][1] - row_starts[row][1]);
    }



    inline
    const unsigned int *
    DoFInfo::begin_indices_plain (const unsigned int row) const
    {
      // if we have no constraints, should take the data from dof_indices
      if (row_length_indicators(row) == 0)
        {
          Assert (row_starts_plain_indices[row]==numbers::invalid_unsigned_int,
                  ExcInternalError());
          return begin_indices(row);
        }
      else
        {
          AssertDimension (row_starts.size(), row_starts_plain_indices.size());
          const unsigned int index = row_starts_plain_indices[row];
          AssertIndexRange(index, plain_dof_indices.size()+1);
          return plain_dof_indices.empty() ? 0 : &plain_dof_indices[0] + index;
        }
    }



    inline
    const unsigned int *
    DoFInfo::end_indices_plain (const unsigned int row) const
    {
      return begin_indices_plain(row) +
             dofs_per_cell[(cell_active_fe_index.size()==0)?
                           0:cell_active_fe_index[row]];
    }



    inline
    unsigned int
    DoFInfo::fe_index_from_degree (const unsigned int fe_degree) const
    {
      const unsigned int n_indices = fe_index_conversion.size();
      for (unsigned int i=0; i<n_indices; ++i)
        if (fe_index_conversion[i].first == fe_degree)
          return i;
      return n_indices;
    }



    inline
    unsigned int
    DoFInfo::fe_index_from_dofs_per_cell (const unsigned int dofs_per_cell) const
    {
      for (unsigned int i=0; i<fe_index_conversion.size(); ++i)
        if (fe_index_conversion[i].second == dofs_per_cell)
          return i;
      return 0;
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

#endif  // ifndef DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
