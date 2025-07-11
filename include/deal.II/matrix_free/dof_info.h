// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_dof_info_h
#define dealii_matrix_free_dof_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/shape_info.h>

#include <array>
#include <memory>


DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN

// forward declarations

namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int dim>
    class HangingNodes;

    struct TaskInfo;

    template <typename Number>
    struct ConstraintValues;

    namespace VectorDataExchange
    {
      class Base;
    }
  } // namespace MatrixFreeFunctions
} // namespace internal

template <typename>
class AffineConstraints;

class DynamicSparsityPattern;

template <typename>
class TriaIterator;

template <int, int, bool>
class DoFCellAccessor;

namespace Utilities
{
  namespace MPI
  {
    class Partitioner;
  }
} // namespace Utilities

#endif

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * Type of the 8-bit representation of the refinement configuration that
     * is in hanging_nodes_internal.h.
     */
    using compressed_constraint_kind = std::uint8_t;

    /**
     * The class that stores the indices of the degrees of freedom for all the
     * cells. Essentially, this is a smart number cache in the style of a
     * DoFHandler that also embeds the description of constraints directly on
     * the cell level without the need to refer to the external
     * AffineConstraints object.
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
     * @ingroup matrixfree
     */
    struct DoFInfo
    {
      /**
       * This value is used to define subranges in the vectors which we can
       * zero inside the MatrixFree::loop() call. The goal is to only clear a
       * part of the vector at a time to keep the values that are zeroed in
       * caches, saving one global vector access for the case where this is
       * applied rather than `vector = 0.;`.
       *
       * We set the granularity to 64 - that is a number sufficiently large
       * to minimize loop peel overhead within the work (and compatible with
       * vectorization lengths of up to 16) and small enough to not waste on
       * the size of the individual chunks.
       */
      static const unsigned int chunk_size_zero_vector = 64;

      /**
       * Default empty constructor.
       */
      DoFInfo();

      /**
       * Copy constructor.
       */
      DoFInfo(const DoFInfo &) = default;

      /**
       * Move constructor.
       */
      DoFInfo(DoFInfo &&) noexcept = default;

      /**
       * Destructor.
       */
      ~DoFInfo() = default;

      /**
       * Copy assignment operator.
       */
      DoFInfo &
      operator=(const DoFInfo &) = default;

      /**
       * Move assignment operator.
       */
      DoFInfo &
      operator=(DoFInfo &&) noexcept = default;

      /**
       * Clear all data fields in this class.
       */
      void
      clear();

      /**
       * Return the FE index for a given finite element degree. If not in
       * hp-mode, this function always returns index 0. If an index is not found
       * in hp-mode, it returns numbers::invalid_unsigned_int.
       */
      unsigned int
      fe_index_from_degree(const unsigned int first_selected_component,
                           const unsigned int fe_degree) const;

      /**
       * Populate the vector @p local_indices with locally owned degrees of freedom
       * stored on the cell batch @p cell_batch.
       * If @p with_constraints is `true`, then the returned vector will contain indices
       * required to resolve constraints.
       *
       * The image below illustrates the output of this function for cell blocks
       * zero and one with zero Dirichlet boundary conditions at the bottom of
       * the domain. Note that due to the presence of constraints, the DoFs
       * returned by this function for the case `with_constraints = true` are
       * not a simple union
       * of per cell DoFs on the cell block @p cell.
       *
       * @image html dofinfo_get_dof_indices.png
       *
       * @note The returned indices may contain duplicates. The unique set can be
       * obtain using `std::sort()` followed by `std::unique()` and
       * `std::vector::erase()`.
       */
      void
      get_dof_indices_on_cell_batch(std::vector<unsigned int> &local_indices,
                                    const unsigned int         cell_batch,
                                    const bool with_constraints = true) const;

      /**
       * This internal method takes the local indices on a cell (two versions:
       * hanging-node constraints resolved if possible and plain, i.e., not
       * resolved) and fills them into this class. It resolves the constraints
       * and distributes the results. Ghost indices, i.e., indices that are
       * located on another processor, get a temporary number by this function,
       * and will later be assigned the correct index after all the ghost
       * indices have been collected by the call to @p assign_ghosts.
       */
      template <typename number>
      void
      read_dof_indices(
        const std::vector<types::global_dof_index> &local_indices_resolved,
        const std::vector<types::global_dof_index> &local_indices,
        const bool                                  cell_has_hanging_nodes,
        const dealii::AffineConstraints<number>    &constraints,
        const unsigned int                          cell_number,
        ConstraintValues<double>                   &constraint_values,
        bool                                       &cell_at_boundary);

      /**
       * For a given cell, determine if it has hanging node constraints. If yes,
       * adjust the dof indices, store the mask, and return true as indication.
       */
      template <int dim>
      bool
      process_hanging_node_constraints(
        const HangingNodes<dim>                      &hanging_nodes,
        const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
        const unsigned int                            cell_number,
        const TriaIterator<DoFCellAccessor<dim, dim, false>> &cell,
        std::vector<types::global_dof_index>                 &dof_indices);

      /**
       * This method assigns the correct indices to ghost indices from the
       * temporary numbering employed by the @p read_dof_indices function. The
       * numbers are localized with respect to the MPI process, and ghosts
       * start at the end of the locally owned range. This way, we get direct
       * access to all vector entries.
       */
      void
      assign_ghosts(const std::vector<unsigned int> &boundary_cells,
                    const MPI_Comm                   communicator_sm,
                    const bool use_vector_data_exchanger_full);

      /**
       * This method reorders the way cells are gone through based on a given
       * renumbering of the cells. It also takes @p vectorization_length cells
       * together and interprets them as one cell only, as is needed for
       * vectorization.
       */
      void
      reorder_cells(const TaskInfo                   &task_info,
                    const std::vector<unsigned int>  &renumbering,
                    const std::vector<unsigned int>  &constraint_pool_row_index,
                    const std::vector<unsigned char> &irregular_cells);

      /**
       * Finds possible compression for the cell indices that we can apply for
       * increased efficiency. Run at the end of reorder_cells.
       */
      void
      compute_cell_index_compression(
        const std::vector<unsigned char> &irregular_cells);

      /**
       * Finds possible compression for the face indices that we can apply for
       * increased efficiency. Run at the end of reorder_cells.
       */
      template <int length>
      void
      compute_face_index_compression(
        const std::vector<FaceToCellTopology<length>> &faces,
        bool hold_all_faces_to_owned_cells);

      /**
       * This function computes the connectivity of the currently stored
       * indices in terms of connections between the individual cells and
       * fills the structure into a sparsity pattern.
       */
      void
      make_connectivity_graph(const TaskInfo                  &task_info,
                              const std::vector<unsigned int> &renumbering,
                              DynamicSparsityPattern &connectivity) const;

      /**
       * In case face integrals are enabled, find out whether certain loops
       * over the unknowns only access a subset of all the ghost dofs we keep
       * in the main partitioner.
       */
      void
      compute_tight_partitioners(
        const Table<2, ShapeInfo<double>>        &shape_info,
        const unsigned int                        n_owned_cells,
        const unsigned int                        n_lanes,
        const std::vector<FaceToCellTopology<1>> &inner_faces,
        const std::vector<FaceToCellTopology<1>> &ghosted_faces,
        const bool                                fill_cell_centric,
        const MPI_Comm                            communicator_sm,
        const bool use_vector_data_exchanger_full);

      /**
       * Given @p cell_indices_contiguous_sm containing the local index of
       * cells of face batches (inner/outer) and cell batches compute
       * dof_indices_contiguous_sm.
       */
      void
      compute_shared_memory_contiguous_indices(
        std::array<std::vector<std::pair<unsigned int, unsigned int>>, 3>
          &cell_indices_contiguous_sm);

      /**
       * Compute a renumbering of the degrees of freedom to improve the data
       * access patterns for this class that can be utilized by the categories
       * in the IndexStorageVariants enum. For example, the index ordering can
       * be improved for typical DG elements by interleaving the degrees of
       * freedom from batches of cells, which avoids the explicit data
       * transposition in IndexStorageVariants::contiguous. Currently, these
       * more advanced features are not implemented, so there is only limited
       * value of this function.
       */
      void
      compute_dof_renumbering(
        std::vector<types::global_dof_index> &renumbering);

      /**
       * Fills the array that defines how to zero selected ranges in the result
       * vector within the cell loop, filling the two member variables @p
       * vector_zero_range_list_index and @p vector_zero_range_list.
       *
       * The intent of this pattern is to zero the vector entries in close
       * temporal proximity to the first access and thus keeping the vector
       * entries in cache.
       */
      template <int length>
      void
      compute_vector_zero_access_pattern(
        const TaskInfo                                &task_info,
        const std::vector<FaceToCellTopology<length>> &faces);

      /**
       * Return the memory consumption in bytes of this class.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Prints a detailed summary of memory consumption in the different
       * structures of this class to the given output stream.
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType     &out,
                               const TaskInfo &size_info) const;

      /**
       * Prints a representation of the indices in the class to the given
       * output stream.
       */
      template <typename Number>
      void
      print(const std::vector<Number>       &constraint_pool_data,
            const std::vector<unsigned int> &constraint_pool_row_index,
            std::ostream                    &out) const;

      /**
       * Enum for various storage variants of the indices. This storage format
       * is used to implement more efficient indexing schemes in case the
       * underlying data structures allow for them, and to inform the access
       * functions in FEEvaluationBase::read_write_operation() on which array
       * to get the data from. One example of more efficient storage is the
       * enum value IndexStorageVariants::contiguous, which means that one can
       * get the indices to all degrees of freedom of a cell by reading only
       * the first index for each cell, whereas all subsequent indices are
       * merely an offset from the first index.
       */
      enum class IndexStorageVariants : unsigned char
      {
        /**
         * This value indicates that no index compression was found and the
         * only valid storage is to access all indices present on the cell,
         * possibly including constraints. For a cell/face of this index type,
         * the data access in FEEvaluationBase is directed to the array @p
         * dof_indices with the index
         * `row_starts[cell_index*n_vectorization*n_components].first`.
         */
        full,
        /**
         * This value indicates that the indices are interleaved for access
         * with vectorized gather and scatter operation. This storage variant
         * is possible in case there are no constraints on the cell and the
         * indices in the batch of cells are not pointing to the same global
         * index in different slots of a vectorized array (in order to support
         * scatter operations). For a cell/face of this index type, the data
         * access in FEEvaluationBase is directed to the array
         * `dof_indices_interleaved` with the index
         * `row_starts[cell_index*n_vectorization*n_components].first`.
         */
        interleaved,
        /**
         * This value indicates that the indices within a cell are all
         * contiguous, and one can get the index to the cell by reading that
         * single value for each of the cells in the cell batch. For a
         * cell/face of this index type, the data access in FEEvaluationBase
         * is directed to the array `dof_indices_contiguous` with the index
         * `cell_index*n_vectorization*n_components`.
         */
        contiguous,
        /**
         * This value indicates that the indices with a cell are contiguous and
         * interleaved for vectorization, i.e., the first DoF index on a cell
         * to the four or eight cells in the vectorization batch come first,
         * than the second DoF index, and so on. Furthermore, the interleaving
         * between cells implies that only the batches for vectorization can be
         * accessed efficiently, whereas there is a strided access for getting
         * only some of the entries.
         *
         * The two additional categories `interleaved_contiguous_strided` and
         * `interleaved_contiguous_mixed_strides` are a consequence of this
         * storage type. The former is for faces where at least one of the two
         * adjacent sides will break with the interleaved storage. We then have
         * to make a strided access as described in the next category. The last
         * category `interleaved_contiguous_mixed_strides` appears in the ghost
         * layer, see the more detailed description of that category below.
         * Again, this is something that cannot be avoided in general once we
         * interleave the indices between cells.
         *
         * For a cell/face of this index type, the data access in
         * FEEvaluationBase is directed to the array `dof_indices_contiguous`
         * with the index `cell_index*n_vectorization*n_components`.
         */
        interleaved_contiguous,
        /**
         * Similar to interleaved_contiguous storage, but for the case when the
         * interleaved indices are only contiguous within the degrees of
         * freedom, but not also over the components of a vectorized array.
         * This happens typically on faces with DG where the cells have
         * `interleaved_contiguous` storage but the faces' numbering is not the
         * same as the cell's numbering. For a
         * cell/face of this index type, the data access in FEEvaluationBase
         * is directed to the array `dof_indices_contiguous` with the index
         * `cell_index*n_vectorization*n_components`.
         */
        interleaved_contiguous_strided,
        /**
         * Similar to interleaved_contiguous_separate storage, but for the case
         * when the interleaved indices are not `n_vectorization apart`. This
         * happens typically within the ghost layer of DG where the remote
         * owner has applied an interleaved storage and the current processor
         * only sees some of the cells. For a
         * cell/face of this index type, the data access in FEEvaluationBase
         * is directed to the array `dof_indices_contiguous` with the index
         * `cell_index*n_vectorization*n_components`, including the array
         * `dof_indices_interleave_strides` for the information about the
         * actual stride.
         */
        interleaved_contiguous_mixed_strides
      };

      /**
       * Enum used to distinguish the data arrays for the vectorization type
       * in cells and faces.
       */
      enum DoFAccessIndex : unsigned char
      {
        /**
         * The data index for the faces designated as interior
         */
        dof_access_face_interior = 0,
        /**
         * The data index for the faces designated as exterior
         */
        dof_access_face_exterior = 1,
        /**
         * The data index for the cells
         */
        dof_access_cell = 2
      };

      /**
       * For efficiency reasons, always keep a fixed number of cells with
       * similar properties together. This variable controls the number of
       * cells batched together. As opposed to the other classes which are
       * templated on the number type, this class as a pure index container is
       * not templated, so we need to keep the information otherwise contained
       * in VectorizedArrayType::size().
       */
      unsigned int vectorization_length;

      /**
       * Stores the index storage variant of all cell and face batches.
       *
       * The three arrays given here address the types for the faces decorated
       * as interior (0), the faces decorated with as exterior (1), and the
       * cells (2) according to CellOrFaceAccess.
       */
      std::array<std::vector<IndexStorageVariants>, 3> index_storage_variants;

      /**
       * Stores the rowstart indices of the compressed row storage in the @p
       * dof_indices and @p constraint_indicator fields. These two fields are
       * always accessed together, so it is simpler to keep just one variable
       * for them. This also obviates keeping two rowstart vectors in sync.
       */
      std::vector<std::pair<unsigned int, unsigned int>> row_starts;

      /**
       * Stores the indices of the degrees of freedom for each cell. These
       * indices are computed in MPI-local index space, i.e., each processor
       * stores the locally owned indices as numbers between <tt>0</tt> and
       * <tt>n_locally_owned_dofs-1</tt> and ghost indices in the range
       * <tt>n_locally_owned_dofs</tt> to
       * <tt>n_locally_owned_dofs+n_ghost_dofs</tt>. The translation between
       * this MPI-local index space and the global numbering of degrees of
       * freedom is stored in the @p vector_partitioner data structure.  This
       * array also includes the indirect contributions from constraints,
       * which are described by the @p constraint_indicator field. Because of
       * variable lengths of rows, this would be a vector of a vector.
       * However, we use one contiguous memory region and store the rowstart
       * in the variable @p row_starts.
       */
      std::vector<unsigned int> dof_indices;

      /**
       * Supported components of all entries of the hp::FECollection object of
       * the given DoFHandler.
       */
      std::vector<std::vector<bool>> hanging_node_constraint_masks_comp;

      /**
       * Masks indicating for each cell and component if the optimized
       * hanging-node constraint is applicable and if yes which type.
       */
      std::vector<compressed_constraint_kind> hanging_node_constraint_masks;

      /**
       * This variable describes the position of constraints in terms of the
       * local numbering of degrees of freedom on a cell. The first number
       * stores the distance from one constrained degree of freedom to the
       * next. This allows to identify the position of constrained DoFs as we
       * loop through the local degrees of freedom of the cell when reading
       * from or writing to a vector. The second number stores the index of
       * the constraint weights, stored in the variable constraint_pool_data.
       */
      std::vector<std::pair<unsigned short, unsigned short>>
        constraint_indicator;

      /**
       * Reordered index storage for `IndexStorageVariants::interleaved`.
       */
      std::vector<unsigned int> dof_indices_interleaved;

      /**
       * Compressed index storage for faster access than through @p
       * dof_indices used according to the description in IndexStorageVariants.
       *
       * The three arrays given here address the types for the faces decorated
       * as interior (0), the faces decorated with as exterior (1), and the
       * cells (2) according to CellOrFaceAccess.
       */
      std::array<std::vector<unsigned int>, 3> dof_indices_contiguous;

      /**
       * The same as above but for shared-memory usage. The first value of the
       * pair is identifying the owning process and the second the index
       * within that locally-owned data of that process.
       *
       * @note This data structure is only set up if all entries in
       *   index_storage_variants[2] are IndexStorageVariants::contiguous.
       */
      std::array<std::vector<std::pair<unsigned int, unsigned int>>, 3>
        dof_indices_contiguous_sm;

      /**
       * Compressed index storage for faster access than through @p
       * dof_indices used according to the description in IndexStorageVariants.
       *
       * The three arrays given here address the types for the faces decorated
       * as minus (0), the faces decorated with as plus (1), and the cells
       * (2).
       */
      std::array<std::vector<unsigned int>, 3> dof_indices_interleave_strides;

      /**
       * Caches the number of indices filled when vectorizing. This
       * information can implicitly deduced from the row_starts data fields,
       * but this field allows for faster access.
       *
       * The three arrays given here address the types for the faces decorated
       * as interior (0), the faces decorated with as exterior (1), and the
       * cells (2) according to CellOrFaceAccess.
       */
      std::array<std::vector<unsigned char>, 3> n_vectorization_lanes_filled;

      /**
       * This stores the parallel partitioning that can be used to set up
       * vectors. The partitioner includes the description of the local range
       * in the vector, and also includes how the ghosts look like. This
       * enables initialization of vectors based on the DoFInfo field.
       */
      std::shared_ptr<const Utilities::MPI::Partitioner> vector_partitioner;

      /**
       * Vector exchanger compatible with vector_partitioner.
       */
      std::shared_ptr<
        const internal::MatrixFreeFunctions::VectorDataExchange::Base>
        vector_exchanger;

      /**
       * Vector exchanger compatible with partitioners that select a subset of
       * ghost indices to the full
       * vector partitioner stored in @p vector_partitioner. These
       * partitioners are used in specialized loops that only import parts of
       * the ghosted region for reducing the amount of communication. There
       * are five variants of the partitioner initialized:
       * - one that queries only the cell values,
       * - one that additionally describes the indices for
       *   evaluating the function values on relevant faces,
       * - one that describes the indices for evaluation both the function
       *   values and the gradients on relevant faces adjacent to the locally
       *   owned cells,
       * - one that additionally describes the indices for
       *   evaluating the function values on all faces, and
       * - one that describes the indices for evaluation both the function
       *   values and the gradients on all faces adjacent to the locally owned
       *   cells.
       */
      std::array<
        std::shared_ptr<
          const internal::MatrixFreeFunctions::VectorDataExchange::Base>,
        5>
        vector_exchanger_face_variants;

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
       * Stores the offset in terms of the number of base elements over all
       * DoFInfo objects.
       */
      unsigned int global_base_element_offset;

      /**
       * Stores the number of base elements in the DoFHandler where the
       * indices have been read from.
       */
      unsigned int n_base_elements;

      /**
       * Stores the number of components of each base element in the finite
       * element where the indices have been read from.
       */
      std::vector<unsigned int> n_components;

      /**
       * The ith entry of this vector stores the component number of the given
       * base element.
       */
      std::vector<unsigned int> start_components;

      /**
       * For a given component in an FESystem, this variable tells which base
       * element the index belongs to.
       */
      std::vector<unsigned int> component_to_base_index;

      /**
       * For a vector-valued element, this gives the constant offset in the
       * number of degrees of freedom starting at the given component, as the
       * degrees are numbered by degrees of freedom. This data structure does
       * not take possible constraints and thus, shorter or longer lists, into
       * account. This information is encoded in the row_starts variables
       * directly.
       *
       * The outer vector goes through the various FE indices in the hp-case,
       * similarly to the @p dofs_per_cell variable.
       */
      std::vector<std::vector<unsigned int>> component_dof_indices_offset;

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
       * Stores the index of the active finite element in the hp-case.
       */
      std::vector<unsigned int> cell_active_fe_index;

      /**
       * Stores the maximum degree of different finite elements for the
       * hp-case.
       */
      unsigned int max_fe_index;

      /**
       * To each of the slots in an hp-adaptive case, the inner vector stores
       * the corresponding element degree. This is used by the constructor of
       * FEEvaluationBase to identify the correct data slot in the hp-case.
       */
      std::vector<std::vector<unsigned int>> fe_index_conversion;

      /**
       * Temporarily stores the numbers of ghosts during setup. Cleared when
       * calling @p assign_ghosts. Then, all information is collected by the
       * partitioner.
       */
      std::vector<types::global_dof_index> ghost_dofs;

      /**
       * Stores an integer to each partition in TaskInfo that indicates
       * whether to clear certain parts in the result vector if the user
       * requested it with the respective argument in the MatrixFree::loop.
       */
      std::vector<unsigned int> vector_zero_range_list_index;

      /**
       * Stores the actual ranges in the vector to be cleared.
       */
      std::vector<std::pair<unsigned int, unsigned int>> vector_zero_range_list;

      /**
       * Stores an integer to each partition in TaskInfo that indicates when
       * to schedule operations that will be done before any access to vector
       * entries.
       */
      std::vector<unsigned int> cell_loop_pre_list_index;

      /**
       * Stores the actual ranges of the operation before any access to vector
       * entries.
       */
      std::vector<std::pair<unsigned int, unsigned int>> cell_loop_pre_list;

      /**
       * Stores an integer to each partition in TaskInfo that indicates when
       * to schedule operations that will be done after all access to vector
       * entries.
       */
      std::vector<unsigned int> cell_loop_post_list_index;

      /**
       * Stores the actual ranges of the operation after all access to vector
       * entries.
       */
      std::vector<std::pair<unsigned int, unsigned int>> cell_loop_post_list;
    };


    /*-------------------------- Inline functions ---------------------------*/

#ifndef DOXYGEN


    inline unsigned int
    DoFInfo::fe_index_from_degree(const unsigned int first_selected_component,
                                  const unsigned int fe_degree) const
    {
      const unsigned int n_indices = fe_index_conversion.size();
      if (n_indices <= 1)
        return 0;
      for (unsigned int i = 0; i < n_indices; ++i)
        if (fe_index_conversion[i][first_selected_component] == fe_degree)
          return i;
      return numbers::invalid_unsigned_int;
    }

#endif // ifndef DOXYGEN

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
