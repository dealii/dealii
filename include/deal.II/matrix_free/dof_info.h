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


#ifndef dealii_matrix_free_dof_info_h
#define dealii_matrix_free_dof_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/task_info.h>

#include <array>
#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    // Forward declaration
#ifndef DOXYGEN
    template <typename Number>
    struct ConstraintValues;
#endif

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
     *
     * @author Katharina Kormann and Martin Kronbichler, 2010, 2011, 2018
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
       * Clear all data fields in this class.
       */
      void
      clear();

      /**
       * Return the FE index for a given finite element degree. If not in hp
       * mode, this function always returns index 0. If an index is not found
       * in hp mode, it returns numbers::invalid_unsigned_int.
       */
      unsigned int
      fe_index_from_degree(const unsigned int first_selected_component,
                           const unsigned int fe_degree) const;

      /**
       * Populate the vector @p locall_indices with locally owned degrees of freedom
       * stored on the cell block @p cell.
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
      get_dof_indices_on_cell_batch(std::vector<unsigned int> &locall_indices,
                                    const unsigned int         cell,
                                    const bool with_constraints = true) const;

      /**
       * This internal method takes the local indices on a cell and fills them
       * into this class. It resolves the constraints and distributes the
       * results. Ghost indices, i.e., indices that are located on another
       * processor, get a temporary number by this function, and will later be
       * assigned the correct index after all the ghost indices have been
       * collected by the call to @p assign_ghosts.
       */
      template <typename number>
      void
      read_dof_indices(
        const std::vector<types::global_dof_index> &local_indices,
        const std::vector<unsigned int> &           lexicographic_inv,
        const AffineConstraints<number> &           constraints,
        const unsigned int                          cell_number,
        ConstraintValues<double> &                  constraint_values,
        bool &                                      cell_at_boundary);

      /**
       * This method assigns the correct indices to ghost indices from the
       * temporary numbering employed by the @p read_dof_indices function. The
       * numbers are localized with respect to the MPI process, and ghosts
       * start at the end of the locally owned range. This way, we get direct
       * access to all vector entries.
       */
      void
      assign_ghosts(const std::vector<unsigned int> &boundary_cells);

      /**
       * This method reorders the way cells are gone through based on a given
       * renumbering of the cells. It also takes @p vectorization_length cells
       * together and interprets them as one cell only, as is needed for
       * vectorization.
       */
      void
      reorder_cells(const TaskInfo &                  task_info,
                    const std::vector<unsigned int> & renumbering,
                    const std::vector<unsigned int> & constraint_pool_row_index,
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
        const std::vector<FaceToCellTopology<length>> &faces);

      /**
       * This function computes the connectivity of the currently stored
       * indices in terms of connections between the individual cells and
       * fills the structure into a sparsity pattern.
       */
      void
      make_connectivity_graph(const TaskInfo &                 task_info,
                              const std::vector<unsigned int> &renumbering,
                              DynamicSparsityPattern &connectivity) const;

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
        const TaskInfo &                               task_info,
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
      print_memory_consumption(StreamType &    out,
                               const TaskInfo &size_info) const;

      /**
       * Prints a representation of the indices in the class to the given
       * output stream.
       */
      template <typename Number>
      void
      print(const std::vector<Number> &      constraint_pool_data,
            const std::vector<unsigned int> &constraint_pool_row_index,
            std::ostream &                   out) const;

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
       * Stores the dimension of the underlying DoFHandler. Since the indices
       * are not templated, this is the variable that makes the dimension
       * accessible in the (rare) cases it is needed inside this class.
       */
      unsigned int dimension;

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
      std::vector<IndexStorageVariants> index_storage_variants[3];

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
      std::vector<unsigned int> dof_indices_contiguous[3];

      /**
       * Compressed index storage for faster access than through @p
       * dof_indices used according to the description in IndexStorageVariants.
       *
       * The three arrays given here address the types for the faces decorated
       * as minus (0), the faces decorated with as plus (1), and the cells
       * (2).
       */
      std::vector<unsigned int> dof_indices_interleave_strides[3];

      /**
       * Caches the number of indices filled when vectorizing. This
       * information can implicitly deduced from the row_starts data fields,
       * but this field allows for faster access.
       *
       * The three arrays given here address the types for the faces decorated
       * as interior (0), the faces decorated with as exterior (1), and the
       * cells (2) according to CellOrFaceAccess.
       */
      std::vector<unsigned char> n_vectorization_lanes_filled[3];

      /**
       * This stores the parallel partitioning that can be used to set up
       * vectors. The partitioner includes the description of the local range
       * in the vector, and also includes how the ghosts look like. This
       * enables initialization of vectors based on the DoFInfo field.
       */
      std::shared_ptr<const Utilities::MPI::Partitioner> vector_partitioner;

      /**
       * This partitioning selects a subset of ghost indices to the full
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
      std::array<std::shared_ptr<const Utilities::MPI::Partitioner>, 5>
        vector_partitioner_face_variants;

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
       * The outer vector goes through the various fe indices in the hp case,
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
       * Stores the index of the active finite element in the hp case.
       */
      std::vector<unsigned int> cell_active_fe_index;

      /**
       * Stores the maximum degree of different finite elements for the hp
       * case.
       */
      unsigned int max_fe_index;

      /**
       * To each of the slots in an hp adaptive case, the inner vector stores
       * the corresponding element degree. This is used by the constructor of
       * FEEvaluationBase to identify the correct data slot in the hp case.
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
