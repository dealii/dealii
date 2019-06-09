// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_block_csr_matrix_h
#define dealii_block_csr_matrix_h

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/blas_extension_templates.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector_operations_internal.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/base/array_view.h>

// we will need to serialize pairs and maps
#include <boost/serialization/utility.hpp>
#include <boost/serialization/map.hpp>

#include <iostream>

DEAL_II_NAMESPACE_OPEN
using namespace dealii;

  // forward declaration
  template <typename NumberType>
  class BlockCSRMatrix;


  /**
   * A class similar to internal::MatrixFreeFunctions::DoFInfo
   * that stores and precomputes certain DoF information to speedup
   * application of FE cell operators to BlockCSRMatrix.
   */
  struct DoFInfo
  {
    /**
     * Initialize internal data structures
     */
    template <int dim, typename NumberType>
    void initialize(
      const DoFHandler<dim> &dof_handler,
      const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> &partitioner,
      const std::shared_ptr<const MatrixFree<dim,NumberType>> &matrix_free,
      const std::shared_ptr<const BlockIndices> &row_blocks,
      const std::vector<unsigned int> &cell_index_permutation);

    /**
     * Initialize this object for a single _cell_ with rows @p my_rows.
     */
    void initialize(const std::vector<unsigned int> &my_rows,
                    const std::shared_ptr<const BlockIndices> &row_blocks);

    /**
     * Clear all data fields in this class.
     */
    void clear();

    /**
     * Stores the rowstart indices of the compressed row storage in the dof_indices and block_indices fields
     * for each cell.
     */
    std::vector<std::pair<unsigned int, unsigned int>> row_starts;

    /**
     * Stores the indices of the degrees of freedom for each cell.
     *
     * The first element is the row index within the block, whereas the second one
     * is the index of this degree of freedom within the cell.
     */
    std::vector<std::pair<unsigned int, unsigned int>> dof_indices;

    /**
     * stores the row block number and the number of elements within this
     * block for all cells.
     */
    std::vector<std::pair<unsigned int, unsigned int>> block_indices;

  private:
    using TUP = std::tuple<unsigned int, unsigned int, unsigned int>;

    /**
     * Internal function to process DoFs stored in @p cell_dof_indices_local
     * by sorting and adding to the data structures.
     */
    void process_cell_indices(std::vector<TUP> &cell_dof_indices_local);
  };

  /**
   * A namespace in which we declare iterators over the elements of block sparse
   * matrices.
   */
  namespace BlockCSRMatrixIterators
  {
    /**
     * Declare type for container size.
     */
    using size_type = types::global_dof_index;

    // forward declaration
    template <typename NumberType, bool Constness>
    class Iterator;

    /**
     * General accessor class for matrices. The first template
     * argument denotes the underlying number type, the second the constness of
     * the matrix.
     *
     * This class builds on the accessor classes used for sparsity patterns
     * to loop over all nonzero blocks, and only adds the accessor functions to
     * gain access to the pointer representing the beginning of the non-empty block.
     */
    template <typename NumberType, bool Constness>
    class Accessor : public SparsityPatternIterators::Accessor
    {
    public:
      /**
       * Typedef for the type (including constness) of the matrix to be used
       * here.
       */
      using MatrixType = typename std::conditional<Constness,
                                  const BlockCSRMatrix<NumberType>,
                                  BlockCSRMatrix<NumberType>>::type;

      /*
       * Pointer to the row data (including constness)
       */
      using pointer = typename std::conditional<Constness,
                                const NumberType *,
                                NumberType *>::type;

      /**
       * Constructor.
       */
      Accessor(MatrixType *matrix,
               const std::size_t index_within_matrix);

      /**
       * Constructor. Construct the end accessor for the given matrix.
       */
      Accessor(MatrixType *matrix);

      /**
       * Default constructor.
       *
       * This constructor is provided only to be able to store accessors in STL
       * containers such as `std::vector`.
       */
      Accessor();

      /**
       * Copy constructor to get from a non-const accessor to a const accessor.
       */
      Accessor(const Accessor<NumberType, false> &a);

      /**
       * Pointer to the beginning of the current block
       */
      pointer data() const;

      /**
       * Return a reference to the matrix into which this accessor points. Note
       * that in the present case, this is a constant reference.
       */
      MatrixType &get_matrix() const;

    private:
      /**
       * Pointer to the matrix we use.
       */
      MatrixType *matrix;

      /**
       * Pointer to the actual data
       */
      pointer current_data;

      /**
       * Move the accessor to the next non-empty block in the block matrix.
       */
      void advance();

      /**
       * Make iterator class a friend.
       */
      template <typename, bool>
      friend class Iterator;
    };

    /**
     * Iterator for constant and non-constant matrices.
     *
     * The typical use for these iterators is to iterate over the elements of a
     * sparse matrix or over the elements of individual rows.
     *
     * The first template argument denotes the underlying number type, the
     * second the constness of the matrix.
     */
    template <typename NumberType, bool Constness>
    class Iterator
    {
    public:
      /**
       * Typedef for the matrix type (including constness) we are to operate on.
       */
      using MatrixType = typename Accessor<NumberType, Constness>::MatrixType;

      /**
       * An alias for the type you get when you dereference an iterator of the
       * current kind.
       */
      using value_type = const Accessor<NumberType, Constness> &;

      /**
       * Constructor. Create an iterator into the matrix @p matrix for the given
       * @p row and @p index_within_row (counting from the first entry in this
       * row).
       */
      Iterator(MatrixType *matrix,
               const std::size_t index_within_matrix);

      /**
       * Constructor. Create the end iterator for the given matrix.
       */
      Iterator(MatrixType *matrix);

      /**
       * Default constructor.
       *
       * This constructor is provided only to be able to store accessors in STL
       * containers such as `std::vector`.
       */
      Iterator();

      /**
       * Conversion constructor to get from a non-const iterator to a const
       * iterator.
       */
      Iterator(const Iterator<NumberType, false> &i);

      /**
       * Prefix increment.
       */
      Iterator &operator++();

      /**
       * Postfix increment.
       */
      Iterator operator++(int);

      /**
       * Dereferencing operator.
       */
      const Accessor<NumberType, Constness> &operator*() const;

      /**
       * Dereferencing operator.
       */
      const Accessor<NumberType, Constness> *operator->() const;

      /**
       * Comparison. True, if both iterators point to the same matrix position.
       */
      bool operator==(const Iterator &) const;

      /**
       * Inverse of <tt>==</tt>.
       */
      bool operator!=(const Iterator &) const;

      /**
       * Comparison operator. Result is true if either the first row number is
       * smaller or if the row numbers are equal and the first index is smaller.
       *
       * This function is only valid if both iterators point into the same
       * matrix.
       */
      bool operator<(const Iterator &) const;

      /**
       * Comparison operator. Works in the same way as above operator, just the
       * other way round.
       */
      bool operator>(const Iterator &) const;

      /**
       * Return an iterator that is @p n ahead of the current one.
       */
      Iterator operator+(const size_type n) const;

    private:
      /**
       * Store an object of the accessor class.
       */
      Accessor<NumberType, Constness> accessor;
    };



    /**
     * A class to provide read/write access
     * to certain rows in BlockCSRMatrix. The iteration over
     * columns is done via the advance() method.
     */
    template <typename NumberType, bool Constness>
    class RowsBlockAccessor
    {
    public:
      /**
       * Typedef for the matrix type (including constness) we are to operate on.
       */
      using MatrixType = typename Accessor<NumberType, Constness>::MatrixType;

      /**
       * A pointer to vectorized array (including constness)
       */
      using vectorized_pointer =
        typename std::conditional<Constness,
                                  const VectorizedArray<NumberType> *,
                                  VectorizedArray<NumberType> *>::type;

      /**
       * A constructor to give access to the BCSR matrix.
       */
      RowsBlockAccessor(MatrixType *matrix, const DoFInfo& dof_info);

      /**
       * Return the current block column.
       */
      types::global_dof_index get_current_block_column() const;

      /**
       * Return the block size of the current column.
       */
      types::global_dof_index get_col_block_size() const;

      /**
       * Advance the accessor to the next global column with at least
       * one non-zero element among the chosen rows.
       *
       * The function will return the new column block or
       * numbers::invalid_size_type if all column blocks have been processed.
       */
      types::global_dof_index advance();

      /**
       * Reinitialize this object to provide access to dofs/rows
       * stored on the @p subcell within the Matrix-free macro cell @p cell.
       */
      types::global_dof_index reinit(const unsigned int cell,
                                     const unsigned int subcell = 0);

      /**
       * Clear content of this object. The state is the same as the one after
       * calling the default constructor.
       */
      void clear();

      /**
       * Return number of row blocks used to access rows/DoFs, provided to
       * reinit() function.
       */
      unsigned int n_row_blocks() const;

      /**
       * Process all active rows for the current block column using @p function
       * with three arguments:
       * - a view to DoFInfo
       * - a pointer to start of the block `VectorizedArray<Number> *`
       * - stride to access vectorized rows
       */
      void process_active_rows_vectorized(const std::function<void(
        const ArrayView<const std::pair<unsigned int, unsigned int>>&,
        vectorized_pointer const block_start,
        const unsigned int stride)>
        &func) const;

    private:

      /**
       * Using the current value of `col_block`,
       * decide on active/disabled row blocks and set
       * pointers in each row within the block.
       */
      void set_pointers_and_active_flags();

      /**
       * A class to store information about rows that belong to the same group.
       */
      struct RowBlock
      {
        /**
         * Current iterator to the block matrix
         */
        Iterator<NumberType, Constness> it;

        /**
         * End of the row.
         */
        Iterator<NumberType, Constness> end;

        /**
         * Vectorized array pointer to the beginning of this block.
         */
        vectorized_pointer pointer;

        /**
         * A flag to indicate whether or this block is active,
         * meaning `it!=end` and `it->column()` matches the current
         * column block in RowsBlockAccessor.
         */
        bool active;

        /**
         * Row block number
         */
        unsigned int block_row;

        /**
         * A view to DoFs of this block, stored in DoFInfo::block_indices
         */
        ArrayView<const std::pair<unsigned int, unsigned int>> dof_view;
      };

      /**
       * A vector which collects row blocks
       */
      std::vector<RowBlock> row_blocks;

      /**
       * Pointer to the matrix we use.
       */
      MatrixType *matrix;

      /**
       * Current block column.
       */
      types::global_dof_index col_block;

      /**
       * Block size for the current column.
       */
      types::global_dof_index col_block_size;

      /**
       * Possibly padded block size for the current column.
       */
      types::global_dof_index col_block_size_p;

      /**
       * Stride size for vectorized pointer.
       */
      unsigned int stride;

      /**
       * Auxiliary class that stores dof information
       */
      const DoFInfo& dof_info;
    };



    /**
     * A class to provide relatively fast read/write access
     * to certain rows in BlockCSRMatrix. The iteration over
     * columns is done via the advance() method.
     */
    template <typename NumberType, bool Constness>
    class RowsAccessorBase
    {
    public:
      /**
       * Typedef for the matrix type (including constness) we are to operate on.
       */
      using MatrixType = typename Accessor<NumberType, Constness>::MatrixType;

      /**
       * Return value at the current global column.
       *
       * @p local_index is allowed to be only one of those indices provided in
       * the constructor.
       */
      NumberType local_element(const unsigned int local_index) const;

      /**
       * Same as above but for global DoF/row @p index.
       */
      NumberType operator()(const types::global_dof_index index) const;

      /**
       * Return the current global column.
       */
      typename MatrixType::size_type current_column() const;

      /**
       * Advance the accessor to the next global column with at least
       * one non-zero element among the desired rows.
       *
       * The function will return `false` if all columns have been processed.
       * At this point the accessor is in the invalid state and shall not be
       * used.
       */
      bool advance();

      /**
       * Provide access to the column @p col.
       *
       * @p col should be larger than the current active column.
       */
      void advance(const unsigned int col);

      /**
       * Reinitialize this object to provide access to rows specified in @p
       * local_rows.
       *
       * @note @p local_rows don't have to be unique. This is intentional
       * to have a better support for internal::MatrixFreeFunctions::DoFInfo.
       */
      void reinit(const std::vector<unsigned int> &local_rows);

      /**
       * Reinitialize this object to provide access to rows within
       * given blocks. @p local_row_blocks should contain unique
       * block numbers.
       */
      void reinit_blocks(const std::vector<unsigned int> &local_row_blocks);

      /**
       * Return number of row blocks used to access rows/DoFs, provided to
       * reinit() function.
       */
      unsigned int n_row_blocks() const;

      /**
       * Return the total number of rows in the matrix.
       */
      unsigned int size() const;

      struct DataType
      {
        /**
         * Current iterator
         */
        Iterator<NumberType, Constness> it;

        /**
         * End of the row
         */
        Iterator<NumberType, Constness> end;

        /**
         * A flag to indicate whether or this struct is active,
         * meaning `it!=end` and `it->column()` matches the current
         * column block.
         */
        bool active;

        /**
         * First row for this block.
         */
        unsigned int block_start;

        /**
         * Row block size of this iterator.
         */
        unsigned int block_size;

        /**
         * Pointer to the data
         */
        typename Accessor<NumberType, Constness>::pointer data;
      };

    protected:
      /**
       * A constructor which will provide access to some rows matrix @p matrix.
       *
       * @note The constructor is protected as we do not want to support
       * creation of objects of the base class.
       */
      RowsAccessorBase(MatrixType *matrix);

      /**
       * Move iterators @p iterators selected by @p active_row_blocks to
       * point to the current block @p block_col .
       */
      void move_iterators();

      /**
       * Pointer to the matrix we use.
       */
      MatrixType *matrix;

      /**
       * A collection of unique row blocks that are active
       */
      std::vector<DataType> active_row_blocks;

      /**
       * Current block column.
       */
      types::global_dof_index block_col;

      /**
       * Current column within the block
       */
      unsigned int col_within_block;
    };

    /**
     * Read and write row accessor.
     */
    template <typename NumberType, bool Constness>
    class RowsAccessor : public RowsAccessorBase<NumberType, Constness>
    {
      /**
       * A constructor which will provide access to matrix @p matrix.
       */
      RowsAccessor(
        typename RowsAccessorBase<NumberType, false>::MatrixType *matrix);

      /**
       * A constructor which will provide access to rows specified in @p
       * local_rows of the matrix @p matrix.
       */
      RowsAccessor(
        typename RowsAccessorBase<NumberType, false>::MatrixType *matrix,
        const std::vector<unsigned int> &local_rows);

      /**
       * Set local row @p local_index to the value @p val.
       *
       * If the matrix sparsity does not hold this value, the function will
       * throw an error in debug mode unless @p val is zero.
       */
      void set_local_element(const unsigned int local_index,
                             const NumberType &val) const;
    };

    /**
     * Read-only row accessor.
     */
    template <typename NumberType>
    class RowsAccessor<NumberType, true>
      : public RowsAccessorBase<NumberType, true>
    {
    public:
      /**
       * Type of the matrix entries. This alias is analogous to
       * <tt>value_type</tt> in the standard library containers.
       */
      using value_type = NumberType;

      /**
       * A constructor which will provide read access to matrix @p matrix.
       */
      RowsAccessor(
        typename RowsAccessorBase<NumberType, true>::MatrixType *matrix);

      /**
       * A constructor which will provide access to rows specified in @p
       * local_rows of the matrix @p matrix.
       *
       * Equivalent to calling a constructor above followed by
       * RowsAccessorBase::reinit().
       */
      RowsAccessor(
        typename RowsAccessorBase<NumberType, true>::MatrixType *matrix,
        const std::vector<unsigned int> &local_rows);
    };

    /**
     * Read and write row accessor.
     */
    template <typename NumberType>
    class RowsAccessor<NumberType, false>
      : public RowsAccessorBase<NumberType, false>
    {
    public:
      /**
       * Type of the matrix entries. This alias is analogous to
       * <tt>value_type</tt> in the standard library containers.
       */
      using value_type = NumberType;

      /**
       * A constructor which will provide read access to matrix @p matrix.
       */
      RowsAccessor(
        typename RowsAccessorBase<NumberType, false>::MatrixType *matrix);

      /**
       * A constructor which will provide access to rows specified in @p
       * local_rows of the matrix @p matrix.
       *
       * Equivalent to calling a constructor above followed by
       * RowsAccessorBase::reinit().
       */
      RowsAccessor(
        typename RowsAccessorBase<NumberType, false>::MatrixType *matrix,
        const std::vector<unsigned int> &local_rows);

      /**
       * Set local row @p local_index to the value @p val.
       *
       * If the matrix sparsity does not hold this value, the function will
       * throw an error in debug mode unless @p val is zero.
       */
      void set_local_element(const unsigned int local_index,
                             const NumberType &val) const;

      /**
       * Same as above but add to the element
       */
      void add_local_element(const unsigned int local_index,
                             const NumberType &val) const;

      /**
       * Same as above but add based on global index
       */
      void add(const types::global_dof_index index,
               const NumberType &val) const;

    private:

      /**
       * Private function that implements addition or setting based on template
       */
      template <bool addition=false>
      void modify_local_element(const unsigned int local_index,
                                const NumberType &val) const;

    };

  } // namespace BlockCSRMatrixIterators



  /**
   * A function that provides a similar functionality to
   * FEEValuation::read_values() but tailored to sparse vectors.
   *
   * @param src_row_accessor read accessor
   * @param c column number within the column block
   * @param n_items number of items to read starting from @p c
   * @param phi FEEValuation class into which the data will be read.
   * @param touched a vector with flags that will be set for zero'th column
   * according to the active rows within the @p src_row_accessor.
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  void
  read_values(const BlockCSRMatrixIterators::RowsBlockAccessor<Number, true>
                &src_row_accessor,
              const unsigned int c,
              const unsigned int n_items,
              FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi,
              std::vector<bool> &touched);

  /**
   * A function that provides a similar functionality to
   * FEEValuation::distribute_local_to_global() but tailored to sparse vectors.
   *
   * @param dst_row_accessor write accessor
   * @param c column number within the column block
   * @param n_items number of items to read starting from @p c
   * @param phi FEEValuation class from which the data will be taken
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  void distribute_local_to_global(
    BlockCSRMatrixIterators::RowsBlockAccessor<Number, false> &dst_row_accessor,
    const unsigned int c,
    const unsigned int n_items,
    const FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi);

  /**
   * FIXME: should be moved upstream
   */
  class SparsityPatternStandard : public SparsityPatternBase
  {
  public:
    using size_type = SparsityPatternBase::size_type;

    static const size_type invalid_entry;

    using SparsityPatternBase::reinit;

    SparsityPatternStandard();

    virtual void
    reinit(const size_type m,
           const size_type n,
           const ArrayView<const unsigned int> &row_lengths) override;

    size_type operator()(const size_type i, const size_type j) const;

    void copy_from(const DynamicSparsityPattern &dsp);

    template <typename number>
    friend class BlockCSRMatrix;
  };

  template <typename NumberType>
  class BlockCSRMatrix : public virtual Subscriptor
  {
  public:

    /**
     * A divisor for column block sizes used when allocating
     * memory for dense blocks stored in row-major format.
     */
    static const unsigned int col_block_size_factor;

    /**
     * Declare type for container size.
     */
    using size_type = types::global_dof_index;

    /**
     * Type of the matrix entries. This alias is analogous to
     * <tt>value_type</tt> in the standard library containers.
     */
    using value_type = NumberType;

    /**
     * Declare a type that has holds real-valued numbers with the same precision
     * as the template argument to this class. If the template argument of this
     * class is a real data type, then real_type equals the template argument.
     * If the template argument is a std::complex type then real_type equals the
     * type underlying the complex numbers.
     *
     * This alias is used to represent the return type of norms.
     */
    using real_type =
      typename numbers::NumberTraits<value_type>::real_type;

    /**
     * Typedef of an iterator class walking over all the nonzero blocks of this
     * matrix. This iterator cannot change the values of the matrix.
     */
    using const_iterator = BlockCSRMatrixIterators::Iterator<NumberType, true>;

    /**
     * Typedef of an iterator class walking over all the nonzero entries of this
     * matrix. This iterator @em can change the values in the blocks, but of
     * course can't change the sparsity pattern as this is fixed once a block sparse
     * matrix is attached to it.
     */
    using iterator = BlockCSRMatrixIterators::Iterator<NumberType, false>;

    /**
     * @name Constructors and initialization
     */
    //@{
    /**
     * Constructor; initializes the matrix to be empty, without any structure,
     * i.e.  the matrix is not usable at all. This constructor is therefore only
     * useful for matrices which are members of a class. All other matrices
     * should be created at a point in the data flow where all necessary
     * information is available.
     *
     * You have to initialize the matrix before usage with reinit(const
     * SparsityPattern&).
     */
    BlockCSRMatrix();

    /**
     * Destructor. Free all memory, but do not release the memory of the
     * sparsity structure.
     */
    virtual ~BlockCSRMatrix() override;

    /**
     * Reinitialize the sparse matrix with the given sparsity pattern. The
     * latter tells the matrix how many nonzero elements there need to be
     * reserved for locally owned rows.
     *
     * The elements of the matrix are set to zero by this function.
     *
     * @parameter row_blocks defines row blocks of this matrix, its size should be
     * equal to 'sparsity.m()'
     * @parameter col_blocks defines column blocks of this matrix, its size should be
     * equal to 'sparsity.n()`.
     * @parameter symmetric if `true`, the matrix is considered to be
     * symmetric/hermitian. In this case row blocks and column blocks should match.
     * @parameter row_partitioner specifies 1D partitioner. This can be both row block partitioner and row partitioner.
     */
    void reinit(const DynamicSparsityPattern &sparsity,
                const std::shared_ptr<const BlockIndices> &row_blocks,
                const std::shared_ptr<const BlockIndices> &col_blocks,
                const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> row_partitioner,
                const bool symmetric = false);

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor. It also forgets the sparsity pattern it was
     * previously tied to.
     */
    void clear();

    //@}
    /**
     * @name Modifying entries
     */
    //@{

    /**
     * Read access to the element of the matrix at row @p i and column @p j.
     * This function throws an exception if the required element does not exist
     * in the matrix.
     */
    const value_type &operator()(const size_type i, const size_type j) const;

    /**
     * Same as above, but will return zero if the provided @p i and @p j index
     * is outside of the sparsity.
     */
    value_type el(const size_type i, const size_type j) const;

    /**
     * Same as above, but for locally owned row index @p i.
     */
    value_type local_el(const unsigned int i, const size_type j) const;

    /**
     * In contrast to the one above, this function allows modifying the value.
     */
    value_type &operator()(const size_type i, const size_type j);

    /**
     * This operator assigns a scalar to a matrix. To avoid confusion with
     * constructors, zero (when cast to the number type) is the only value
     * allowed for d.
     */
    void operator=(const value_type d);

    /**
     * Subtract the matrix @p B from the present one.
     * This is only allowed if both matrices have the same sparsity.
     */
    BlockCSRMatrix<NumberType> &
    operator-=(const BlockCSRMatrix<NumberType> &B);

    //@}
    /**
     * @name Matrix multiplication
     */
    //@{

    /**
     * Perform the matrix-matrix multiplication <tt>C = A * B</tt>, or, if an
     * optional vector argument is given, <tt>C = A * diag(V) * B</tt>, where
     * <tt>diag(V)</tt> defines a diagonal matrix with the vector entries.
     *
     * This function assumes that the calling matrix @p A and the argument @p B
     * have compatible sizes. By default, the output matrix @p C will be
     * resized appropriately.
     *
     * By default, i.e., if the optional argument @p rebuild_sparsity_pattern
     * is @p true, the sparsity pattern of the matrix C will be
     * changed to ensure that all entries that result from the product $AB$
     * can be stored in $C$. This is an expensive operation, and if there is
     * a way to predict the sparsity pattern up front, you should probably
     * build it yourself before calling this function with @p false as last
     * argument. In this case, the rebuilding of the sparsity pattern is
     * bypassed.
     *
     * When setting @p rebuild_sparsity_pattern to @p true (i.e., leaving it
     * at the default value), it is important to realize that the matrix
     * @p C passed as first argument still has to be initialized with a
     * sparsity pattern (either at the time of creation of the SparseMatrix
     * object, or via the SparseMatrix::reinit() function). This is because
     * we could create a sparsity pattern inside the current function, and
     * then associate @p C with it, but there would be no way to transfer
     * ownership of this sparsity pattern to anyone once the current function
     * finishes. Consequently, the function requires that @p C be already
     * associated with a sparsity pattern object, and this object is then
     * reset to fit the product of @p A and @p B.
     *
     * As a consequence of this, however, it is also important to realize
     * that the sparsity pattern of @p C is modified and that this would
     * render invalid <i>all other SparseMatrix objects</i> that happen
     * to <i>also</i> use that sparsity pattern object.
     */
    void mmult(BlockCSRMatrix<NumberType> &C,
               const BlockCSRMatrix<NumberType> &B,
               const bool rebuild_sparsity_pattern = true) const;

    /**
     * Perform the matrix-matrix multiplication with the transpose of
     * <tt>this</tt>, i.e., <tt>C = A<sup>T</sup> * B</tt>, or, if an optional
     * vector argument is given, <tt>C = A<sup>T</sup> * diag(V) * B</tt>, where
     * <tt>diag(V)</tt> defines a diagonal matrix with the vector entries.
     *
     * This function assumes that the calling matrix <tt>A</tt> and <tt>B</tt>
     * have compatible sizes. The size of <tt>C</tt> will be set within this
     * function.
     *
     * The content as well as the sparsity pattern of the matrix C will be
     * changed by this function, so make sure that the sparsity pattern is not
     * used somewhere else in your program. This is an expensive operation, so
     * think twice before you use this function.
     *
     * There is an optional flag <tt>rebuild_sparsity_pattern</tt> that can be
     * used to bypass the creation of a new sparsity pattern and instead uses
     * the sparsity pattern stored in <tt>C</tt>. In that case, make sure that
     * it really fits. The default is to rebuild the sparsity pattern.
     *
     * @note Rebuilding the sparsity pattern requires changing it. This means
     * that all other matrices that are associated with this sparsity pattern
     * will then have invalid entries.
     */
    void Tmmult(BlockCSRMatrix<NumberType> &C,
                const BlockCSRMatrix<NumberType> &B,
                const bool rebuild_sparsity_pattern = true) const;

    /**
     * Evaluates trace of Tmmult operation.
     * To that end, of course, only diagonal elements have to be calculated.
     * This operation makes sense only if the resulting matrix is square,
     * i.e. the number of columns in @p B and in this object match.
     */
    value_type Tr_Tmmult(const BlockCSRMatrix<NumberType> &B) const;

    //@}
    /**
     * @name Iterators
     */
    //@{

    /**
     * Return an iterator pointing to the first element of the matrix.
     *
     * Note the discussion in the general documentation of this class about the
     * order in which elements are accessed.
     */
    const_iterator begin() const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator begin();

    /**
     * Return an iterator pointing the element past the last one of this matrix.
     */
    const_iterator end() const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator end();

    /**
     * Return an iterator pointing to the first element of block row @p r.
     *
     * Note that if the given row is empty, i.e. does not contain any nonzero
     * entries, then the iterator returned by this function equals
     * <tt>end(r)</tt>. The returned iterator may not be dereferenceable in that
     * case if neither row @p r nor any of the following rows contain any
     * nonzero entries.
     */
    const_iterator begin(const size_type r) const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator begin(const size_type r);

    /**
     * Return an iterator pointing the element past the last one of block row @p
     * r , or past the end of the entire sparsity pattern if none of the rows
     * after
     * @p r contain any entries at all.
     *
     * Note that the end iterator is not necessarily dereferenceable. This is in
     * particular the case if it is the end iterator for the last block row of a
     * matrix.
     */
    const_iterator end(const size_type r) const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator end(const size_type r);

    /**
     * Same as above but for local row @p r.
     */
    const_iterator begin_local(const unsigned int r) const;

    /**
     * Same as above but for local row @p r.
     */
    iterator begin_local(const unsigned int r);

    /**
     * Same as above but for local row @p r
     */
    const_iterator end_local(const unsigned int r) const;

    /**
     * Same as above but for local row @p r
     */
    iterator end_local(const unsigned int r);

    //@}
    /**
     * @name Input/Output
     */
    //@{

    /**
     * Copy the contents of this matrix into @p matrix .
     *
     * @note This function should only be used for relatively small matrix
     * dimensions. It is primarily intended for debugging purposes.
     */
    template <typename DenseMatrix>
    void copy_to(DenseMatrix &matrix) const;

    //@}

    /**
     * @name Parallel data exchange
     */
    //@{

    /**
     * This function copies the data that has accumulated in the data buffer
     * for ghost indices to the owning processor. For the meaning of the
     * argument @p operation, see the entry on
     * @ref GlossCompress "Compressing distributed vectors and matrices"
     * in the glossary.
     */
    void
    compress(::dealii::VectorOperation::values operation);

    /**
     * Fills ghost rows with data stored on respective owing processes.
     * This function is needed
     * before reading from ghosts. The function is @p const even though
     * ghost data is changed. This is needed to allow functions with a @p
     * const vector to perform the data exchange without creating
     * temporaries.
     *
     * After calling this method, write access to ghost elements of the
     * vector is forbidden and an exception is thrown. Only read access to
     * ghost elements is allowed in this state.
     */
    void update_ghost_values() const;

    /**
     * Initiates communication for the @p compress() function with non-
     * blocking communication. This function does not wait for the transfer
     * to finish, in order to allow for other computations during the time
     * it takes until all data arrives.
     */
    void
      compress_start(
        const unsigned int                communication_channel = 0,
        ::dealii::VectorOperation::values operation = VectorOperation::add);

    /**
     * For all requests that have been initiated in compress_start, wait for
     * the communication to finish. Once it is finished, add or set the data
     * (depending on the flag operation) to the respective positions in the
     * owning processor, and clear the contents in the ghost data fields.
     * The meaning of this argument is the same as in compress().
     */
    void
    compress_finish(::dealii::VectorOperation::values operation);

    /**
     * Initiates communication for the @p update_ghost_values() function
     * with non-blocking communication. This function does not wait for the
     * transfer to finish, in order to allow for other computations during
     * the time it takes until all data arrives.
     *
     * Before the data is actually exchanged, the function must be followed
     * by a call to @p update_ghost_values_finish().
     *
     * In case this function is called for more than one matrix before @p
     * update_ghost_values_finish() is invoked, it is mandatory to specify a
     * unique communication channel to each such call, in order to avoid
     * several messages with the same ID that will corrupt this operation.
     */
    void
    update_ghost_values_start(
      const unsigned int communication_channel = 0) const;

    /**
     * For all requests that have been started in update_ghost_values_start,
     * wait for the communication to finish.
     *
     * Must follow a call to the @p update_ghost_values_start function
     * before reading data from ghost indices.
     */
    void
    update_ghost_values_finish() const;

    /**
     * This method zeros the entries on ghost rows, but does not touch
     * locally owned rows.
     *
     * After calling this method, read access to ghost elements of the
     * matrix is forbidden and an exception is thrown. Only write access to
     * ghost elements is allowed in this state.
     */
    void
    zero_out_ghosts() const;

    /**
     * Return whether the matrix currently is in a state where ghost values
     * can be read or not. This is the same functionality as parallel
     * vectors have. If this method returns false, this only means that
     * read-access to ghost elements is prohibited whereas write access is
     * still possible (to those entries specified as ghosts during
     * initialization), not that there are no ghost elements at all.
     *
     * @see
     * @ref GlossGhostedVector "vectors with ghost elements"
     */
    bool
    has_ghost_elements() const;

    /**
     * Returns `true` if the matrix was initialized with MPI partitioner for
     * blocks
     */
    bool is_block_partitioned() const;

    /**
     * Return partitioner of this matrix.
     */
    const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> &
    get_partitioner() const;

    //@}

    /**
     * @name Miscellaneous
     */
    //@{

    /**
     * Output of the matrix in user-defined format given by the specified
     * precision and width. This function saves width and precision of the
     * stream before setting these given values for output, and restores the
     * previous values after output.
     */
    void print(std::ostream &s,
               const unsigned int width = 5,
               const unsigned int precision = 2) const;

    /**
     * Return the memory consumption of this class in bytes.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Return true if matrix is symmetric
     */
    bool is_symmetric() const;

    /**
     * Return number of rows represented by this matrix.
     *
     * @note this shall not be confused with the number of row blocks.
     */
    size_type m() const;

    /**
     * Return number of columns represented by this matrix.
     *
     * @note this shall not be confused with the number of column blocks.
     */
    size_type n() const;

    /**
     * Return number of locally owned row blocks.
     */
    unsigned int n_local_row_blocks() const;

    /**
     * Return range of locally owned rows.
     */
    std::pair<size_type, size_type> local_range() const;

    /**
     * Return a (constant) reference to the underlying sparsity pattern of this
     * matrix.
     */
    const SparsityPatternStandard &get_sparsity_pattern() const;

    /**
     * Return row blocks used in this matrix
     */
    const std::shared_ptr<const BlockIndices> get_row_blocks() const;

    /**
     * Return columns blocks used in this matrix
     */
    const std::shared_ptr<const BlockIndices> &get_col_blocks() const;

    /**
     * Return starting index of this block row @p r and ranges `[begin,end)`
     * relative to first element in the pair.
     */
    std::pair<types::global_dof_index,
              std::vector<std::pair<unsigned int, unsigned int>>>
    get_block_data(const unsigned int r) const;

    /**
     * Compute the Frobenius norm of the matrix. Return value is the root of the
     * square sum of all matrix entries.
     */
    NumberType frobenius_norm() const;

    //@}

    /**
     * Exception for accessing element outside of sparsity parttern
     */
    DeclException2(
      ExcInvalidIndex,
      int,
      int,
      << "You are trying to access the matrix entry with index <" << arg1 << ','
      << arg2
      << ">, but this entry does not exist in the sparsity pattern "
         "of this matrix."
         "\n\n"
         "The most common cause for this problem is that you used "
         "a method to build the sparsity pattern that did not "
         "(completely) take into account all of the entries you "
         "will later try to write into. An example would be "
         "building a sparsity pattern that does not include "
         "the entries you will write into due to constraints "
         "on degrees of freedom such as hanging nodes or periodic "
         "boundary conditions. In such cases, building the "
         "sparsity pattern will succeed, but you will get errors "
         "such as the current one at one point or other when "
         "trying to write into the entries of the matrix.");

    /**
     * Exception for reading a ghost element.
     */
    DeclException1(
      ExcReadGhost,
      unsigned int,
      << "You are trying to read a ghost block " << arg1
      << "of this matrix, but it has not imported its ghost values yet. "
      << "Call update_ghost_values() to get read access.");

    /**
     * Exception for writing into a ghost element.
     */
    DeclException1(
      ExcWriteGhost,
      unsigned int,
      << "You are trying to get write access into ghost block " << arg1
      << "of this matrix, but it has already imported its ghost values. "
      << "Call zero_out_ghosts() prior to writing into ghost blocks.");

    /**
     * Return row/column major linear index for the element (i,j) within the block (M,N).
     */
    static
    unsigned int local_index(const unsigned int i, const unsigned int j, const unsigned int M, const unsigned int N);

  private:

    /**
     * Internal function to setup sparsity for ghost rows or row blocks
     */
    void setup_ghosts(const DynamicSparsityPattern &locally_owned_sparsity);

    /**
     * Auxiliary function to unify read/write from/to imported data during
     * parallel data exchange.
     */
    template <bool Constness, bool Addition = false>
    static void read_write_import(
      typename std::conditional<Constness,
                                const BlockCSRMatrix<NumberType> &,
                                BlockCSRMatrix<NumberType> &>::type matrix);

    /**
     * Sparsity pattern used for this matrix.
     */
    SparsityPatternStandard sp;

    /**
     * Values of all nonzero entries.
     */
    AlignedVector<NumberType> values;

    /**
     * For each linear index points to the start of the block in `values`.
     * Its size is the number of non-zero blocks plus one.
     * `data_start.back()` is the used size of `values` array.
     */
    std::vector<size_type> data_start;

    /**
     * Auxiliary vector used in Tr_Tmmult
     */
    mutable AlignedVector<NumberType> diagonal;

    /**
     * Auxiliary vector used in Tr_Tmmult
     */
    mutable std::vector<size_type> diagonal_data_start;

    /**
     * a flag for symmetric/hermitian matrices
     */
    bool symmetric;

    /**
     * Total number of row blocks across MPI communicator
     */
    unsigned int n_row_blocks;

    /**
     * Number of locally owned row blocks on this MPI process.
     */
    unsigned int n_owned_row_blocks;

    /**
     * Number of locally owned rows (same as local_size() for
     * LinearAlgebra::distributed::Vector)
     */
    unsigned int n_owned_rows;

    /**
     * Start of the locally owned range of rows
     */
    types::global_dof_index owned_row_start;

    /**
     * Total number of rows across MPI communicator
     */
    types::global_dof_index n_rows;

    /**
     * A shared pointer to the row block indices of this matrix.
     */
    std::shared_ptr<BlockIndices> row_blocks;

    /**
     * A shared pointer to the column block indices of this matrix.
     */
    std::shared_ptr<const BlockIndices> col_blocks;

    /**
     * MPI partitioner for rows of this matrix. This can be both
     * block partitioner or row partitioner.
     */
    std::shared_ptr<const dealii::Utilities::MPI::Partitioner> partitioner;

    /**
     * A flag to denote if the row partitioner corresponds to blocks or to the
     * total number of rows.
     */
    bool block_partitioner;

    /**
     * Stores whether the matrix currently allows for reading ghost elements
     * or not. Note that this is to ensure consistent ghost data and does
     * not indicate whether the matrix actually can store ghost elements.
     */
    mutable bool matrix_is_ghosted;

    /**
     * Similar to second element in MPI::Partitioner::ghost_targets() contains
     * information on how many elements we expect from each process during
     * update_ghosts() call.
     */
    std::vector<unsigned int> ghost_targets;

    /**
     * Similar to second element in MPI::Partitioner::import_targets() contains
     * information on how many elements we import during compress() from each
     * process.
     */
    std::vector<unsigned int> import_targets;

    /**
     * Sum over all elements in `import_targets`, i.e. total size of
     * import_data.
     */
    unsigned int n_import_data;

    /**
     * Temporary storage that holds the data that is sent to this processor
     * in @p compress() or sent from this processor in @p
     * update_ghost_values.
     */
    mutable std::unique_ptr<NumberType[]> import_data;

    /**
     * A reworked import_indices() from Utilities::MPI::Partitioner.
     * Namely, we store a vector of `[begin,end)` ranges of import indices
     * separately for each block row they belong to.
     *
     * Note that import_indices() are NOT unique, i.e. Utilities::MPI::Partitioner
     * may have the same `[begin,end)` range for different import targets.
     */
    std::vector<std::pair<unsigned int,
                          std::vector<std::pair<unsigned int, unsigned int>>>>
      import_indices;

    /**
     * Similar to import_indices, for each ghost block, store the start
     * of the block and its ranges within this block.
     */
    std::vector<std::pair<types::global_dof_index,
                          std::vector<std::pair<unsigned int, unsigned int>>>>
      ghost_data;

#ifdef DEAL_II_WITH_MPI
    /**
     * A vector that collects all requests from @p compress() operations.
     * This class uses persistent MPI communicators, i.e., the communication
     * channels are stored during successive calls to a given function. This
     * reduces the overhead involved with setting up the MPI machinery, but
     * it does not remove the need for a receive operation to be posted
     * before the data can actually be sent.
     */
    std::vector<MPI_Request>   compress_requests;

    /**
     * A vector that collects all requests from @p update_ghost_values()
     * operations. This class uses persistent MPI communicators.
     */
    mutable std::vector<MPI_Request>   update_ghost_values_requests;
#endif

    /**
     * For parallel loops with TBB, this member variable stores the affinity
     * information of loops.
     */
    mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
      thread_loop_partitioner;

    /**
     * A lock that makes sure that the @p compress and @p
     * update_ghost_values functions give reasonable results also when used
     * with several threads.
     */
    mutable Threads::Mutex mutex;

    /**
     * A helper function that clears the compress_requests and
     * update_ghost_values_requests field. Used in reinit functions.
     */
    void clear_mpi_requests ();

    /**
     * Give access to internal details to the iterator/accessor classes.
     */
    template <typename, bool>
    friend class BlockCSRMatrixIterators::Iterator;
    template <typename, bool>
    friend class BlockCSRMatrixIterators::Accessor;
    template <typename, bool>
    friend class BlockCSRMatrixIterators::RowsAccessor;
    template <typename, bool>
    friend class BlockCSRMatrixIterators::RowsAccessorBase;
    template <typename, bool>
    friend class BlockCSRMatrixIterators::RowsBlockAccessor;
  };

  // -------------------  inline and template functions ----------------
#ifndef DOXYGEN

  // We need to have a separate declaration for static const members
  template <typename NumberType>
  const unsigned int
    BlockCSRMatrix<NumberType>::col_block_size_factor = 64 / sizeof(NumberType);

  namespace internal
  {
    /**
     * Computes the smallest integer value divisible by @factor not less than @p
     * number.
     */
    inline unsigned int ceil_divisible_by(const unsigned int number,
                                          const unsigned int factor)
    {
      return factor * (number / factor + (number % factor > 0));
    }

    /**
     * Given @p size, calculate lowest integers such that
     * `size * sizeof(NumberType)` is 64 byte aligned.
     * This also means that `size * sizeof(NumberType)`
     * is divisable by SIMD registers size.
     */
    template <typename NumberType>
    inline unsigned int padded_size(const unsigned int size)
    {
      Assert(64 % sizeof(NumberType) == 0, ExcInternalError());
      Assert(
        ceil_divisible_by(size,
                          VectorizedArray<NumberType>::n_array_elements) <=
          ceil_divisible_by(size, 64 / sizeof(NumberType)),
        ExcMessage("We assume that SIMD width is less or equal to 64 bytes."));

      return ceil_divisible_by(
        size, BlockCSRMatrix<NumberType>::col_block_size_factor);
    }
  } // namespace internal



  template <typename NumberType>
  inline unsigned int
  BlockCSRMatrix<NumberType>::local_index(const unsigned int i, const unsigned int j, const unsigned int , const unsigned int N)
  {
    // row-major is better for matrix-free operators
    return i*internal::padded_size<NumberType>(N)+j;
  }



  template <typename NumberType>
  inline std::pair<typename BlockCSRMatrix<NumberType>::size_type,
                   typename BlockCSRMatrix<NumberType>::size_type>
  BlockCSRMatrix<NumberType>::local_range() const
  {
    return {owned_row_start, owned_row_start + n_owned_rows};
  }



  template <typename NumberType>
  inline const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> &
  BlockCSRMatrix<NumberType>::get_partitioner() const
  {
    return partitioner;
  }



  template <typename NumberType>
  inline unsigned int
  BlockCSRMatrix<NumberType>::n_local_row_blocks() const
  {
    return n_owned_row_blocks;
  }



  template <typename NumberType>
  inline const std::shared_ptr<const BlockIndices>
  BlockCSRMatrix<NumberType>::get_row_blocks() const
  {
    return row_blocks;
  }



  template <typename NumberType>
  inline const std::shared_ptr<const BlockIndices> &
  BlockCSRMatrix<NumberType>::get_col_blocks() const
  {
    return col_blocks;
  }



  template <typename NumberType>
  inline const SparsityPatternStandard &
  BlockCSRMatrix<NumberType>::get_sparsity_pattern() const
  {
    return sp;
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::size_type
  BlockCSRMatrix<NumberType>::m() const
  {
    return n_rows;
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::size_type
  BlockCSRMatrix<NumberType>::n() const
  {
    return col_blocks->total_size();
  }



  template <typename NumberType>
  inline bool BlockCSRMatrix<NumberType>::is_symmetric() const
  {
    return symmetric;
  }



  template <typename NumberType>
  inline bool
  BlockCSRMatrix<NumberType>::has_ghost_elements() const
  {
    return matrix_is_ghosted;
  }



  template <typename NumberType>
  inline bool
  BlockCSRMatrix<NumberType>::is_block_partitioned() const
  {
    return block_partitioner;
  }



  template <typename NumberType>
  inline
  std::pair<types::global_dof_index,
            std::vector<std::pair<unsigned int, unsigned int>>>
  BlockCSRMatrix<NumberType>::get_block_data(const unsigned int r) const
  {
    Assert(r < row_blocks->size(), ExcMessage("AssertRange"));
    if (r < n_owned_row_blocks)
      return {owned_row_start + row_blocks->block_start(r), {{0, row_blocks->block_size(r)}}};
    else
      return ghost_data[r - n_owned_row_blocks];
  }



  template <typename NumberType>
  inline
  void BlockCSRMatrix<NumberType>::
  operator=(const typename BlockCSRMatrix<NumberType>::value_type d)
  {
    (void)d;
    Assert(d == NumberType(0), ExcScalarAssignmentOnlyForZeroValue());

    values.fill();
  }



  template <typename NumberType>
  template <typename DenseMatrix>
  void BlockCSRMatrix<NumberType>::copy_to(DenseMatrix &matrix) const
  {
    matrix.reinit(this->m(), this->n());
    matrix = NumberType();

    if (sp.n_nonzero_elements() == 0)
      return;

    // go through all non-empty blocks and copy content inefficiently
    for (types::global_dof_index i = 0; i < n_owned_row_blocks; ++i)
      {
        const auto row_start = owned_row_start + row_blocks->block_start(i);
        const auto &row_size = row_blocks->block_size(i);
        const auto end = this->end_local(i);
        for (auto it = this->begin_local(i); it != end; ++it)
          {
            const auto &col_start = col_blocks->block_start(it->column());
            const auto &col_size = col_blocks->block_size(it->column());
            for (unsigned int ii = 0; ii < row_size; ++ii)
              for (unsigned int jj = 0; jj < col_size; ++jj)
                matrix(row_start + ii, col_start + jj) =
                  *(it->data() + BlockCSRMatrix<NumberType>::local_index(
                                   ii, jj, row_size, col_size));
          }
      }

    // do block reduction on serial matrix
    dealii::Utilities::MPI::sum(
      matrix, partitioner->get_mpi_communicator(), matrix);
  }



  template <typename NumberType>
  inline void BlockCSRMatrix<NumberType>::clear()
  {
    sp.reinit(0,0,0);
    values.clear();
    row_blocks.reset();
    col_blocks.reset();
  }


  //---------------------------- ITERATORS ---------------------------------

  namespace BlockCSRMatrixIterators
  {

    //--------------------------
    //-------- ACCESSOR --------
    //--------------------------

    template <typename NumberType, bool Constness>
    inline Accessor<NumberType, Constness>::Accessor(
      MatrixType *matrix, const std::size_t linear_index)
      : SparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern(),
                                           linear_index),
        matrix(matrix),
        current_data(matrix->values.data() + matrix->data_start[linear_index])
    {
    }



    template <typename NumberType, bool Constness>
    inline Accessor<NumberType, Constness>::Accessor(MatrixType *matrix)
      : SparsityPatternIterators::Accessor(
          &matrix->get_sparsity_pattern()),
        matrix(matrix),
        current_data(nullptr)
    {
    }



    template <typename NumberType, bool Constness>
    inline Accessor<NumberType, Constness>::Accessor()
      : SparsityPatternIterators::Accessor(),
        matrix(nullptr),
        current_data(nullptr)
    {
    }



    template <typename NumberType, bool Constness>
    inline Accessor<NumberType, Constness>::Accessor(
      const BlockCSRMatrixIterators::Accessor<NumberType, false> &a)
      : SparsityPatternIterators::Accessor(a),
        matrix(&a.get_matrix()),
        current_data(a.data())
    {
    }



    template <typename NumberType, bool Constness>
    inline
    typename Accessor<NumberType, Constness>::pointer
    Accessor<NumberType, Constness>::data() const
    {
      return current_data;
    }



    template <typename NumberType, bool Constness>
    inline typename Accessor<NumberType, Constness>::MatrixType &
    Accessor<NumberType, Constness>::get_matrix() const
    {
      Assert(matrix != nullptr, ExcInternalError());
      return *matrix;
    }



    template <typename NumberType, bool Constness>
    inline void Accessor<NumberType, Constness>::advance()
    {
      // first move data
      Assert (this->global_index() + 1 < matrix->data_start.size(), ExcMessage("AssertRange"));
      current_data += matrix->data_start[this->global_index() + 1] -  matrix->data_start[this->global_index()];

      // then advance parent class and its linear index
      SparsityPatternIterators::Accessor::advance();
    }


    //--------------------------
    //------- ITERATOR ---------
    //--------------------------

    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness>::Iterator(MatrixType *matrix, const std::size_t linear_index)
      : accessor(matrix, linear_index)
    {
    }



    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness>::Iterator(MatrixType *matrix)
      : accessor(matrix)
    {
    }



    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness>::Iterator()
      : accessor()
    {
    }



    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness>::Iterator(
      const BlockCSRMatrixIterators::Iterator<NumberType, false> &i)
      : accessor(*i)
    {
    }



    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness> &Iterator<NumberType, Constness>::
    operator++()
    {
      accessor.advance();
      return *this;
    }



    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness> Iterator<NumberType, Constness>::
    operator++(int)
    {
      const Iterator iter = *this;
      accessor.advance();
      return iter;
    }



    template <typename NumberType, bool Constness>
    inline const Accessor<NumberType, Constness> &
      Iterator<NumberType, Constness>::operator*() const
    {
      return accessor;
    }



    template <typename NumberType, bool Constness>
    inline const Accessor<NumberType, Constness> *
      Iterator<NumberType, Constness>::operator->() const
    {
      return &accessor;
    }



    template <typename NumberType, bool Constness>
    inline bool Iterator<NumberType, Constness>::
    operator==(const Iterator &other) const
    {
      return (accessor == other.accessor);
    }



    template <typename NumberType, bool Constness>
    inline bool Iterator<NumberType, Constness>::
    operator!=(const Iterator &other) const
    {
      return !(*this == other);
    }



    template <typename NumberType, bool Constness>
    inline bool Iterator<NumberType, Constness>::
    operator<(const Iterator &other) const
    {
      Assert(&accessor.get_matrix() == &other.accessor.get_matrix(),
             ExcInternalError());

      return (accessor < other.accessor);
    }



    template <typename NumberType, bool Constness>
    inline bool Iterator<NumberType, Constness>::
    operator>(const Iterator &other) const
    {
      return (other < *this);
    }



    template <typename NumberType, bool Constness>
    inline Iterator<NumberType, Constness> Iterator<NumberType, Constness>::
    operator+(const size_type n) const
    {
      Iterator x = *this;
      for (size_type i = 0; i < n; ++i)
        ++x;

      return x;
    }

    //--------------------------------
    //----- RowsBlockAccessor --------
    //--------------------------------

    template <typename NumberType, bool Constness>
    inline unsigned int
    RowsBlockAccessor<NumberType, Constness>::n_row_blocks() const
    {
      return row_blocks.size();
    }



    template <typename NumberType, bool Constness>
    inline types::global_dof_index
    RowsBlockAccessor<NumberType, Constness>::get_current_block_column() const
    {
      return col_block;
    }



    template <typename NumberType, bool Constness>
    inline types::global_dof_index
    RowsBlockAccessor<NumberType, Constness>::get_col_block_size() const
    {
      return col_block_size;
    }



    template <typename NumberType, bool Constness>
    inline void
    RowsBlockAccessor<NumberType, Constness>::process_active_rows_vectorized(
      const std::function<
        void(const ArrayView<const std::pair<unsigned int, unsigned int>> &,
             vectorized_pointer const block_start,
             const unsigned int stride)> &func) const
    {
      for (const auto &block : row_blocks)
        if (block.active)
          func(block.dof_view, block.pointer, stride);
    }

    //--------------------------
    //----- ROWACCESSOR --------
    //--------------------------

    template <typename NumberType, bool Constness>
    inline unsigned int RowsAccessorBase<NumberType, Constness>::size() const
    {
      return matrix->m();
    }



    template <typename NumberType, bool Constness>
    inline unsigned int RowsAccessorBase<NumberType, Constness>::n_row_blocks() const
    {
      return active_row_blocks.size();
    }



    template <typename NumberType, bool Constness>
    inline
    NumberType
    RowsAccessorBase<NumberType, Constness>::operator()(const types::global_dof_index index) const
    {
      Assert(this->matrix->block_partitioner == false, ExcNotImplemented());
      return local_element(this->matrix->partitioner->global_to_local(index));
    }


// #define BINARY_SEARCH

    template <typename NumberType, bool Constness>
    inline
    NumberType
    RowsAccessorBase<NumberType, Constness>::local_element(
      const unsigned int local_index) const
    {
      Assert (active_row_blocks.size() > 0, ExcMessage("Not initialized"));
#ifdef DEBUG
      const auto pair = matrix->row_blocks->global_to_local(local_index);
#endif

#ifdef BINARY_SEARCH
      // get iterator using binary search with custom comparator. wee need
      // last element with block_start <= local_index
      // see https://stackoverflow.com/a/50225224/888478
      const auto d_it = std::lower_bound(
        active_row_blocks.rbegin(),
        active_row_blocks.rend(),
        local_index,
        [](const DataType &d, const unsigned int &row) -> bool {
          return d.block_start > row;
        });

      const auto &d = *d_it;

      Assert(d_it != active_row_blocks.rend(),
             ExcMessage(
               "The RowsAccessor was not initialized to access row " +
               std::to_string(local_index) + " at block (" + std::to_string(pair.first) + ", " + std::to_string(pair.second) + ")"));
#else
      unsigned int ind = active_row_blocks.size()-1;
      while (active_row_blocks[ind].block_start > local_index)
        --ind;
      const auto &d = active_row_blocks[ind];
#endif

      // make sure what we found actually is the same as going through the global_to_local route:
      AssertDimension ( d.block_start, matrix->row_blocks->block_start(pair.first) );

      // if iterator is not at the end and it's column matches current
      // column, we can return the value
      if (d.active)
        {
          Assert(d.it != d.end && d.it->column() == block_col,
                 ExcInternalError());

          Assert(local_index - d.block_start < d.block_size,
                 ExcMessage(
                   "The RowsAccessor was not initialized to access row " +
                   std::to_string(local_index)));

          return *(d.data + BlockCSRMatrix<NumberType>::local_index(local_index - d.block_start,
                                                                    col_within_block,
                                                                    /*does not matter*/0,
                                                                    matrix->get_col_blocks()->block_size(block_col)));
        }
      else
        return 0.;
    }



    // same as local element, but with write access:
    template <typename NumberType>
    inline void RowsAccessor<NumberType, false>::set_local_element(
      const unsigned int local_index,
      const NumberType &val) const
    {
      modify_local_element<false>(local_index, val);
    }



    template <typename NumberType>
    inline void RowsAccessor<NumberType, false>::add_local_element(
      const unsigned int local_index,
      const NumberType &val) const
    {
      modify_local_element<true>(local_index, val);
    }



    template <typename NumberType>
    inline void
    RowsAccessor<NumberType, false>::add(const types::global_dof_index index,
                                         const NumberType &val) const
    {
      Assert(this->matrix->block_partitioner == false, ExcNotImplemented());
      modify_local_element<true>(
        this->matrix->partitioner->global_to_local(index), val);
    }



    template <typename NumberType>
    template <bool addition>
    inline void RowsAccessor<NumberType, false>::modify_local_element(
      const unsigned int local_index,
      const NumberType &val) const
    {
      Assert (this->active_row_blocks.size() > 0, ExcMessage("Not initialized"));

#ifdef DEBUG
      const auto pair = this->matrix->row_blocks->global_to_local(local_index);
#endif

#ifdef BINARY_SEARCH
      // get iterator using binary search with custom comparator. wee need
      // last element with block_start <= local_index
      // see https://stackoverflow.com/a/50225224/888478
      const auto d_it = std::lower_bound(
        this->active_row_blocks.rbegin(),
        this->active_row_blocks.rend(),
        local_index,
        [](const typename RowsAccessorBase<NumberType, false>::DataType &d,
           const unsigned int &row) -> bool {
          return d.block_start > row;
        });

      Assert(d_it != this->active_row_blocks.rend(),
             ExcMessage(
               "The RowsAccessor was not initialized to access row " +
               std::to_string(local_index) + " at block (" + std::to_string(pair.first) + ", " + std::to_string(pair.second) + ")"));

      const auto &d = *d_it;
#else
      unsigned int ind = this->active_row_blocks.size()-1;
      while (this->active_row_blocks[ind].block_start > local_index)
        --ind;
      const auto &d = this->active_row_blocks[ind];
#endif

      // make sure what we found actually is the same as going through the global_to_local route:
      AssertDimension ( d.block_start, this->matrix->row_blocks->block_start(pair.first) );

      // if iterator is not at the end and it's column matches current
      // column, we can return the value
      if (d.active)
        {
          Assert(d.it != d.end && d.it->column() == this->block_col,
                 ExcInternalError());

          Assert(local_index - d.block_start < d.block_size,
                 ExcMessage(
                   "The RowsAccessor was not initialized to access row " +
                   std::to_string(local_index)));

          if (addition)
            *(d.data + BlockCSRMatrix<NumberType>::local_index(local_index - d.block_start,
                                                               this->col_within_block,
                                                               /*does not matter*/0,
                                                               this->matrix->get_col_blocks()->block_size(this->block_col))) += val;
          else
            *(d.data + BlockCSRMatrix<NumberType>::local_index(local_index - d.block_start,
                                                               this->col_within_block,
                                                               /*does not matter*/0,
                                                               this->matrix->get_col_blocks()->block_size(this->block_col))) = val;
        }
      else
        Assert(val == NumberType(),
               typename BlockCSRMatrix<NumberType>::ExcInvalidIndex(
                 local_index, this->current_column()));
    }



    template <typename NumberType, bool Constness>
    inline
      typename RowsAccessorBase<NumberType, Constness>::MatrixType::size_type
      RowsAccessorBase<NumberType, Constness>::current_column() const
    {
      return matrix->col_blocks->local_to_global(block_col, col_within_block);
    }



    template <typename NumberType, bool Constness>
    inline
    void RowsAccessorBase<NumberType, Constness>::move_iterators()
    {
      // advance iterators if necessary to point
      // to the same block_col (if possible)
      for (auto &v : active_row_blocks)
        {
          for (;; ++v.it)
            {
              if (v.it == v.end)
                {
                  v.active = false;
                  break;
                }
              else if (v.it->column() == block_col)
                {
                  v.active = true;
                  break;
                }
              else if (v.it->column() > block_col)
                {
                  v.active = false;
                  break;
                }
            }

          if (v.active)
            {
              v.data = v.it->data();
            }
          else
            v.data = nullptr;
        }
    }



    template <typename NumberType, bool Constness>
    inline
    void RowsAccessorBase<NumberType, Constness>::advance(const unsigned int col)
    {
      const std::pair<unsigned int, unsigned int> col_pair = matrix->col_blocks->global_to_local(col);

      col_within_block = col_pair.second;
      Assert(col_pair.first >= block_col,
             ExcMessage(
               "Access only in ascending order of column blocks is possible: " +
               std::to_string(block_col) + " > " +
               std::to_string(col_pair.first)));

      // move iterators and adjust pointers
      block_col = col_pair.first;
      move_iterators();
    }



    template <typename NumberType, bool Constness>
    bool RowsAccessorBase<NumberType, Constness>::advance()
    {
      if (col_within_block + 1 < matrix->col_blocks->block_size(block_col))
        {
          // if we stay within the same block, just advance the local column
          // index
          ++col_within_block;

          move_iterators();

          return true;
        }
      else if (block_col == matrix->col_blocks->size() - 1)
        {
          // if we are at the last block, we are done
          return false;
        }
      else
        {
          // reset local column index and increment the block index
          col_within_block = 0;
          ++block_col;

          // move_iterators() will find iterators pointing to block_col
          // but we might not have any non-empty blocks for the
          // incremented block_col. Therefore first try to
          // increment iterators by one
          for (auto &v : active_row_blocks)
            if (v.it != v.end && v.it->column() < block_col)
              ++v.it;

          // and get the minimum column they point to
          types::global_dof_index current_min_column =
            std::numeric_limits<unsigned int>::max();
          for (auto &v : active_row_blocks)
            if (v.it != v.end)
              current_min_column = std::min(current_min_column, v.it->column());

          if (current_min_column < std::numeric_limits<unsigned int>::max())
            {
              // if we update block_col, it can't decrease
              Assert(current_min_column >= block_col, ExcInternalError());
              block_col = current_min_column;
            }

          // now we can set active flags and update pointer to data
          move_iterators();

          // check if there is any iterator that is not at the end
          for (auto &d : active_row_blocks)
            if (d.it != d.end)
              return true;

          // otherwise we are done (all iterators are at the end)
          return false;
        }
    }

  } // namespace BlockCSRMatrixIterators



  //--------------------------
  //---- BCSR begin/end ------
  //--------------------------



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::const_iterator
  BlockCSRMatrix<NumberType>::begin() const
  {
    return const_iterator(this, 0);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::const_iterator
  BlockCSRMatrix<NumberType>::end() const
  {
    return const_iterator(this);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::iterator
  BlockCSRMatrix<NumberType>::begin()
  {
    return iterator(this,0);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::iterator
  BlockCSRMatrix<NumberType>::end()
  {
    // FIXME: SparseMatrix does
    // iterator(this, cols->rowstart[cols->rows]);
    return iterator(this);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::const_iterator
  BlockCSRMatrix<NumberType>::begin(const size_type r) const
  {
    Assert(block_partitioner, ExcNotImplemented());
    return begin_local(partitioner->global_to_local(r));
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::const_iterator
  BlockCSRMatrix<NumberType>::end(const size_type r) const
  {
    Assert(block_partitioner, ExcNotImplemented());
    return end_local(partitioner->global_to_local(r));
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::iterator
  BlockCSRMatrix<NumberType>::begin(const size_type r)
  {
    Assert(block_partitioner, ExcNotImplemented());
    return begin_local(partitioner->global_to_local(r));
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::iterator
  BlockCSRMatrix<NumberType>::end(const size_type r)
  {
    Assert(block_partitioner, ExcNotImplemented());
    return end_local(partitioner->global_to_local(r));
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::const_iterator
  BlockCSRMatrix<NumberType>::begin_local(const unsigned int r) const
  {
    // do not allow reading a matrix which is not in ghost mode
    Assert(r < n_owned_row_blocks || matrix_is_ghosted == true,
           ExcReadGhost(r));
    Assert(r < sp.n_rows(), ExcIndexRange(r, 0, sp.n_rows()));
    return const_iterator(this, sp.rowstart[r]);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::const_iterator
  BlockCSRMatrix<NumberType>::end_local(const unsigned int r) const
  {
    // do not allow reading a matrix which is not in ghost mode
    Assert(r < n_owned_row_blocks || matrix_is_ghosted == true,
           ExcReadGhost(r));
    Assert(r < sp.n_rows(), ExcIndexRange(r, 0, sp.n_rows()));
    return const_iterator(this, sp.rowstart[r+1]);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::iterator
  BlockCSRMatrix<NumberType>::begin_local(const unsigned int r)
  {
    Assert(r < n_owned_row_blocks || matrix_is_ghosted == false,
           ExcWriteGhost(r));
    Assert(r < sp.n_rows(), ExcIndexRange(r, 0, sp.n_rows()));
    return iterator(this, sp.rowstart[r]);
  }



  template <typename NumberType>
  inline typename BlockCSRMatrix<NumberType>::iterator
  BlockCSRMatrix<NumberType>::end_local(const unsigned int r)
  {
    Assert(r < n_owned_row_blocks || matrix_is_ghosted == false,
           ExcWriteGhost(r));
    Assert(r < sp.n_rows(), ExcIndexRange(r, 0, sp.n_rows()));
    return iterator(this, sp.rowstart[r+1]);
  }



  //--------------------------
  //-- FEEvaluation alike ----
  //--------------------------


  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline
  void read_values(
    const BlockCSRMatrixIterators::RowsBlockAccessor<Number, true>
      &src_row_accessor,
    const unsigned int c,
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi,
    std::bitset<FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number>::
                  static_dofs_per_cell> &touched)
  {
    // make_vectorized_array is not constexpr
    static const VectorizedArray<Number> zero =
      make_vectorized_array<Number>(0);

    if (c == 0)
      {
        // read data and set touched
        src_row_accessor.process_active_rows_vectorized(
          [&](const ArrayView<const std::pair<unsigned int, unsigned int>>
                &dof_view,
              typename BlockCSRMatrixIterators::
                RowsBlockAccessor<Number, true>::vectorized_pointer const val,
              const unsigned int stride) {
            // loop over all rows within this block
            for (unsigned int i = 0; i < dof_view.size(); ++i)
              {
                // note: c == 0
                phi.submit_dof_value(val[dof_view[i].first * stride],
                                     dof_view[i].second);
                touched.set(dof_view[i].second);
              }
          });
      }
    else
      // just read data
      {
        src_row_accessor.process_active_rows_vectorized(
          [&](const ArrayView<const std::pair<unsigned int, unsigned int>>
                &dof_view,
              typename BlockCSRMatrixIterators::
                RowsBlockAccessor<Number, true>::vectorized_pointer const val,
              const unsigned int stride) {
            // loop over all rows within this block
            for (unsigned int i = 0; i < dof_view.size(); ++i)
              {
                phi.submit_dof_value(val[dof_view[i].first * stride + c],
                                     dof_view[i].second);
              }
          });
      }

    // flash all untouched
    for (unsigned int i = 0; i < touched.size(); ++i)
      if (!touched[i])
        phi.submit_dof_value(zero, i);
  }



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline void distribute_local_to_global(
    BlockCSRMatrixIterators::RowsBlockAccessor<Number, false> &dst_row_accessor,
    const unsigned int c,
    const FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi)
  {
    dst_row_accessor.process_active_rows_vectorized(
      [&](
        const ArrayView<const std::pair<unsigned int, unsigned int>> &dof_view,
        typename BlockCSRMatrixIterators::RowsBlockAccessor<Number, false>::
          vectorized_pointer const val,
        const unsigned int stride) {
        // loop over all rows within this block
        for (unsigned int i = 0; i < dof_view.size(); ++i)
          val[dof_view[i].first * stride + c] +=
            phi.get_dof_value(dof_view[i].second);
      });
  }

  //--------------------------
  //-- FEEvaluation fused ----
  //--------------------------

  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  void read_values_fused(
    const BlockCSRMatrixIterators::RowsBlockAccessor<Number, true>
      &src_row_accessor,
    const unsigned int c,
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi1,
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi2,
    std::bitset<FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number>::
                  static_dofs_per_cell> &touched)
  {
    static const VectorizedArray<Number> zero =
      make_vectorized_array<Number>(0);

    src_row_accessor.process_active_rows_vectorized(
      [&](
        const ArrayView<const std::pair<unsigned int, unsigned int>> &dof_view,
        typename BlockCSRMatrixIterators::RowsBlockAccessor<Number, true>::
          vectorized_pointer const val,
        const unsigned int stride) {
        // loop over all rows within this block
        for (unsigned int i = 0; i < dof_view.size(); ++i)
          {
            phi1.submit_dof_value(val[dof_view[i].first * stride + c],
                                  dof_view[i].second);
            phi2.submit_dof_value(val[dof_view[i].first * stride + c + 1],
                                  dof_view[i].second);
            touched.set(dof_view[i].second);
          }
      });

    // flash all untouched
    for (unsigned int i = 0; i < touched.size(); ++i)
      if (!touched[i])
        {
          phi1.submit_dof_value(zero, i);
          phi2.submit_dof_value(zero, i);
        }
  }



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  void distribute_local_to_global_fused(
    BlockCSRMatrixIterators::RowsBlockAccessor<Number, false> &dst_row_accessor,
    const unsigned int c,
    const FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi1,
    const FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> &phi2)
  {
    dst_row_accessor.process_active_rows_vectorized(
      [&](
        const ArrayView<const std::pair<unsigned int, unsigned int>> &dof_view,
        typename BlockCSRMatrixIterators::RowsBlockAccessor<Number, false>::
          vectorized_pointer const val,
        const unsigned int stride) {
        // loop over all rows within this block
        for (unsigned int i = 0; i < dof_view.size(); ++i)
          {
            val[dof_view[i].first * stride + c] +=
              phi1.get_dof_value(dof_view[i].second);
            val[dof_view[i].first * stride + c + 1] +=
              phi2.get_dof_value(dof_view[i].second);
          }
      });
  }

#endif // doxygen

  /**
   * Declare BlockCSRMatrixIterators::RowsAccessor as distributed vector.
   */
  template <typename Number, bool Constness>
  struct is_serial_vector<BlockCSRMatrixIterators::RowsAccessor<Number, Constness>>
    : std::false_type
  {};

  /**
   * Declare BlockCSRMatrixIterators::RowsAccessor as distributed vector.
   */
  template <typename Number>
  struct is_serial_vector<BlockCSRMatrix<Number>>
    : std::false_type
  {};

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_block_csr_matrix_h

