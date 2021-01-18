// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_trilinos_sparse_matrix_h
#  define dealii_trilinos_sparse_matrix_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_TRILINOS

#    include <deal.II/base/index_set.h>
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/full_matrix.h>
#    include <deal.II/lac/trilinos_epetra_vector.h>
#    include <deal.II/lac/trilinos_index_access.h>
#    include <deal.II/lac/trilinos_tpetra_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/vector_memory.h>
#    include <deal.II/lac/vector_operation.h>

#    include <Epetra_Comm.h>
#    include <Epetra_CrsGraph.h>
#    include <Epetra_Export.h>
#    include <Epetra_FECrsMatrix.h>
#    include <Epetra_Map.h>
#    include <Epetra_MultiVector.h>
#    include <Epetra_Operator.h>

#    include <cmath>
#    include <memory>
#    include <type_traits>
#    include <vector>
#    ifdef DEAL_II_WITH_MPI
#      include <Epetra_MpiComm.h>
#      include <mpi.h>
#    else
#      include <Epetra_SerialComm.h>
#    endif

DEAL_II_NAMESPACE_OPEN

// forward declarations
#    ifndef DOXYGEN
template <typename MatrixType>
class BlockMatrixBase;

template <typename number>
class SparseMatrix;
class SparsityPattern;
class DynamicSparsityPattern;

namespace TrilinosWrappers
{
  class SparseMatrix;
  class SparsityPattern;

  namespace SparseMatrixIterators
  {
    template <bool Constness>
    class Iterator;
  }
} // namespace TrilinosWrappers
#    endif

namespace TrilinosWrappers
{
  /**
   * Iterators for Trilinos matrices
   */
  namespace SparseMatrixIterators
  {
    /**
     * Exception
     */
    DeclException0(ExcBeyondEndOfMatrix);

    /**
     * Exception
     */
    DeclException3(ExcAccessToNonlocalRow,
                   std::size_t,
                   std::size_t,
                   std::size_t,
                   << "You tried to access row " << arg1
                   << " of a distributed sparsity pattern, "
                   << " but only rows " << arg2 << " through " << arg3
                   << " are stored locally and can be accessed.");

    /**
     * Handling of indices for both constant and non constant Accessor objects
     *
     * For a regular dealii::SparseMatrix, we would use an accessor for the
     * sparsity pattern. For Trilinos matrices, this does not seem so simple,
     * therefore, we write a little base class here.
     */
    class AccessorBase
    {
    public:
      /**
       * Declare the type for container size.
       */
      using size_type = dealii::types::global_dof_index;

      /**
       * Constructor.
       */
      AccessorBase(SparseMatrix *  matrix,
                   const size_type row,
                   const size_type index);

      /**
       * Row number of the element represented by this object.
       */
      size_type
      row() const;

      /**
       * Index in row of the element represented by this object.
       */
      size_type
      index() const;

      /**
       * Column number of the element represented by this object.
       */
      size_type
      column() const;

    protected:
      /**
       * Pointer to the matrix object. This object should be handled as a
       * const pointer or non-const by the appropriate derived classes. In
       * order to be able to implement both, it is not const here, so handle
       * with care!
       */
      mutable SparseMatrix *matrix;
      /**
       * Current row number.
       */
      size_type a_row;

      /**
       * Current index in row.
       */
      size_type a_index;

      /**
       * Discard the old row caches (they may still be used by other
       * accessors) and generate new ones for the row pointed to presently by
       * this accessor.
       */
      void
      visit_present_row();

      /**
       * Cache where we store the column indices of the present row. This is
       * necessary, since Trilinos makes access to the elements of its
       * matrices rather hard, and it is much more efficient to copy all
       * column entries of a row once when we enter it than repeatedly asking
       * Trilinos for individual ones. This also makes some sense since it is
       * likely that we will access them sequentially anyway.
       *
       * In order to make copying of iterators/accessor of acceptable
       * performance, we keep a shared pointer to these entries so that more
       * than one accessor can access this data if necessary.
       */
      std::shared_ptr<std::vector<size_type>> colnum_cache;

      /**
       * Cache for the values of this row.
       */
      std::shared_ptr<std::vector<TrilinosScalar>> value_cache;
    };

    /**
     * General template for sparse matrix accessors. The first template
     * argument denotes the underlying numeric type, the second the constness
     * of the matrix.
     *
     * The general template is not implemented, only the specializations for
     * the two possible values of the second template argument. Therefore, the
     * interface listed here only serves as a template provided since doxygen
     * does not link the specializations.
     */
    template <bool Constess>
    class Accessor : public AccessorBase
    {
      /**
       * Value of this matrix entry.
       */
      TrilinosScalar
      value() const;

      /**
       * Value of this matrix entry.
       */
      TrilinosScalar &
      value();
    };

    /**
     * The specialization for a const Accessor.
     */
    template <>
    class Accessor<true> : public AccessorBase
    {
    public:
      /**
       * Typedef for the type (including constness) of the matrix to be used
       * here.
       */
      using MatrixType = const SparseMatrix;

      /**
       * Constructor. Since we use accessors only for read access, a const
       * matrix pointer is sufficient.
       */
      Accessor(MatrixType *matrix, const size_type row, const size_type index);

      /**
       * Copy constructor to get from a const or non-const accessor to a const
       * accessor.
       */
      template <bool Other>
      Accessor(const Accessor<Other> &a);

      /**
       * Value of this matrix entry.
       */
      TrilinosScalar
      value() const;

    private:
      // Make iterator class a friend.
      template <bool>
      friend class Iterator;
    };

    /**
     * The specialization for a mutable Accessor.
     */
    template <>
    class Accessor<false> : public AccessorBase
    {
      class Reference
      {
      public:
        /**
         * Constructor.
         */
        Reference(const Accessor<false> &accessor);

        /**
         * Conversion operator to the data type of the matrix.
         */
        operator TrilinosScalar() const;

        /**
         * Set the element of the matrix we presently point to to @p n.
         */
        const Reference &
        operator=(const TrilinosScalar n) const;

        /**
         * Add @p n to the element of the matrix we presently point to.
         */
        const Reference &
        operator+=(const TrilinosScalar n) const;

        /**
         * Subtract @p n from the element of the matrix we presently point to.
         */
        const Reference &
        operator-=(const TrilinosScalar n) const;

        /**
         * Multiply the element of the matrix we presently point to by @p n.
         */
        const Reference &
        operator*=(const TrilinosScalar n) const;

        /**
         * Divide the element of the matrix we presently point to by @p n.
         */
        const Reference &
        operator/=(const TrilinosScalar n) const;

      private:
        /**
         * Pointer to the accessor that denotes which element we presently
         * point to.
         */
        Accessor &accessor;
      };

    public:
      /**
       * Typedef for the type (including constness) of the matrix to be used
       * here.
       */
      using MatrixType = SparseMatrix;

      /**
       * Constructor. Since we use accessors only for read access, a const
       * matrix pointer is sufficient.
       */
      Accessor(MatrixType *matrix, const size_type row, const size_type index);

      /**
       * Value of this matrix entry.
       */
      Reference
      value() const;

    private:
      // Make iterator class a friend.
      template <bool>
      friend class Iterator;

      // Make Reference object a friend.
      friend class Reference;
    };

    /**
     * This class acts as an iterator walking over the elements of Trilinos
     * matrices. The implementation of this class is similar to the one for
     * PETSc matrices.
     *
     * Note that Trilinos stores the elements within each row in ascending
     * order. This is opposed to the deal.II sparse matrix style where the
     * diagonal element (if it exists) is stored before all other values, and
     * the PETSc sparse matrices, where one can't guarantee a certain order of
     * the elements.
     *
     * @ingroup TrilinosWrappers
     */
    template <bool Constness>
    class Iterator
    {
    public:
      /**
       * Declare type for container size.
       */
      using size_type = dealii::types::global_dof_index;

      /**
       * Typedef for the matrix type (including constness) we are to operate
       * on.
       */
      using MatrixType = typename Accessor<Constness>::MatrixType;

      /**
       * Constructor. Create an iterator into the matrix @p matrix for the
       * given row and the index within it.
       */
      Iterator(MatrixType *matrix, const size_type row, const size_type index);

      /**
       * Copy constructor with optional change of constness.
       */
      template <bool Other>
      Iterator(const Iterator<Other> &other);

      /**
       * Prefix increment.
       */
      Iterator<Constness> &
      operator++();

      /**
       * Postfix increment.
       */
      Iterator<Constness>
      operator++(int);

      /**
       * Dereferencing operator.
       */
      const Accessor<Constness> &operator*() const;

      /**
       * Dereferencing operator.
       */
      const Accessor<Constness> *operator->() const;

      /**
       * Comparison. True, if both iterators point to the same matrix
       * position.
       */
      bool
      operator==(const Iterator<Constness> &) const;

      /**
       * Inverse of <tt>==</tt>.
       */
      bool
      operator!=(const Iterator<Constness> &) const;

      /**
       * Comparison operator. Result is true if either the first row number is
       * smaller or if the row numbers are equal and the first index is
       * smaller.
       */
      bool
      operator<(const Iterator<Constness> &) const;

      /**
       * Comparison operator. The opposite of the previous operator
       */
      bool
      operator>(const Iterator<Constness> &) const;

      /**
       * Exception
       */
      DeclException2(ExcInvalidIndexWithinRow,
                     size_type,
                     size_type,
                     << "Attempt to access element " << arg2 << " of row "
                     << arg1 << " which doesn't have that many elements.");

    private:
      /**
       * Store an object of the accessor class.
       */
      Accessor<Constness> accessor;

      template <bool Other>
      friend class Iterator;
    };

  } // namespace SparseMatrixIterators


  /**
   * This class implements a wrapper to use the Trilinos distributed sparse
   * matrix class Epetra_FECrsMatrix. This is precisely the kind of matrix we
   * deal with all the time - we most likely get it from some assembly
   * process, where also entries not locally owned might need to be written
   * and hence need to be forwarded to the owner process.  This class is
   * designed to be used in a distributed memory architecture with an MPI
   * compiler on the bottom, but works equally well also for serial processes.
   * The only requirement for this class to work is that Trilinos has been
   * installed with the same compiler as is used for generating deal.II.
   *
   * The interface of this class is modeled after the existing SparseMatrix
   * class in deal.II. It has almost the same member functions, and is often
   * exchangeable. However, since Trilinos only supports a single scalar type
   * (double), it is not templated, and only works with doubles.
   *
   * Note that Trilinos only guarantees that operations do what you expect if
   * the functions @p GlobalAssemble has been called after matrix assembly.
   * Therefore, you need to call SparseMatrix::compress() before you actually
   * use the matrix. This also calls @p FillComplete that compresses the
   * storage format for sparse matrices by discarding unused elements.
   * Trilinos allows to continue with assembling the matrix after calls to
   * these functions, though.
   *
   * <h3>Thread safety of Trilinos matrices</h3>
   *
   * When writing into Trilinos matrices from several threads in shared
   * memory, several things must be kept in mind as there is no built-in locks
   * in this class to prevent data races. Simultaneous access to the same
   * matrix row at the same time can lead to data races and must be explicitly
   * avoided by the user. However, it is possible to access <b>different</b>
   * rows of the matrix from several threads simultaneously under the
   * following three conditions:
   * <ul>
   * <li> The matrix uses only one MPI process.
   * <li> The matrix has been initialized with the reinit() method with a
   * DynamicSparsityPattern (that includes the set of locally relevant rows,
   * i.e., the rows that an assembly routine will possibly write into).
   * <li> The matrix has been initialized from a
   * TrilinosWrappers::SparsityPattern object that in turn has been
   * initialized with the reinit function specifying three index sets, one for
   * the rows, one for the columns and for the larger set of @p
   * writeable_rows, and the operation is an addition. At some point in the
   * future, Trilinos support might be complete enough such that initializing
   * from a TrilinosWrappers::SparsityPattern that has been filled by a
   * function similar to DoFTools::make_sparsity_pattern always results in a
   * matrix that allows several processes to write into the same matrix row.
   * However, Trilinos until version at least 11.12 does not correctly support
   * this feature.
   * </ul>
   *
   * Note that all other reinit methods and constructors of
   * TrilinosWrappers::SparsityPattern will result in a matrix that needs to
   * allocate off-processor entries on demand, which breaks thread-safety. Of
   * course, using the respective reinit method for the block Trilinos
   * sparsity pattern and block matrix also results in thread-safety.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Matrix1
   */
  class SparseMatrix : public Subscriptor
  {
  public:
    /**
     * Declare the type for container size.
     */
    using size_type = dealii::types::global_dof_index;

    /**
     * Exception
     */
    DeclException1(ExcAccessToNonlocalRow,
                   std::size_t,
                   << "You tried to access row " << arg1
                   << " of a non-contiguous locally owned row set."
                   << " The row " << arg1
                   << " is not stored locally and can't be accessed.");

    /**
     * A structure that describes some of the traits of this class in terms of
     * its run-time behavior. Some other classes (such as the block matrix
     * classes) that take one or other of the matrix classes as its template
     * parameters can tune their behavior based on the variables in this
     * class.
     */
    struct Traits
    {
      /**
       * It is safe to elide additions of zeros to individual elements of this
       * matrix.
       */
      static const bool zero_addition_can_be_elided = true;
    };

    /**
     * Declare an alias for the iterator class.
     */
    using iterator = SparseMatrixIterators::Iterator<false>;

    /**
     * Declare an alias for the const iterator class.
     */
    using const_iterator = SparseMatrixIterators::Iterator<true>;

    /**
     * Declare an alias in analogy to all the other container classes.
     */
    using value_type = TrilinosScalar;

    /**
     * @name Constructors and initialization.
     */
    //@{
    /**
     * Default constructor. Generates an empty (zero-size) matrix.
     */
    SparseMatrix();

    /**
     * Generate a matrix that is completely stored locally, having #m rows and
     * #n columns.
     *
     * The number of columns entries per row is specified as the maximum
     * number of entries argument.
     */
    SparseMatrix(const size_type    m,
                 const size_type    n,
                 const unsigned int n_max_entries_per_row);

    /**
     * Generate a matrix that is completely stored locally, having #m rows and
     * #n columns.
     *
     * The vector <tt>n_entries_per_row</tt> specifies the number of entries
     * in each row.
     */
    SparseMatrix(const size_type                  m,
                 const size_type                  n,
                 const std::vector<unsigned int> &n_entries_per_row);

    /**
     * Generate a matrix from a Trilinos sparsity pattern object.
     */
    SparseMatrix(const SparsityPattern &InputSparsityPattern);

    /**
     * Move constructor. Create a new sparse matrix by stealing the internal
     * data.
     */
    SparseMatrix(SparseMatrix &&other) noexcept;

    /**
     * Copy constructor is deleted.
     */
    SparseMatrix(const SparseMatrix &) = delete;

    /**
     * operator= is deleted.
     */
    SparseMatrix &
    operator=(const SparseMatrix &) = delete;

    /**
     * Destructor. Made virtual so that one can use pointers to this class.
     */
    virtual ~SparseMatrix() override = default;

    /**
     * This function initializes the Trilinos matrix with a deal.II sparsity
     * pattern, i.e. it makes the Trilinos Epetra matrix know the position of
     * nonzero entries according to the sparsity pattern. This function is
     * meant for use in serial programs, where there is no need to specify how
     * the matrix is going to be distributed among different processors. This
     * function works in %parallel, too, but it is recommended to manually
     * specify the %parallel partitioning of the matrix using an Epetra_Map.
     * When run in %parallel, it is currently necessary that each processor
     * holds the sparsity_pattern structure because each processor sets its
     * rows.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template <typename SparsityPatternType>
    void
    reinit(const SparsityPatternType &sparsity_pattern);

    /**
     * This function reinitializes the Trilinos sparse matrix from a (possibly
     * distributed) Trilinos sparsity pattern.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     *
     * If you want to write to the matrix from several threads and use MPI,
     * you need to use this reinit method with a sparsity pattern that has
     * been created with explicitly stating writeable rows. In all other
     * cases, you cannot mix MPI with multithreaded writing into the matrix.
     */
    void
    reinit(const SparsityPattern &sparsity_pattern);

    /**
     * This function copies the layout of @p sparse_matrix to the calling
     * matrix. The values are not copied, but you can use copy_from() for
     * this.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    void
    reinit(const SparseMatrix &sparse_matrix);

    /**
     * This function initializes the Trilinos matrix using the deal.II sparse
     * matrix and the entries stored therein. It uses a threshold to copy only
     * elements with modulus larger than the threshold (so zeros in the
     * deal.II matrix can be filtered away).
     *
     * The optional parameter <tt>copy_values</tt> decides whether only the
     * sparsity structure of the input matrix should be used or the matrix
     * entries should be copied, too.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a deadlock.
     *
     * @note If a different sparsity pattern is given in the last argument
     * (i.e., one that differs from the one used in the sparse matrix given in
     * the first argument), then the resulting Trilinos matrix will have the
     * sparsity pattern so given. This of course also means that all entries
     * in the given matrix that are not part of this separate sparsity pattern
     * will in fact be dropped.
     */
    template <typename number>
    void
    reinit(const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
           const double                          drop_tolerance    = 1e-13,
           const bool                            copy_values       = true,
           const ::dealii::SparsityPattern *     use_this_sparsity = nullptr);

    /**
     * This reinit function takes as input a Trilinos Epetra_CrsMatrix and
     * copies its sparsity pattern. If so requested, even the content (values)
     * will be copied.
     */
    void
    reinit(const Epetra_CrsMatrix &input_matrix, const bool copy_values = true);
    //@}

    /**
     * @name Constructors and initialization using an IndexSet description
     */
    //@{
    /**
     * Constructor using an IndexSet and an MPI communicator to describe the
     * %parallel partitioning. The parameter @p n_max_entries_per_row sets the
     * number of nonzero entries in each row that will be allocated. Note that
     * this number does not need to be exact, and it is even allowed that the
     * actual matrix structure has more nonzero entries than specified in the
     * constructor. However it is still advantageous to provide good estimates
     * here since this will considerably increase the performance of the
     * matrix setup. However, there is no effect in the performance of matrix-
     * vector products, since Trilinos reorganizes the matrix memory prior to
     * use (in the compress() step).
     */
    SparseMatrix(const IndexSet &   parallel_partitioning,
                 const MPI_Comm &   communicator          = MPI_COMM_WORLD,
                 const unsigned int n_max_entries_per_row = 0);

    /**
     * Same as before, but now set the number of nonzeros in each matrix row
     * separately. Since we know the number of elements in the matrix exactly
     * in this case, we can already allocate the right amount of memory, which
     * makes the creation process including the insertion of nonzero elements
     * by the respective SparseMatrix::reinit call considerably faster.
     */
    SparseMatrix(const IndexSet &                 parallel_partitioning,
                 const MPI_Comm &                 communicator,
                 const std::vector<unsigned int> &n_entries_per_row);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different IndexSet partitions for row and columns. This interface is
     * meant to be used for generating rectangular matrices, where the first
     * index set describes the %parallel partitioning of the degrees of
     * freedom associated with the matrix rows and the second one the
     * partitioning of the matrix columns. The second index set specifies the
     * partitioning of the vectors this matrix is to be multiplied with, not
     * the distribution of the elements that actually appear in the matrix.
     *
     * The parameter @p n_max_entries_per_row defines how much memory will be
     * allocated for each row. This number does not need to be accurate, as
     * the structure is reorganized in the compress() call.
     */
    SparseMatrix(const IndexSet &row_parallel_partitioning,
                 const IndexSet &col_parallel_partitioning,
                 const MPI_Comm &communicator          = MPI_COMM_WORLD,
                 const size_type n_max_entries_per_row = 0);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different Epetra maps for rows and columns. This interface is meant to
     * be used for generating rectangular matrices, where one map specifies
     * the %parallel distribution of degrees of freedom associated with matrix
     * rows and the second one specifies the %parallel distribution the dofs
     * associated with columns in the matrix. The second map also provides
     * information for the internal arrangement in matrix vector products
     * (i.e., the distribution of vector this matrix is to be multiplied
     * with), but is not used for the distribution of the columns &ndash;
     * rather, all column elements of a row are stored on the same processor
     * in any case. The vector <tt>n_entries_per_row</tt> specifies the number
     * of entries in each row of the newly generated matrix.
     */
    SparseMatrix(const IndexSet &                 row_parallel_partitioning,
                 const IndexSet &                 col_parallel_partitioning,
                 const MPI_Comm &                 communicator,
                 const std::vector<unsigned int> &n_entries_per_row);

    /**
     * This function is initializes the Trilinos Epetra matrix according to
     * the specified sparsity_pattern, and also reassigns the matrix rows to
     * different processes according to a user-supplied index set and
     * %parallel communicator. In programs following the style of the tutorial
     * programs, this function (and the respective call for a rectangular
     * matrix) are the natural way to initialize the matrix size, its
     * distribution among the MPI processes (if run in %parallel) as well as
     * the location of non-zero elements. Trilinos stores the sparsity pattern
     * internally, so it won't be needed any more after this call, in contrast
     * to the deal.II own object. The optional argument @p exchange_data can
     * be used for reinitialization with a sparsity pattern that is not fully
     * constructed. This feature is only implemented for input sparsity
     * patterns of type DynamicSparsityPattern. If the flag is not set, each
     * processor just sets the elements in the sparsity pattern that belong to
     * its rows.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template <typename SparsityPatternType>
    void
    reinit(const IndexSet &           parallel_partitioning,
           const SparsityPatternType &sparsity_pattern,
           const MPI_Comm &           communicator  = MPI_COMM_WORLD,
           const bool                 exchange_data = false);

    /**
     * This function is similar to the other initialization function above,
     * but now also reassigns the matrix rows and columns according to two
     * user-supplied index sets.  To be used for rectangular matrices. The
     * optional argument @p exchange_data can be used for reinitialization
     * with a sparsity pattern that is not fully constructed. This feature is
     * only implemented for input sparsity patterns of type
     * DynamicSparsityPattern.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template <typename SparsityPatternType>
    typename std::enable_if<
      !std::is_same<SparsityPatternType,
                    dealii::SparseMatrix<double>>::value>::type
    reinit(const IndexSet &           row_parallel_partitioning,
           const IndexSet &           col_parallel_partitioning,
           const SparsityPatternType &sparsity_pattern,
           const MPI_Comm &           communicator  = MPI_COMM_WORLD,
           const bool                 exchange_data = false);

    /**
     * This function initializes the Trilinos matrix using the deal.II sparse
     * matrix and the entries stored therein. It uses a threshold to copy only
     * elements with modulus larger than the threshold (so zeros in the
     * deal.II matrix can be filtered away). In contrast to the other reinit
     * function with deal.II sparse matrix argument, this function takes a
     * %parallel partitioning specified by the user instead of internally
     * generating it.
     *
     * The optional parameter <tt>copy_values</tt> decides whether only the
     * sparsity structure of the input matrix should be used or the matrix
     * entries should be copied, too.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template <typename number>
    void
    reinit(const IndexSet &                      parallel_partitioning,
           const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
           const MPI_Comm &                      communicator = MPI_COMM_WORLD,
           const double                          drop_tolerance    = 1e-13,
           const bool                            copy_values       = true,
           const ::dealii::SparsityPattern *     use_this_sparsity = nullptr);

    /**
     * This function is similar to the other initialization function with
     * deal.II sparse matrix input above, but now takes index sets for both
     * the rows and the columns of the matrix. Chosen for rectangular
     * matrices.
     *
     * The optional parameter <tt>copy_values</tt> decides whether only the
     * sparsity structure of the input matrix should be used or the matrix
     * entries should be copied, too.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template <typename number>
    void
    reinit(const IndexSet &                      row_parallel_partitioning,
           const IndexSet &                      col_parallel_partitioning,
           const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
           const MPI_Comm &                      communicator = MPI_COMM_WORLD,
           const double                          drop_tolerance    = 1e-13,
           const bool                            copy_values       = true,
           const ::dealii::SparsityPattern *     use_this_sparsity = nullptr);
    //@}
    /**
     * @name Information on the matrix
     */
    //@{

    /**
     * Return the number of rows in this matrix.
     */
    size_type
    m() const;

    /**
     * Return the number of columns in this matrix.
     */
    size_type
    n() const;

    /**
     * Return the local dimension of the matrix, i.e. the number of rows
     * stored on the present MPI process. For sequential matrices, this number
     * is the same as m(), but for %parallel matrices it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range().
     */
    unsigned int
    local_size() const;

    /**
     * Return a pair of indices indicating which rows of this matrix are
     * stored locally. The first number is the index of the first row stored,
     * the second the index of the one past the last one that is stored
     * locally. If this is a sequential matrix, then the result will be the
     * pair (0,m()), otherwise it will be a pair (i,i+n), where
     * <tt>n=local_size()</tt>.
     */
    std::pair<size_type, size_type>
    local_range() const;

    /**
     * Return whether @p index is in the local range or not, see also
     * local_range().
     */
    bool
    in_local_range(const size_type index) const;

    /**
     * Return the total number of nonzero elements of this matrix (summed
     * over all MPI processes).
     */
    size_type
    n_nonzero_elements() const;

    /**
     * Number of entries in a specific row.
     */
    unsigned int
    row_length(const size_type row) const;

    /**
     * Return the state of the matrix, i.e., whether compress() needs to be
     * called after an operation requiring data exchange. A call to compress()
     * is also needed when the method set() has been called (even when working
     * in serial).
     */
    bool
    is_compressed() const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object. Note that only the memory reserved on the current processor is
     * returned in case this is called in an MPI-based program.
     */
    size_type
    memory_consumption() const;

    /**
     * Return the MPI communicator object in use with this matrix.
     */
    MPI_Comm
    get_mpi_communicator() const;

    //@}
    /**
     * @name Modifying entries
     */
    //@{

    /**
     * This operator assigns a scalar to a matrix. Since this does usually not
     * make much sense (should we set all matrix entries to this value?  Only
     * the nonzero entries of the sparsity pattern?), this operation is only
     * allowed if the actual value to be assigned is zero. This operator only
     * exists to allow for the obvious notation <tt>matrix=0</tt>, which sets
     * all elements of the matrix to zero, but keeps the sparsity pattern
     * previously used.
     */
    SparseMatrix &
    operator=(const double d);

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    void
    clear();

    /**
     * This command does two things:
     * <ul>
     * <li> If the matrix was initialized without a sparsity pattern, elements
     * have been added manually using the set() command. When this process is
     * completed, a call to compress() reorganizes the internal data
     * structures (sparsity pattern) so that a fast access to data is possible
     * in matrix-vector products.
     * <li> If the matrix structure has already been fixed (either by
     * initialization with a sparsity pattern or by calling compress() during
     * the setup phase), this command does the %parallel exchange of data.
     * This is necessary when we perform assembly on more than one (MPI)
     * process, because then some non-local row data will accumulate on nodes
     * that belong to the current's processor element, but are actually held
     * by another. This command is usually called after all elements have been
     * traversed.
     * </ul>
     *
     * In both cases, this function compresses the data structures and allows
     * the resulting matrix to be used in all other operations like matrix-
     * vector products. This is a collective operation, i.e., it needs to be
     * run on all processors when used in %parallel.
     *
     * See
     * @ref GlossCompress "Compressing distributed objects"
     * for more information.
     */
    void
    compress(::dealii::VectorOperation::values operation);

    /**
     * Set the element (<i>i,j</i>) to @p value.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. When compress() is called for the first time (or in case the
     * matrix is initialized from a sparsity pattern), no new elements can be
     * added and an insertion of elements at positions which have not been
     * initialized will throw an exception.
     *
     * For the case that the matrix is constructed without a sparsity pattern
     * and new matrix entries are added on demand, please note the following
     * behavior imposed by the underlying Epetra_FECrsMatrix data structure:
     * If the same matrix entry is inserted more than once, the matrix entries
     * will be added upon calling compress() (since Epetra does not track
     * values to the same entry before the final compress() is called), even
     * if VectorOperation::insert is specified as argument to compress(). In
     * the case you cannot make sure that matrix entries are only set once,
     * initialize the matrix with a sparsity pattern to fix the matrix
     * structure before inserting elements.
     */
    void
    set(const size_type i, const size_type j, const TrilinosScalar value);

    /**
     * Set all elements given in a FullMatrix<double> into the sparse matrix
     * locations given by <tt>indices</tt>. In other words, this function
     * writes the elements in <tt>full_matrix</tt> into the calling matrix,
     * using the local-to-global indexing specified by <tt>indices</tt> for
     * both the rows and the columns of the matrix. This function assumes a
     * quadratic sparse matrix and a quadratic full_matrix, the usual
     * situation in FE calculations.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. After compress() has been called for the first time or the
     * matrix has been initialized from a sparsity pattern, extending the
     * sparsity pattern is no longer possible and an insertion of elements at
     * positions which have not been initialized will throw an exception.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     *
     * For the case that the matrix is constructed without a sparsity pattern
     * and new matrix entries are added on demand, please note the following
     * behavior imposed by the underlying Epetra_FECrsMatrix data structure:
     * If the same matrix entry is inserted more than once, the matrix entries
     * will be added upon calling compress() (since Epetra does not track
     * values to the same entry before the final compress() is called), even
     * if VectorOperation::insert is specified as argument to compress(). In
     * the case you cannot make sure that matrix entries are only set once,
     * initialize the matrix with a sparsity pattern to fix the matrix
     * structure before inserting elements.
     */
    void
    set(const std::vector<size_type> &    indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = false);

    /**
     * Same function as before, but now including the possibility to use
     * rectangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     */
    void
    set(const std::vector<size_type> &    row_indices,
        const std::vector<size_type> &    col_indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = false);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by <tt>col_indices</tt> to the respective value.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. After compress() has been called for the first time or the
     * matrix has been initialized from a sparsity pattern, extending the
     * sparsity pattern is no longer possible and an insertion of elements at
     * positions which have not been initialized will throw an exception.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     *
     * For the case that the matrix is constructed without a sparsity pattern
     * and new matrix entries are added on demand, please note the following
     * behavior imposed by the underlying Epetra_FECrsMatrix data structure:
     * If the same matrix entry is inserted more than once, the matrix entries
     * will be added upon calling compress() (since Epetra does not track
     * values to the same entry before the final compress() is called), even
     * if VectorOperation::insert is specified as argument to compress(). In
     * the case you cannot make sure that matrix entries are only set once,
     * initialize the matrix with a sparsity pattern to fix the matrix
     * structure before inserting elements.
     */
    void
    set(const size_type                    row,
        const std::vector<size_type> &     col_indices,
        const std::vector<TrilinosScalar> &values,
        const bool                         elide_zero_values = false);

    /**
     * Set several elements to values given by <tt>values</tt> in a given row
     * in columns given by col_indices into the sparse matrix.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. After compress() has been called for the first time or the
     * matrix has been initialized from a sparsity pattern, extending the
     * sparsity pattern is no longer possible and an insertion of elements at
     * positions which have not been initialized will throw an exception.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     *
     * For the case that the matrix is constructed without a sparsity pattern
     * and new matrix entries are added on demand, please note the following
     * behavior imposed by the underlying Epetra_FECrsMatrix data structure:
     * If the same matrix entry is inserted more than once, the matrix entries
     * will be added upon calling compress() (since Epetra does not track
     * values to the same entry before the final compress() is called), even
     * if VectorOperation::insert is specified as argument to compress(). In
     * the case you cannot make sure that matrix entries are only set once,
     * initialize the matrix with a sparsity pattern to fix the matrix
     * structure before inserting elements.
     */
    template <typename Number>
    void
    set(const size_type  row,
        const size_type  n_cols,
        const size_type *col_indices,
        const Number *   values,
        const bool       elide_zero_values = false);

    /**
     * Add @p value to the element (<i>i,j</i>).
     *
     * Just as the respective call in deal.II SparseMatrix<Number> class (but
     * in contrast to the situation for PETSc based matrices), this function
     * throws an exception if an entry does not exist in the sparsity pattern.
     * Moreover, if <tt>value</tt> is not a finite number an exception is
     * thrown.
     */
    void
    add(const size_type i, const size_type j, const TrilinosScalar value);

    /**
     * Add all elements given in a FullMatrix<double> into sparse matrix
     * locations given by <tt>indices</tt>. In other words, this function adds
     * the elements in <tt>full_matrix</tt> to the respective entries in
     * calling matrix, using the local-to-global indexing specified by
     * <tt>indices</tt> for both the rows and the columns of the matrix. This
     * function assumes a quadratic sparse matrix and a quadratic full_matrix,
     * the usual situation in FE calculations.
     *
     * Just as the respective call in deal.II SparseMatrix<Number> class (but
     * in contrast to the situation for PETSc based matrices), this function
     * throws an exception if an entry does not exist in the sparsity pattern.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     */
    void
    add(const std::vector<size_type> &    indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = true);

    /**
     * Same function as before, but now including the possibility to use
     * rectangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     */
    void
    add(const std::vector<size_type> &    row_indices,
        const std::vector<size_type> &    col_indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = true);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by <tt>col_indices</tt> to the respective value.
     *
     * Just as the respective call in deal.II SparseMatrix<Number> class (but
     * in contrast to the situation for PETSc based matrices), this function
     * throws an exception if an entry does not exist in the sparsity pattern.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     */
    void
    add(const size_type                    row,
        const std::vector<size_type> &     col_indices,
        const std::vector<TrilinosScalar> &values,
        const bool                         elide_zero_values = true);

    /**
     * Add an array of values given by <tt>values</tt> in the given global
     * matrix row at columns specified by col_indices in the sparse matrix.
     *
     * Just as the respective call in deal.II SparseMatrix<Number> class (but
     * in contrast to the situation for PETSc based matrices), this function
     * throws an exception if an entry does not exist in the sparsity pattern.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     */
    void
    add(const size_type       row,
        const size_type       n_cols,
        const size_type *     col_indices,
        const TrilinosScalar *values,
        const bool            elide_zero_values      = true,
        const bool            col_indices_are_sorted = false);

    /**
     * Multiply the entire matrix by a fixed factor.
     */
    SparseMatrix &
    operator*=(const TrilinosScalar factor);

    /**
     * Divide the entire matrix by a fixed factor.
     */
    SparseMatrix &
    operator/=(const TrilinosScalar factor);

    /**
     * Copy the given (Trilinos) matrix (sparsity pattern and entries).
     */
    void
    copy_from(const SparseMatrix &source);

    /**
     * Add <tt>matrix</tt> scaled by <tt>factor</tt> to this matrix, i.e. the
     * matrix <tt>factor*matrix</tt> is added to <tt>this</tt>. If the
     * sparsity pattern of the calling matrix does not contain all the
     * elements in the sparsity pattern of the input matrix, this function
     * will throw an exception.
     */
    void
    add(const TrilinosScalar factor, const SparseMatrix &matrix);

    /**
     * Remove all elements from this <tt>row</tt> by setting them to zero. The
     * function does not modify the number of allocated nonzero entries, it
     * only sets the entries to zero.
     *
     * This operation is used in eliminating constraints (e.g. due to hanging
     * nodes) and makes sure that we can write this modification to the matrix
     * without having to read entries (such as the locations of non-zero
     * elements) from it &mdash; without this operation, removing constraints
     * on %parallel matrices is a rather complicated procedure.
     *
     * The second parameter can be used to set the diagonal entry of this row
     * to a value different from zero. The default is to set it to zero.
     *
     * @note If the matrix is stored in parallel across multiple processors
     * using MPI, this function only touches rows that are locally stored and
     * simply ignores all other row indices. Further, in the context of
     * parallel computations, you will get into trouble if you clear a row
     * while other processors still have pending writes or additions into the
     * same row. In other words, if another processor still wants to add
     * something to an element of a row and you call this function to zero out
     * the row, then the next time you call compress() may add the remote
     * value to the zero you just created. Consequently, you will want to call
     * compress() after you made the last modifications to a matrix and before
     * starting to clear rows.
     */
    void
    clear_row(const size_type row, const TrilinosScalar new_diag_value = 0);

    /**
     * Same as clear_row(), except that it works on a number of rows at once.
     *
     * The second parameter can be used to set the diagonal entries of all
     * cleared rows to something different from zero. Note that all of these
     * diagonal entries get the same value -- if you want different values for
     * the diagonal entries, you have to set them by hand.
     *
     * @note If the matrix is stored in parallel across multiple processors
     * using MPI, this function only touches rows that are locally stored and
     * simply ignores all other row indices. Further, in the context of
     * parallel computations, you will get into trouble if you clear a row
     * while other processors still have pending writes or additions into the
     * same row. In other words, if another processor still wants to add
     * something to an element of a row and you call this function to zero out
     * the row, then the next time you call compress() may add the remote
     * value to the zero you just created. Consequently, you will want to call
     * compress() after you made the last modifications to a matrix and before
     * starting to clear rows.
     */
    void
    clear_rows(const std::vector<size_type> &rows,
               const TrilinosScalar          new_diag_value = 0);

    /**
     * Sets an internal flag so that all operations performed by the matrix,
     * i.e., multiplications, are done in transposed order. However, this does
     * not reshape the matrix to transposed form directly, so care should be
     * taken when using this flag.
     *
     * @note Calling this function any even number of times in succession will
     * return the object to its original state.
     */
    void
    transpose();

    //@}
    /**
     * @name Entry Access
     */
    //@{

    /**
     * Return the value of the entry (<i>i,j</i>).  This may be an expensive
     * operation and you should always take care where to call this function.
     * As in the deal.II sparse matrix class, we throw an exception if the
     * respective entry doesn't exist in the sparsity pattern of this class,
     * which is requested from Trilinos. Moreover, an exception will be thrown
     * when the requested element is not saved on the calling process.
     */
    TrilinosScalar
    operator()(const size_type i, const size_type j) const;

    /**
     * Return the value of the matrix entry (<i>i,j</i>). If this entry does
     * not exist in the sparsity pattern, then zero is returned. While this
     * may be convenient in some cases, note that it is simple to write
     * algorithms that are slow compared to an optimal solution, since the
     * sparsity of the matrix is not used.  On the other hand, if you want to
     * be sure the entry exists, you should use operator() instead.
     *
     * The lack of error checking in this function can also yield surprising
     * results if you have a parallel matrix. In that case, just because you
     * get a zero result from this function does not mean that either the
     * entry does not exist in the sparsity pattern or that it does but has a
     * value of zero. Rather, it could also be that it simply isn't stored on
     * the current processor; in that case, it may be stored on a different
     * processor, and possibly so with a nonzero value.
     */
    TrilinosScalar
    el(const size_type i, const size_type j) const;

    /**
     * Return the main diagonal element in the <i>i</i>th row. This function
     * throws an error if the matrix is not quadratic and it also throws an
     * error if <i>(i,i)</i> is not element of the local matrix.  See also the
     * comment in trilinos_sparse_matrix.cc.
     */
    TrilinosScalar
    diag_element(const size_type i) const;

    //@}
    /**
     * @name Multiplications
     */
    //@{

    /**
     * Matrix-vector multiplication: let <i>dst = M*src</i> with <i>M</i>
     * being this matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several types of vector objects,
     * namely @p VectorType can be
     * <ul>
     * <li> TrilinosWrappers::MPI::Vector,
     * <li> LinearAlgebra::EpetraWrappers::Vector,
     * <li> LinearAlgebra::TpetraWrappers::Vector,
     * <li> Vector<double>,
     * <li> LinearAlgebra::distributed::Vector<double>.
     * </ul>
     *
     * When using vectors of type TrilinosWrappers::MPI::Vector, the vector
     * @p dst has to be initialized with the same IndexSet that was used for
     * the row indices of the matrix and the vector @p src has to be
     * initialized with the same IndexSet that was used for the column indices
     * of the matrix.
     *
     * This function will be called when the underlying number type for the
     * matrix object and the one for the vector object are the same.
     * Despite looking complicated, the return type is just `void`.
     *
     * In case of a serial vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, an exception will be thrown.
     */
    template <typename VectorType>
    typename std::enable_if<std::is_same<typename VectorType::value_type,
                                         TrilinosScalar>::value>::type
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Same as the function above for the case that the underlying number type
     * for the matrix object and the one for the vector object do not coincide.
     * This case is not implemented. Calling it will result in a runtime error.
     * Despite looking complicated, the return type is just `void`.
     */
    template <typename VectorType>
    typename std::enable_if<!std::is_same<typename VectorType::value_type,
                                          TrilinosScalar>::value>::type
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Matrix-vector multiplication: let <i>dst = M<sup>T</sup>*src</i> with
     * <i>M</i> being this matrix. This function does the same as vmult() but
     * takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several types of vector objects,
     * see the discussion about @p VectorType in vmult().
     *
     * This function will be called when the underlying number type for the
     * matrix object and the one for the vector object are the same.
     * Despite looking complicated, the return type is just `void`.
     */
    template <typename VectorType>
    typename std::enable_if<std::is_same<typename VectorType::value_type,
                                         TrilinosScalar>::value>::type
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * Same as the function above for the case that the underlying number type
     * for the matrix object and the one for the vector object do not coincide.
     * This case is not implemented. Calling it will result in a runtime error.
     * Despite looking complicated, the return type is just `void`.
     */
    template <typename VectorType>
    typename std::enable_if<!std::is_same<typename VectorType::value_type,
                                          TrilinosScalar>::value>::type
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * Adding matrix-vector multiplication. Add <i>M*src</i> on <i>dst</i>
     * with <i>M</i> being this matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several types of vector objects,
     * see the discussion about @p VectorType in vmult().
     */
    template <typename VectorType>
    void
    vmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * Adding matrix-vector multiplication. Add <i>M<sup>T</sup>*src</i> to
     * <i>dst</i> with <i>M</i> being this matrix. This function does the same
     * as vmult_add() but takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several types of vector objects,
     * see the discussion about @p VectorType in vmult().
     */
    template <typename VectorType>
    void
    Tvmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * Return the square of the norm of the vector $v$ with respect to the
     * norm induced by this matrix, i.e., $\left(v,Mv\right)$. This is useful,
     * e.g. in the finite element context, where the $L_2$ norm of a function
     * equals the matrix norm with respect to the mass matrix of the vector
     * representing the nodal values of the finite element function.
     *
     * Obviously, the matrix needs to be quadratic for this operation.
     *
     * The implementation of this function is not as efficient as the one in
     * the @p SparseMatrix class used in deal.II (i.e. the original one, not
     * the Trilinos wrapper class) since Trilinos doesn't support this
     * operation and needs a temporary vector.
     *
     * The vector has to be initialized with the same IndexSet the matrix
     * was initialized with.
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, an exception will be thrown.
     */
    TrilinosScalar
    matrix_norm_square(const MPI::Vector &v) const;

    /**
     * Compute the matrix scalar product $\left(u,Mv\right)$.
     *
     * The implementation of this function is not as efficient as the one in
     * the @p SparseMatrix class used in deal.II (i.e. the original one, not
     * the Trilinos wrapper class) since Trilinos doesn't support this
     * operation and needs a temporary vector.
     *
     * The vector @p u has to be initialized with the same IndexSet that
     * was used for the row indices of the matrix and the vector @p v has
     * to be initialized with the same IndexSet that was used for the
     * column indices of the matrix.
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, an exception will be thrown.
     *
     * This function is only implemented for square matrices.
     */
    TrilinosScalar
    matrix_scalar_product(const MPI::Vector &u, const MPI::Vector &v) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * Source <i>x</i> and destination <i>dst</i> must not be the same vector.
     *
     * The vectors @p dst and @p b have to be initialized with the same
     * IndexSet that was used for the row indices of the matrix and the vector
     * @p x has to be initialized with the same IndexSet that was used for the
     * column indices of the matrix.
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, an exception will be thrown.
     */
    TrilinosScalar
    residual(MPI::Vector &      dst,
             const MPI::Vector &x,
             const MPI::Vector &b) const;

    /**
     * Perform the matrix-matrix multiplication <tt>C = A * B</tt>, or, if an
     * optional vector argument is given, <tt>C = A * diag(V) * B</tt>, where
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
     */
    void
    mmult(SparseMatrix &      C,
          const SparseMatrix &B,
          const MPI::Vector & V = MPI::Vector()) const;


    /**
     * Perform the matrix-matrix multiplication with the transpose of
     * <tt>this</tt>, i.e., <tt>C = A<sup>T</sup> * B</tt>, or, if an optional
     * vector argument is given, <tt>C = A<sup>T</sup> * diag(V) * B</tt>,
     * where <tt>diag(V)</tt> defines a diagonal matrix with the vector
     * entries.
     *
     * This function assumes that the calling matrix <tt>A</tt> and <tt>B</tt>
     * have compatible sizes. The size of <tt>C</tt> will be set within this
     * function.
     *
     * The content as well as the sparsity pattern of the matrix C will be
     * changed by this function, so make sure that the sparsity pattern is not
     * used somewhere else in your program. This is an expensive operation, so
     * think twice before you use this function.
     */
    void
    Tmmult(SparseMatrix &      C,
           const SparseMatrix &B,
           const MPI::Vector & V = MPI::Vector()) const;

    //@}
    /**
     * @name Matrix norms
     */
    //@{

    /**
     * Return the <i>l</i><sub>1</sub>-norm of the matrix, that is $|M|_1=
     * \max_{\mathrm{all\ columns\ } j} \sum_{\mathrm{all\ rows\ } i}
     * |M_{ij}|$, (max. sum of columns).  This is the natural matrix norm that
     * is compatible to the l1-norm for vectors, i.e.  $|Mv|_1 \leq |M|_1
     * |v|_1$.  (cf. Haemmerlin-Hoffmann: Numerische Mathematik)
     */
    TrilinosScalar
    l1_norm() const;

    /**
     * Return the linfty-norm of the matrix, that is
     * $|M|_\infty=\max_{\mathrm{all\ rows\ } i}\sum_{\mathrm{all\ columns\ }
     * j} |M_{ij}|$, (max. sum of rows).  This is the natural matrix norm that
     * is compatible to the linfty-norm of vectors, i.e.  $|Mv|_\infty \leq
     * |M|_\infty |v|_\infty$.  (cf. Haemmerlin-Hoffmann: Numerische
     * Mathematik)
     */
    TrilinosScalar
    linfty_norm() const;

    /**
     * Return the frobenius norm of the matrix, i.e. the square root of the
     * sum of squares of all entries in the matrix.
     */
    TrilinosScalar
    frobenius_norm() const;

    //@}
    /**
     * @name Access to underlying Trilinos data
     */
    //@{

    /**
     * Return a const reference to the underlying Trilinos Epetra_CrsMatrix
     * data.
     */
    const Epetra_CrsMatrix &
    trilinos_matrix() const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_CrsGraph
     * data that stores the sparsity pattern of the matrix.
     */
    const Epetra_CrsGraph &
    trilinos_sparsity_pattern() const;

    //@}

    /**
     * @name Partitioners
     */
    //@{

    /**
     * Return the partitioning of the domain space of this matrix, i.e., the
     * partitioning of the vectors this matrix has to be multiplied with.
     */
    IndexSet
    locally_owned_domain_indices() const;

    /**
     * Return the partitioning of the range space of this matrix, i.e., the
     * partitioning of the vectors that are result from matrix-vector
     * products.
     */
    IndexSet
    locally_owned_range_indices() const;

    //@}

    /**
     * @name Iterators
     */
    //@{

    /**
     * Return an iterator pointing to the first element of the matrix.
     *
     * The elements accessed by iterators within each row are ordered in the
     * way in which Trilinos stores them, though the implementation guarantees
     * that all elements of one row are accessed before the elements of the
     * next row. If your algorithm relies on visiting elements within one row,
     * you will need to consult with the Trilinos documentation on the order
     * in which it stores data. It is, however, generally not a good and long-
     * term stable idea to rely on the order in which receive elements if you
     * iterate over them.
     *
     * When you iterate over the elements of a parallel matrix, you will only
     * be able to access the locally owned rows. (You can access the other
     * rows as well, but they will look empty.) In that case, you probably
     * want to call the begin() function that takes the row as an argument to
     * limit the range of elements to loop over.
     */
    const_iterator
    begin() const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator
    begin();

    /**
     * Return an iterator pointing the element past the last one of this
     * matrix.
     */
    const_iterator
    end() const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator
    end();

    /**
     * Return an iterator pointing to the first element of row @p r.
     *
     * Note that if the given row is empty, i.e. does not contain any nonzero
     * entries, then the iterator returned by this function equals
     * <tt>end(r)</tt>. The returned iterator may not be dereferenceable in
     * that case if neither row @p r nor any of the following rows contain any
     * nonzero entries.
     *
     * The elements accessed by iterators within each row are ordered in the
     * way in which Trilinos stores them, though the implementation guarantees
     * that all elements of one row are accessed before the elements of the
     * next row. If your algorithm relies on visiting elements within one row,
     * you will need to consult with the Trilinos documentation on the order
     * in which it stores data. It is, however, generally not a good and long-
     * term stable idea to rely on the order in which receive elements if you
     * iterate over them.
     *
     * @note When you access the elements of a parallel matrix, you can only
     * access the elements of rows that are actually stored locally. (You can
     * access the other rows as well, but they will look empty.) Even then, if
     * another processor has since written into, or added to, an element of
     * the matrix that is stored on the current processor, then you will still
     * see the old value of this entry unless you have called compress()
     * between modifying the matrix element on the remote processor and
     * accessing it on the current processor. See the documentation of the
     * compress() function for more information.
     */
    const_iterator
    begin(const size_type r) const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator
    begin(const size_type r);

    /**
     * Return an iterator pointing the element past the last one of row @p r ,
     * or past the end of the entire sparsity pattern if none of the rows
     * after @p r contain any entries at all.
     *
     * Note that the end iterator is not necessarily dereferenceable. This is
     * in particular the case if it is the end iterator for the last row of a
     * matrix.
     */
    const_iterator
    end(const size_type r) const;

    /**
     * Like the function above, but for non-const matrices.
     */
    iterator
    end(const size_type r);

    //@}
    /**
     * @name Input/Output
     */
    //@{

    /**
     * Abstract Trilinos object that helps view in ASCII other Trilinos
     * objects. Currently this function is not implemented.  TODO: Not
     * implemented.
     */
    void
    write_ascii();

    /**
     * Print the matrix to the given stream, using the format <tt>(line,col)
     * value</tt>, i.e. one nonzero entry of the matrix per line. The optional
     * flag outputs the sparsity pattern in Trilinos style, where the data is
     * sorted according to the processor number when printed to the stream, as
     * well as a summary of the matrix like the global size.
     */
    void
    print(std::ostream &out,
          const bool    write_extended_trilinos_info = false) const;

    //@}
    /**
     * @addtogroup Exceptions
     */
    //@{
    /**
     * Exception
     */
    DeclException1(ExcTrilinosError,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while calling a Trilinos function");

    /**
     * Exception
     */
    DeclException2(ExcInvalidIndex,
                   size_type,
                   size_type,
                   << "The entry with index <" << arg1 << ',' << arg2
                   << "> does not exist.");

    /**
     * Exception
     */
    DeclExceptionMsg(ExcSourceEqualsDestination,
                     "You are attempting an operation on two matrices that "
                     "are the same object, but the operation requires that the "
                     "two objects are in fact different.");

    /**
     * Exception
     */
    DeclException0(ExcMatrixNotCompressed);

    /**
     * Exception
     */
    DeclException4(ExcAccessToNonLocalElement,
                   size_type,
                   size_type,
                   size_type,
                   size_type,
                   << "You tried to access element (" << arg1 << "/" << arg2
                   << ")"
                   << " of a distributed matrix, but only rows in range ["
                   << arg3 << "," << arg4
                   << "] are stored locally and can be accessed.");

    /**
     * Exception
     */
    DeclException2(ExcAccessToNonPresentElement,
                   size_type,
                   size_type,
                   << "You tried to access element (" << arg1 << "/" << arg2
                   << ")"
                   << " of a sparse matrix, but it appears to not"
                   << " exist in the Trilinos sparsity pattern.");
    //@}



  protected:
    /**
     * For some matrix storage formats, in particular for the PETSc
     * distributed blockmatrices, set and add operations on individual
     * elements can not be freely mixed. Rather, one has to synchronize
     * operations when one wants to switch from setting elements to adding to
     * elements.  BlockMatrixBase automatically synchronizes the access by
     * calling this helper function for each block.  This function ensures
     * that the matrix is in a state that allows adding elements; if it
     * previously already was in this state, the function does nothing.
     */
    void
    prepare_add();

    /**
     * Same as prepare_add() but prepare the matrix for setting elements if
     * the representation of elements in this class requires such an
     * operation.
     */
    void
    prepare_set();



  private:
    /**
     * Pointer to the user-supplied Epetra Trilinos mapping of the matrix
     * columns that assigns parts of the matrix to the individual processes.
     */
    std::unique_ptr<Epetra_Map> column_space_map;

    /**
     * A sparse matrix object in Trilinos to be used for finite element based
     * problems which allows for assembling into non-local elements.  The
     * actual type, a sparse matrix, is set in the constructor.
     */
    std::unique_ptr<Epetra_FECrsMatrix> matrix;

    /**
     * A sparse matrix object in Trilinos to be used for collecting the non-
     * local elements if the matrix was constructed from a Trilinos sparsity
     * pattern with the respective option.
     */
    std::unique_ptr<Epetra_CrsMatrix> nonlocal_matrix;

    /**
     * An export object used to communicate the nonlocal matrix.
     */
    std::unique_ptr<Epetra_Export> nonlocal_matrix_exporter;

    /**
     * Trilinos doesn't allow to mix additions to matrix entries and
     * overwriting them (to make synchronization of %parallel computations
     * simpler). The way we do it is to, for each access operation, store
     * whether it is an insertion or an addition. If the previous one was of
     * different type, then we first have to flush the Trilinos buffers;
     * otherwise, we can simply go on. Luckily, Trilinos has an object for
     * this which does already all the %parallel communications in such a
     * case, so we simply use their model, which stores whether the last
     * operation was an addition or an insertion.
     */
    Epetra_CombineMode last_action;

    /**
     * A boolean variable to hold information on whether the vector is
     * compressed or not.
     */
    bool compressed;

    // To allow calling protected prepare_add() and prepare_set().
    friend class BlockMatrixBase<SparseMatrix>;
  };



  // forwards declarations
  class SolverBase;
  class PreconditionBase;

  namespace internal
  {
    inline void
    check_vector_map_equality(const Epetra_CrsMatrix &  mtrx,
                              const Epetra_MultiVector &src,
                              const Epetra_MultiVector &dst,
                              const bool                transpose)
    {
      if (transpose == false)
        {
          Assert(src.Map().SameAs(mtrx.DomainMap()) == true,
                 ExcMessage(
                   "Column map of matrix does not fit with vector map!"));
          Assert(dst.Map().SameAs(mtrx.RangeMap()) == true,
                 ExcMessage("Row map of matrix does not fit with vector map!"));
        }
      else
        {
          Assert(src.Map().SameAs(mtrx.RangeMap()) == true,
                 ExcMessage(
                   "Column map of matrix does not fit with vector map!"));
          Assert(dst.Map().SameAs(mtrx.DomainMap()) == true,
                 ExcMessage("Row map of matrix does not fit with vector map!"));
        }
      (void)mtrx; // removes -Wunused-variable in optimized mode
      (void)src;
      (void)dst;
    }

    inline void
    check_vector_map_equality(const Epetra_Operator &   op,
                              const Epetra_MultiVector &src,
                              const Epetra_MultiVector &dst,
                              const bool                transpose)
    {
      if (transpose == false)
        {
          Assert(src.Map().SameAs(op.OperatorDomainMap()) == true,
                 ExcMessage(
                   "Column map of operator does not fit with vector map!"));
          Assert(dst.Map().SameAs(op.OperatorRangeMap()) == true,
                 ExcMessage(
                   "Row map of operator does not fit with vector map!"));
        }
      else
        {
          Assert(src.Map().SameAs(op.OperatorRangeMap()) == true,
                 ExcMessage(
                   "Column map of operator does not fit with vector map!"));
          Assert(dst.Map().SameAs(op.OperatorDomainMap()) == true,
                 ExcMessage(
                   "Row map of operator does not fit with vector map!"));
        }
      (void)op; // removes -Wunused-variable in optimized mode
      (void)src;
      (void)dst;
    }

    namespace LinearOperatorImplementation
    {
      /**
       * This is an extension class to LinearOperators for Trilinos sparse
       * matrix and preconditioner types. It provides the interface to
       * performing basic operations (<tt>vmult</tt> and <tt>Tvmult</tt>)  on
       * Trilinos vector types. It fulfills the requirements necessary for
       * wrapping a Trilinos solver, which calls Epetra_Operator functions, as a
       * LinearOperator.
       *
       * @note The TrilinosWrappers::SparseMatrix or
       * TrilinosWrappers::PreconditionBase that this payload wraps is passed by
       * reference to the <tt>vmult</tt> and <tt>Tvmult</tt> functions. This
       * object is not thread-safe when the transpose flag is set on it or the
       * Trilinos object to which it refers. See the docuemtation for the
       * TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload::SetUseTranspose()
       * function for further details.
       *
       *
       * @ingroup TrilinosWrappers
       */
      class TrilinosPayload : public Epetra_Operator
      {
      public:
        /**
         * Definition for the internally supported vector type.
         */
        using VectorType = Epetra_MultiVector;

        /**
         * Definition for the vector type for the domain space of the operator.
         */
        using Range = VectorType;

        /**
         * Definition for the vector type for the range space of the operator.
         */
        using Domain = VectorType;

        /**
         * @name Constructors / destructor
         */
        //@{

        /**
         * Default constructor
         *
         * @note By design, the resulting object is inoperable since there is
         * insufficient information with which to construct the domain and
         * range maps.
         */
        TrilinosPayload();

        /**
         * Constructor for a sparse matrix based on an exemplary matrix
         */
        TrilinosPayload(const TrilinosWrappers::SparseMatrix &matrix_exemplar,
                        const TrilinosWrappers::SparseMatrix &matrix);

        /**
         * Constructor for a preconditioner based on an exemplary matrix
         */
        TrilinosPayload(
          const TrilinosWrappers::SparseMatrix &    matrix_exemplar,
          const TrilinosWrappers::PreconditionBase &preconditioner);

        /**
         * Constructor for a preconditioner based on an exemplary preconditioner
         */
        TrilinosPayload(
          const TrilinosWrappers::PreconditionBase &preconditioner_exemplar,
          const TrilinosWrappers::PreconditionBase &preconditioner);

        /**
         * Default copy constructor
         */
        TrilinosPayload(const TrilinosPayload &payload);

        /**
         * Composite copy constructor
         *
         * This is required for PackagedOperations as it sets up the domain and
         * range maps, and composite <tt>vmult</tt> and <tt>Tvmult</tt>
         * operations based on the combined operation of both operations
         */
        TrilinosPayload(const TrilinosPayload &first_op,
                        const TrilinosPayload &second_op);

        /**
         * Destructor
         */
        virtual ~TrilinosPayload() override = default;

        /**
         * Return a payload configured for identity operations
         */
        TrilinosPayload
        identity_payload() const;

        /**
         * Return a payload configured for null operations
         */
        TrilinosPayload
        null_payload() const;

        /**
         * Return a payload configured for transpose operations
         */
        TrilinosPayload
        transpose_payload() const;

        /**
         * Return a payload configured for inverse operations
         *
         * Invoking this factory function will configure two additional
         * functions, namely <tt>inv_vmult</tt> and <tt>inv_Tvmult</tt>, both of
         * which wrap inverse operations. The <tt>vmult</tt> and <tt>Tvmult</tt>
         * operations retain the standard
         * definitions inherited from @p op.
         *
         * @note This function is enabled only if the solver and preconditioner
         * derive from the respective TrilinosWrappers base classes.
         * The C++ compiler will therefore only consider this function if the
         * following criterion are satisfied:
         * 1. the @p Solver derives from TrilinosWrappers::SolverBase, and
         * 2. the @p Preconditioner derives from TrilinosWrappers::PreconditionBase.
         */
        template <typename Solver, typename Preconditioner>
        typename std::enable_if<
          std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
            std::is_base_of<TrilinosWrappers::PreconditionBase,
                            Preconditioner>::value,
          TrilinosPayload>::type
        inverse_payload(Solver &, const Preconditioner &) const;

        /**
         * Return a payload configured for inverse operations
         *
         * Invoking this factory function will configure two additional
         * functions, namely <tt>inv_vmult</tt> and <tt>inv_Tvmult</tt>, both of
         * which
         * are disabled because the @p Solver or @p Preconditioner are not
         * compatible with Epetra_MultiVector.
         * The <tt>vmult</tt> and <tt>Tvmult</tt> operations retain the standard
         * definitions inherited from @p op.
         *
         * @note The C++ compiler will only consider this function if the
         * following criterion are satisfied:
         * 1. the @p Solver does not derive from TrilinosWrappers::SolverBase, and
         * 2. the @p Preconditioner does not derive from
         * TrilinosWrappers::PreconditionBase.
         */
        template <typename Solver, typename Preconditioner>
        typename std::enable_if<
          !(std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
            std::is_base_of<TrilinosWrappers::PreconditionBase,
                            Preconditioner>::value),
          TrilinosPayload>::type
        inverse_payload(Solver &, const Preconditioner &) const;

        //@}

        /**
         * @name LinearOperator functionality
         */
        //@{

        /**
         * Return an IndexSet that defines the partitioning of the domain space
         * of this matrix, i.e., the partitioning of the vectors this matrix has
         * to be multiplied with / operate on.
         */
        IndexSet
        locally_owned_domain_indices() const;

        /**
         * Return an IndexSet that defines the partitioning of the range space
         * of this matrix, i.e., the partitioning of the vectors that result
         * from matrix-vector products.
         */
        IndexSet
        locally_owned_range_indices() const;

        /**
         * Return the MPI communicator object in use with this Payload.
         */
        MPI_Comm
        get_mpi_communicator() const;

        /**
         * Sets an internal flag so that all operations performed by the matrix,
         * i.e., multiplications, are done in transposed order.
         * @note This does not reshape the matrix to transposed form directly,
         * so care should be taken when using this flag.
         */
        void
        transpose();

        /**
         * The standard matrix-vector operation to be performed by the payload
         * when Apply is called.
         *
         * @note This is not called by a LinearOperator, but rather by Trilinos
         * functions that expect this to mimic the action of the LinearOperator.
         */
        std::function<void(VectorType &, const VectorType &)> vmult;

        /**
         * The standard transpose matrix-vector operation to be performed by
         * the payload when Apply is called.
         *
         * @note This is not called by a LinearOperator, but rather by Trilinos
         * functions that expect this to mimic the action of the LinearOperator.
         */
        std::function<void(VectorType &, const VectorType &)> Tvmult;

        /**
         * The inverse matrix-vector operation to be performed by the payload
         * when ApplyInverse is called.
         *
         * @note This is not called by a LinearOperator, but rather by Trilinos
         * functions that expect this to mimic the action of the
         * InverseOperator.
         */
        std::function<void(VectorType &, const VectorType &)> inv_vmult;

        /**
         * The inverse transpose matrix-vector operation to be performed by
         * the payload when ApplyInverse is called.
         *
         * @note This is not called by a LinearOperator, but rather by Trilinos
         * functions that expect this to mimic the action of the
         * InverseOperator.
         */
        std::function<void(VectorType &, const VectorType &)> inv_Tvmult;

        //@}

        /**
         * @name Core Epetra_Operator functionality
         */
        //@{

        /**
         * Return the status of the transpose flag for this operator
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual bool
        UseTranspose() const override;

        /**
         * Sets an internal flag so that all operations performed by the matrix,
         * i.e., multiplications, are done in transposed order.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         *
         * @note This does not reshape the matrix to transposed form directly,
         * so care should be taken when using this flag. When the flag is set to
         * true (either here or directly on the underlying Trilinos object
         * itself), this object is no longer thread-safe. In essence, it is not
         * possible ensure that the transposed state of the LinearOperator and
         * the underlying Trilinos object remain synchronized throughout all
         * operations that may occur on different threads simultaneously.
         */
        virtual int
        SetUseTranspose(bool UseTranspose) override;

        /**
         * Apply the vmult operation on a vector @p X (of internally defined
         * type VectorType) and store the result in the vector @p Y.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         *
         * @note The intended operation depends on the status of the internal
         * transpose flag. If this flag is set to true, the result will be
         * the equivalent of performing a Tvmult operation.
         */
        virtual int
        Apply(const VectorType &X, VectorType &Y) const override;

        /**
         * Apply the vmult inverse operation on a vector @p X (of internally
         * defined type VectorType) and store the result in the vector @p Y.
         *
         * In practise, this function is only called from a Trilinos solver if
         * the wrapped object is to act as a preconditioner.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         *
         * @note This function will only be operable if the payload has been
         * initialized with an InverseOperator, or is a wrapper to a
         * preconditioner. If not, then using this function will lead to an
         * error being thrown.
         * @note The intended operation depends on the status of the internal
         * transpose flag. If this flag is set to true, the result will be
         * the equivalent of performing a Tvmult operation.
         */
        virtual int
        ApplyInverse(const VectorType &Y, VectorType &X) const override;
        //@}

        /**
         * @name Additional Epetra_Operator functionality
         */
        //@{

        /**
         * Return a label to describe this class.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual const char *
        Label() const override;

        /**
         * Return a reference to the underlying MPI communicator for
         * this object.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual const Epetra_Comm &
        Comm() const override;

        /**
         * Return the partitioning of the domain space of this matrix, i.e., the
         * partitioning of the vectors this matrix has to be multiplied with.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual const Epetra_Map &
        OperatorDomainMap() const override;

        /**
         * Return the partitioning of the range space of this matrix, i.e., the
         * partitioning of the vectors that are result from matrix-vector
         * products.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual const Epetra_Map &
        OperatorRangeMap() const override;
        //@}

      private:
        /**
         * A flag recording whether the operator is to perform standard
         * matrix-vector multiplication, or the transpose operation.
         */
        bool use_transpose;

        /**
         * Internal communication pattern in case the matrix needs to be copied
         * from deal.II format.
         */
#    ifdef DEAL_II_WITH_MPI
        Epetra_MpiComm communicator;
#    else
        Epetra_SerialComm communicator;
#    endif

        /**
         * Epetra_Map that sets the partitioning of the domain space of
         * this operator.
         */
        Epetra_Map domain_map;

        /**
         * Epetra_Map that sets the partitioning of the range space of
         * this operator.
         */
        Epetra_Map range_map;

        /**
         * Return a flag that describes whether this operator can return the
         * computation of the infinity norm. Since in general this is not the
         * case, this always returns a negetive result.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual bool
        HasNormInf() const override;

        /**
         * Return the infinity norm of this operator.
         * Throws an error since, in general, we cannot compute this value.
         *
         * This overloads the same function from the Trilinos class
         * Epetra_Operator.
         */
        virtual double
        NormInf() const override;
      };

      /**
       * Return an operator that returns a payload configured to support the
       * addition of two LinearOperators
       */
      TrilinosPayload
      operator+(const TrilinosPayload &first_op,
                const TrilinosPayload &second_op);

      /**
       * Return an operator that returns a payload configured to support the
       * multiplication of two LinearOperators
       */
      TrilinosPayload operator*(const TrilinosPayload &first_op,
                                const TrilinosPayload &second_op);

    } // namespace LinearOperatorImplementation
  }   /* namespace internal */



  // ----------------------- inline and template functions --------------------

#    ifndef DOXYGEN

  namespace SparseMatrixIterators
  {
    inline AccessorBase::AccessorBase(SparseMatrix *matrix,
                                      size_type     row,
                                      size_type     index)
      : matrix(matrix)
      , a_row(row)
      , a_index(index)
    {
      visit_present_row();
    }


    inline AccessorBase::size_type
    AccessorBase::row() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_row;
    }


    inline AccessorBase::size_type
    AccessorBase::column() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*colnum_cache)[a_index];
    }


    inline AccessorBase::size_type
    AccessorBase::index() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_index;
    }


    inline Accessor<true>::Accessor(MatrixType *    matrix,
                                    const size_type row,
                                    const size_type index)
      : AccessorBase(const_cast<SparseMatrix *>(matrix), row, index)
    {}


    template <bool Other>
    inline Accessor<true>::Accessor(const Accessor<Other> &other)
      : AccessorBase(other)
    {}


    inline TrilinosScalar
    Accessor<true>::value() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*value_cache)[a_index];
    }


    inline Accessor<false>::Reference::Reference(const Accessor<false> &acc)
      : accessor(const_cast<Accessor<false> &>(acc))
    {}


    inline Accessor<false>::Reference::operator TrilinosScalar() const
    {
      return (*accessor.value_cache)[accessor.a_index];
    }

    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] = n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator+=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] += n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator-=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] -= n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator*=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] *= n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator/=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] /= n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline Accessor<false>::Accessor(MatrixType *    matrix,
                                     const size_type row,
                                     const size_type index)
      : AccessorBase(matrix, row, index)
    {}


    inline Accessor<false>::Reference
    Accessor<false>::value() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return {*this};
    }



    template <bool Constness>
    inline Iterator<Constness>::Iterator(MatrixType *    matrix,
                                         const size_type row,
                                         const size_type index)
      : accessor(matrix, row, index)
    {}


    template <bool Constness>
    template <bool Other>
    inline Iterator<Constness>::Iterator(const Iterator<Other> &other)
      : accessor(other.accessor)
    {}


    template <bool Constness>
    inline Iterator<Constness> &
    Iterator<Constness>::operator++()
    {
      Assert(accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

      ++accessor.a_index;

      // If at end of line: do one
      // step, then cycle until we
      // find a row with a nonzero
      // number of entries.
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while ((accessor.a_row < accessor.matrix->m()) &&
                 ((accessor.matrix->in_local_range(accessor.a_row) == false) ||
                  (accessor.matrix->row_length(accessor.a_row) == 0)))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }


    template <bool Constness>
    inline Iterator<Constness>
    Iterator<Constness>::operator++(int)
    {
      const Iterator<Constness> old_state = *this;
      ++(*this);
      return old_state;
    }



    template <bool Constness>
    inline const Accessor<Constness> &Iterator<Constness>::operator*() const
    {
      return accessor;
    }



    template <bool Constness>
    inline const Accessor<Constness> *Iterator<Constness>::operator->() const
    {
      return &accessor;
    }



    template <bool Constness>
    inline bool
    Iterator<Constness>::operator==(const Iterator<Constness> &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }



    template <bool Constness>
    inline bool
    Iterator<Constness>::operator!=(const Iterator<Constness> &other) const
    {
      return !(*this == other);
    }



    template <bool Constness>
    inline bool
    Iterator<Constness>::operator<(const Iterator<Constness> &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }


    template <bool Constness>
    inline bool
    Iterator<Constness>::operator>(const Iterator<Constness> &other) const
    {
      return (other < *this);
    }

  } // namespace SparseMatrixIterators



  inline SparseMatrix::const_iterator
  SparseMatrix::begin() const
  {
    return begin(0);
  }



  inline SparseMatrix::const_iterator
  SparseMatrix::end() const
  {
    return const_iterator(this, m(), 0);
  }



  inline SparseMatrix::const_iterator
  SparseMatrix::begin(const size_type r) const
  {
    AssertIndexRange(r, m());
    if (in_local_range(r) && (row_length(r) > 0))
      return const_iterator(this, r, 0);
    else
      return end(r);
  }



  inline SparseMatrix::const_iterator
  SparseMatrix::end(const size_type r) const
  {
    AssertIndexRange(r, m());

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i = r + 1; i < m(); ++i)
      if (in_local_range(i) && (row_length(i) > 0))
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline SparseMatrix::iterator
  SparseMatrix::begin()
  {
    return begin(0);
  }



  inline SparseMatrix::iterator
  SparseMatrix::end()
  {
    return iterator(this, m(), 0);
  }



  inline SparseMatrix::iterator
  SparseMatrix::begin(const size_type r)
  {
    AssertIndexRange(r, m());
    if (in_local_range(r) && (row_length(r) > 0))
      return iterator(this, r, 0);
    else
      return end(r);
  }



  inline SparseMatrix::iterator
  SparseMatrix::end(const size_type r)
  {
    AssertIndexRange(r, m());

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i = r + 1; i < m(); ++i)
      if (in_local_range(i) && (row_length(i) > 0))
        return iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline bool
  SparseMatrix::in_local_range(const size_type index) const
  {
    TrilinosWrappers::types::int_type begin, end;
#      ifndef DEAL_II_WITH_64BIT_INDICES
    begin = matrix->RowMap().MinMyGID();
    end   = matrix->RowMap().MaxMyGID() + 1;
#      else
    begin = matrix->RowMap().MinMyGID64();
    end   = matrix->RowMap().MaxMyGID64() + 1;
#      endif

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline bool
  SparseMatrix::is_compressed() const
  {
    return compressed;
  }



  // Inline the set() and add() functions, since they will be called
  // frequently, and the compiler can optimize away some unnecessary loops
  // when the sizes are given at compile time.
  template <>
  void
  SparseMatrix::set<TrilinosScalar>(const size_type       row,
                                    const size_type       n_cols,
                                    const size_type *     col_indices,
                                    const TrilinosScalar *values,
                                    const bool            elide_zero_values);



  template <typename Number>
  void
  SparseMatrix::set(const size_type  row,
                    const size_type  n_cols,
                    const size_type *col_indices,
                    const Number *   values,
                    const bool       elide_zero_values)
  {
    std::vector<TrilinosScalar> trilinos_values(n_cols);
    std::copy(values, values + n_cols, trilinos_values.begin());
    this->set(
      row, n_cols, col_indices, trilinos_values.data(), elide_zero_values);
  }



  inline void
  SparseMatrix::set(const size_type      i,
                    const size_type      j,
                    const TrilinosScalar value)
  {
    AssertIsFinite(value);

    set(i, 1, &j, &value, false);
  }



  inline void
  SparseMatrix::set(const std::vector<size_type> &    indices,
                    const FullMatrix<TrilinosScalar> &values,
                    const bool                        elide_zero_values)
  {
    Assert(indices.size() == values.m(),
           ExcDimensionMismatch(indices.size(), values.m()));
    Assert(values.m() == values.n(), ExcNotQuadratic());

    for (size_type i = 0; i < indices.size(); ++i)
      set(indices[i],
          indices.size(),
          indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  SparseMatrix::add(const size_type      i,
                    const size_type      j,
                    const TrilinosScalar value)
  {
    AssertIsFinite(value);

    if (value == 0)
      {
        // we have to check after Insert/Add in any case to be consistent
        // with the MPI communication model, but we can save some
        // work if the addend is zero. However, these actions are done in case
        // we pass on to the other function.

        // TODO: fix this (do not run compress here, but fail)
        if (last_action == Insert)
          {
            int ierr;
            ierr = matrix->GlobalAssemble(*column_space_map,
                                          matrix->RowMap(),
                                          false);

            Assert(ierr == 0, ExcTrilinosError(ierr));
            (void)ierr; // removes -Wunused-but-set-variable in optimized mode
          }

        last_action = Add;

        return;
      }
    else
      add(i, 1, &j, &value, false);
  }



  // inline "simple" functions that are called frequently and do only involve
  // a call to some Trilinos function.
  inline SparseMatrix::size_type
  SparseMatrix::m() const
  {
#      ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalRows();
#      else
    return matrix->NumGlobalRows64();
#      endif
  }



  inline SparseMatrix::size_type
  SparseMatrix::n() const
  {
    // If the matrix structure has not been fixed (i.e., we did not have a
    // sparsity pattern), it does not know about the number of columns so we
    // must always take this from the additional column space map
    Assert(column_space_map.get() != nullptr, ExcInternalError());
    return n_global_elements(*column_space_map);
  }



  inline unsigned int
  SparseMatrix::local_size() const
  {
    return matrix->NumMyRows();
  }



  inline std::pair<SparseMatrix::size_type, SparseMatrix::size_type>
  SparseMatrix::local_range() const
  {
    size_type begin, end;
#      ifndef DEAL_II_WITH_64BIT_INDICES
    begin = matrix->RowMap().MinMyGID();
    end   = matrix->RowMap().MaxMyGID() + 1;
#      else
    begin = matrix->RowMap().MinMyGID64();
    end   = matrix->RowMap().MaxMyGID64() + 1;
#      endif

    return std::make_pair(begin, end);
  }



  inline SparseMatrix::size_type
  SparseMatrix::n_nonzero_elements() const
  {
#      ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalNonzeros();
#      else
    return matrix->NumGlobalNonzeros64();
#      endif
  }



  template <typename SparsityPatternType>
  inline void
  SparseMatrix::reinit(const IndexSet &           parallel_partitioning,
                       const SparsityPatternType &sparsity_pattern,
                       const MPI_Comm &           communicator,
                       const bool                 exchange_data)
  {
    reinit(parallel_partitioning,
           parallel_partitioning,
           sparsity_pattern,
           communicator,
           exchange_data);
  }



  template <typename number>
  inline void
  SparseMatrix::reinit(const IndexSet &parallel_partitioning,
                       const ::dealii::SparseMatrix<number> &sparse_matrix,
                       const MPI_Comm &                      communicator,
                       const double                          drop_tolerance,
                       const bool                            copy_values,
                       const ::dealii::SparsityPattern *     use_this_sparsity)
  {
    Epetra_Map map =
      parallel_partitioning.make_trilinos_map(communicator, false);
    reinit(parallel_partitioning,
           parallel_partitioning,
           sparse_matrix,
           drop_tolerance,
           copy_values,
           use_this_sparsity);
  }



  inline const Epetra_CrsMatrix &
  SparseMatrix::trilinos_matrix() const
  {
    return static_cast<const Epetra_CrsMatrix &>(*matrix);
  }



  inline const Epetra_CrsGraph &
  SparseMatrix::trilinos_sparsity_pattern() const
  {
    return matrix->Graph();
  }



  inline IndexSet
  SparseMatrix::locally_owned_domain_indices() const
  {
    return IndexSet(matrix->DomainMap());
  }



  inline IndexSet
  SparseMatrix::locally_owned_range_indices() const
  {
    return IndexSet(matrix->RangeMap());
  }



  inline void
  SparseMatrix::prepare_add()
  {
    // nothing to do here
  }



  inline void
  SparseMatrix::prepare_set()
  {
    // nothing to do here
  }


  namespace internal
  {
    namespace LinearOperatorImplementation
    {
      template <typename Solver, typename Preconditioner>
      typename std::enable_if<
        std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
          std::is_base_of<TrilinosWrappers::PreconditionBase,
                          Preconditioner>::value,
        TrilinosPayload>::type
      TrilinosPayload::inverse_payload(
        Solver &              solver,
        const Preconditioner &preconditioner) const
      {
        const auto &payload = *this;

        TrilinosPayload return_op(payload);

        // Capture by copy so the payloads are always valid

        return_op.inv_vmult = [payload, &solver, &preconditioner](
                                TrilinosPayload::Domain &     tril_dst,
                                const TrilinosPayload::Range &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(payload,
                                              tril_src,
                                              tril_dst,
                                              !payload.UseTranspose());
          solver.solve(payload, tril_dst, tril_src, preconditioner);
        };

        return_op.inv_Tvmult = [payload, &solver, &preconditioner](
                                 TrilinosPayload::Range &       tril_dst,
                                 const TrilinosPayload::Domain &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(payload,
                                              tril_src,
                                              tril_dst,
                                              payload.UseTranspose());

          const_cast<TrilinosPayload &>(payload).transpose();
          solver.solve(payload, tril_dst, tril_src, preconditioner);
          const_cast<TrilinosPayload &>(payload).transpose();
        };

        // If the input operator is already setup for transpose operations, then
        // we must do similar with its inverse.
        if (return_op.UseTranspose() == true)
          std::swap(return_op.inv_vmult, return_op.inv_Tvmult);

        return return_op;
      }

      template <typename Solver, typename Preconditioner>
      typename std::enable_if<
        !(std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
          std::is_base_of<TrilinosWrappers::PreconditionBase,
                          Preconditioner>::value),
        TrilinosPayload>::type
      TrilinosPayload::inverse_payload(Solver &, const Preconditioner &) const
      {
        TrilinosPayload return_op(*this);

        return_op.inv_vmult = [](TrilinosPayload::Domain &,
                                 const TrilinosPayload::Range &) {
          AssertThrow(false,
                      ExcMessage("Payload inv_vmult disabled because of "
                                 "incompatible solver/preconditioner choice."));
        };

        return_op.inv_Tvmult = [](TrilinosPayload::Range &,
                                  const TrilinosPayload::Domain &) {
          AssertThrow(false,
                      ExcMessage("Payload inv_vmult disabled because of "
                                 "incompatible solver/preconditioner choice."));
        };

        return return_op;
      }
    } // namespace LinearOperatorImplementation
  }   // namespace internal

  template <>
  void
  SparseMatrix::set<TrilinosScalar>(const size_type       row,
                                    const size_type       n_cols,
                                    const size_type *     col_indices,
                                    const TrilinosScalar *values,
                                    const bool            elide_zero_values);
#    endif // DOXYGEN

} /* namespace TrilinosWrappers */


DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_TRILINOS


/*-----------------------   trilinos_sparse_matrix.h     --------------------*/

#endif
/*-----------------------   trilinos_sparse_matrix.h     --------------------*/
