// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_sparse_matrix_h
#define dealii_trilinos_tpetra_sparse_matrix_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/enable_observer_pointer.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>
#  include <deal.II/lac/vector.h>

// Tpetra includes
#  include <Tpetra_Core.hpp>
#  include <Tpetra_CrsMatrix.hpp>

#  include <type_traits>


DEAL_II_NAMESPACE_OPEN

// forward declarations

template <typename Number>
class SparseMatrix;

#  ifndef DOXYGEN
template <typename MatrixType>
class BlockMatrixBase;

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename MemorySpace>
    class SparsityPattern;

    namespace SparseMatrixIterators
    {
      template <typename Number, typename MemorySpace, bool Constness>
      class Iterator;
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif

namespace LinearAlgebra
{

  namespace TpetraWrappers
  {
    /**
     * This class implements a wrapper to use the Trilinos distributed sparse
     * matrix class
     * <a
     * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html">Tpetra::CrsMatrix</a>.
     * This is precisely the kind of matrix we deal with all the time - we
     * most likely get it from some assembly process, where also entries not
     * locally owned might need to be written and hence need to be forwarded
     * to the owner process. This class is designed to be used in a distributed
     * memory architecture with an MPI compiler on the bottom, but it works
     * equally well for serial processes. The only requirement for this class to
     * work is that Trilinos has been installed with the same compiler as is
     * used for generating deal.II.
     *
     * Moreover, this class takes an optional template argument for
     * Kokkos::Nodes, allowing the usage of different Kokkos::Nodes.
     * Kokkos allows the writing of portable applications targeting,
     * for example, CUDA, OpenMP, Serial, or Threads, as backends for
     * the execution and memory spaces. The backend is chosen by
     * choosing the corresponding Kokkos Node.
     *
     * The interface of this class is modeled after the existing SparseMatrix
     * class in deal.II. It has almost the same member functions and is often
     * exchangeable. This class is templated and can be used with different
     * scalar types. However, Trilinos need to be installed with complex support
     * for usage with complex scalar types.
     *
     * @note You need to call SparseMatrix::compress() before you actually use
     * the matrix. This calls
     * <a
     * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html#aa985b225a24d2f74602e25b38b4430af">Tpetra::fillComplete</a>
     * that compresses the storage format for sparse matrices by discarding
     * unused elements and prepares the matrix for further usage
     * (e.g., for matrix-vector products).
     * However, to continue assembling the matrix, you need to call
     * SparseMatrix::resume_fill() first. Once you finish modifying
     * the matrix, you must call SparseMatrix::compress() again.
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class SparseMatrix : public EnableObserverPointer
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
      using iterator =
        SparseMatrixIterators::Iterator<Number, MemorySpace, false>;

      /**
       * Declare an alias for the const iterator class.
       */
      using const_iterator =
        SparseMatrixIterators::Iterator<Number, MemorySpace, true>;

      /**
       * Declare an alias for the type used to store matrix elements, in analogy
       * to all the other container classes.
       */
      using value_type = Number;

      /**
       * @name Constructors and initialization.
       */
      /** @{ */
      /**
       * Default constructor. Generates an empty (zero-size) matrix.
       */
      SparseMatrix();

      /**
       * Generate a matrix from a TpetraWrappers::SparsityPattern object.
       */
      SparseMatrix(const SparsityPattern<MemorySpace> &sparsity_pattern);

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
       * Move constructor. Create a new sparse matrix by stealing the internal
       * data of the `other` object.
       */
      SparseMatrix(SparseMatrix<Number, MemorySpace> &&other) noexcept;

      /**
       * Copy constructor is deleted.
       */
      SparseMatrix(const SparseMatrix<Number, MemorySpace> &) = delete;

      /**
       * operator= is deleted.
       */
      SparseMatrix<Number, MemorySpace> &
      operator=(const SparseMatrix<Number, MemorySpace> &) = delete;

      /**
       * Move assignment operator.
       */
      SparseMatrix<Number, MemorySpace> &
      operator=(SparseMatrix<Number, MemorySpace> &&other) noexcept;

      /**
       * Destructor. Made virtual so that one can use pointers to objects of
       * this class.
       */
      virtual ~SparseMatrix() override = default;

      /**
       * This function initializes the Trilinos matrix with a deal.II sparsity
       * pattern, i.e. it makes the underlying Trilinos Tpetra::CrsMatrix know
       * the position of nonzero entries according to the sparsity pattern. This
       * function is meant for use in serial programs, where there is no need to
       * specify how the matrix is going to be distributed among different
       * processors. This function works in %parallel, too, but it is
       * recommended to manually specify the %parallel partitioning of the
       * matrix using a Tpetra::Map. When run in %parallel, it is currently
       * necessary that each processor holds the sparsity_pattern structure
       * because each processor sets its rows.
       *
       * This is a collective operation that needs to be called on all
       * processors in order to avoid a dead lock.
       */
      template <typename SparsityPatternType>
      void
      reinit(const SparsityPatternType &sparsity_pattern);

      /**
       * This function reinitializes the Trilinos sparse matrix from a
       * (possibly distributed) Trilinos sparsity pattern. It also works
       * in parallel. In that case, the partitioning of the Trilinos
       * sparsity pattern is used.
       *
       * This is a collective operation that needs to be called on all
       * processors in order to avoid a dead lock.
       */
      void
      reinit(const SparsityPattern<MemorySpace> &sparsity_pattern);
      /** @} */

      /**
       * @name Constructors and initialization using an IndexSet description
       */
      /** @{ */
      /**
       * Constructor using an IndexSet and an MPI communicator to describe the
       * %parallel partitioning. The parameter @p n_max_entries_per_row sets the
       * number of nonzero entries in each row that will be allocated. Note that
       * this number does not need to be exact, and it is even allowed that the
       * actual matrix structure has more nonzero entries than specified in the
       * constructor. However it is still advantageous to provide good estimates
       * here since this will considerably increase the performance of the
       * matrix setup. However, there is no effect in the performance of
       * matrix-vector products, since Trilinos reorganizes the matrix memory
       * prior to use (in the compress() step).
       */
      SparseMatrix(const IndexSet    &parallel_partitioning,
                   const MPI_Comm     communicator          = MPI_COMM_WORLD,
                   const unsigned int n_max_entries_per_row = 0);

      /**
       * Same as before, but now set the number of non-zero entries in each
       * matrix row separately. Since we know the number of elements in the
       * matrix exactly in this case, we can already allocate the right amount
       * of memory, which makes the creation process including the insertion of
       * nonzero elements by the respective SparseMatrix::reinit call
       * considerably faster.
       */
      SparseMatrix(const IndexSet                  &parallel_partitioning,
                   const MPI_Comm                   communicator,
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
                   const MPI_Comm  communicator          = MPI_COMM_WORLD,
                   const size_type n_max_entries_per_row = 0);

      /**
       * Same as before, but now set the number of non-zero entries in each
       * matrix row separately. Since we know the number of elements in the
       * matrix exactly in this case, we can already allocate the right amount
       * of memory, which makes the creation process including the insertion of
       * nonzero elements by the respective SparseMatrix::reinit call
       * considerably faster.
       */
      SparseMatrix(const IndexSet                  &row_parallel_partitioning,
                   const IndexSet                  &col_parallel_partitioning,
                   const MPI_Comm                   communicator,
                   const std::vector<unsigned int> &n_entries_per_row);

      /**
       * This function is initializes the Trilinos Tpetra matrix according to
       * the specified @p sparsity_pattern, and also reassigns the matrix rows to
       * different processes according to the user-supplied index set @p parallel_partitioning and
       * %parallel communicator. In programs following the style of the tutorial
       * programs, this function (and the respective call for a rectangular
       * matrix) are the natural way to initialize the matrix size, its
       * distribution among the MPI processes (if run in %parallel) as well as
       * the location of non-zero elements. Trilinos stores the sparsity pattern
       * internally, so it won't be needed any more after this call, in contrast
       * to the deal.II own object. The optional argument @p exchange_data can
       * be used for reinitialization with a sparsity pattern that is not fully
       * constructed. If the flag is not set, each
       * processor just sets the elements in the sparsity pattern that belong to
       * its rows.
       *
       * This is a collective operation that needs to be called on all
       * processors in order to avoid a dead lock.
       */
      template <typename SparsityPatternType>
      std::enable_if_t<
        !std::is_same_v<SparsityPatternType, dealii::SparseMatrix<double>>>
      reinit(const IndexSet            &parallel_partitioning,
             const SparsityPatternType &sparsity_pattern,
             const MPI_Comm             communicator  = MPI_COMM_WORLD,
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
      std::enable_if_t<
        !std::is_same_v<SparsityPatternType, dealii::SparseMatrix<double>>>
      reinit(const IndexSet            &row_parallel_partitioning,
             const IndexSet            &col_parallel_partitioning,
             const SparsityPatternType &sparsity_pattern,
             const MPI_Comm             communicator  = MPI_COMM_WORLD,
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
       * This is a @ref GlossCollectiveOperation "collective operation" that needs to be called on all
       * processors in order to avoid a dead lock.
       */
      void
      reinit(const IndexSet                     &row_parallel_partitioning,
             const IndexSet                     &col_parallel_partitioning,
             const dealii::SparseMatrix<Number> &dealii_sparse_matrix,
             const MPI_Comm                      communicator = MPI_COMM_WORLD,
             const double                        drop_tolerance    = 1e-13,
             const bool                          copy_values       = true,
             const dealii::SparsityPattern      *use_this_sparsity = nullptr);

      /** @} */

      /**
       * @name Information on the matrix
       */
      /** @{ */
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
      size_t
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
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const;
      /** @} */

      /**
       * @name Modifying entries
       */
      /** @{ */
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
       * Multiply the entire matrix by a fixed factor.
       */
      SparseMatrix &
      operator*=(const Number factor);

      /**
       * Divide the entire matrix by a fixed factor.
       */
      SparseMatrix &
      operator/=(const Number factor);

      /**
       * Copy the given (Trilinos) matrix (sparsity pattern and entries).
       */
      void
      copy_from(const SparseMatrix<Number, MemorySpace> &source);

      /**
       * Add @p value to the element (<i>i,j</i>).
       * Just as the respective call in deal.II SparseMatrix<Number,
       * MemorySpace> class. Moreover, if <tt>value</tt> is not a finite number
       * an exception is thrown.
       *
       * @note When add is called on a compressed matrix, the matrix is set
       * back to an uncompressed state.
       */
      void
      add(const size_type i, const size_type j, const Number value);

      /**
       * Add an array of values given by <tt>values</tt> in the given global
       * matrix row at columns specified by col_indices in the sparse matrix.
       * Just as the respective call in deal.II SparseMatrix<Number,
       * MemorySpace> class. The optional parameter <tt>elide_zero_values</tt>
       * can be used to specify whether zero values should be added anyway or
       * these should be filtered away and only non-zero data is added. The
       * default value is <tt>true</tt>, i.e., zero values won't be added into
       * the matrix.
       *
       * @note When add is called on a compressed matrix, the matrix is set
       * back to an uncompressed state.
       */
      void
      add(const size_type  row,
          const size_type  n_cols,
          const size_type *col_indices,
          const Number    *values,
          const bool       elide_zero_values      = true,
          const bool       col_indices_are_sorted = false);

      /**
       * Add <tt>matrix</tt> scaled by <tt>factor</tt> to this matrix, i.e. the
       * matrix <tt>factor*matrix</tt> is added to <tt>this</tt>. If the
       * sparsity pattern of the calling matrix does not contain all the
       * elements in the sparsity pattern of the input matrix, this function
       * will throw an exception.
       */
      void
      add(const Number factor, const SparseMatrix<Number, MemorySpace> &matrix);

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
       * behavior imposed by the underlying Tpetra::CrsMatrix data structure:
       * If the same matrix entry is inserted more than once, the matrix entries
       * will be added upon calling compress() (since Tpetra does not track
       * values to the same entry before the final compress() is called), even
       * if VectorOperation::insert is specified as argument to compress(). In
       * the case you cannot make sure that matrix entries are only set once,
       * initialize the matrix with a sparsity pattern to fix the matrix
       * structure before inserting elements.
       */
      void
      set(const size_type i, const size_type j, const Number value);

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
       * behavior imposed by the underlying Tpetra::CrsMatrix data structure:
       * If the same matrix entry is inserted more than once, the matrix entries
       * will be added upon calling compress() (since Epetra does not track
       * values to the same entry before the final compress() is called), even
       * if VectorOperation::insert is specified as argument to compress(). In
       * the case you cannot make sure that matrix entries are only set once,
       * initialize the matrix with a sparsity pattern to fix the matrix
       * structure before inserting elements.
       */
      void
      set(const std::vector<size_type> &indices,
          const FullMatrix<Number>     &full_matrix,
          const bool                    elide_zero_values = false);

      /**
       * Same function as before, but now including the possibility to use
       * rectangular full_matrices and different local-to-global indexing on
       * rows and columns, respectively.
       */
      void
      set(const std::vector<size_type> &row_indices,
          const std::vector<size_type> &col_indices,
          const FullMatrix<Number>     &full_matrix,
          const bool                    elide_zero_values = false);

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
       * behavior imposed by the underlying Tpetra::CrsMatrix data structure:
       * If the same matrix entry is inserted more than once, the matrix entries
       * will be added upon calling compress() (since Epetra does not track
       * values to the same entry before the final compress() is called), even
       * if VectorOperation::insert is specified as argument to compress(). In
       * the case you cannot make sure that matrix entries are only set once,
       * initialize the matrix with a sparsity pattern to fix the matrix
       * structure before inserting elements.
       */
      void
      set(const size_type               row,
          const std::vector<size_type> &col_indices,
          const std::vector<Number>    &values,
          const bool                    elide_zero_values = false);

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
       * behavior imposed by the underlying Tpetra::CrsMatrix data structure:
       * If the same matrix entry is inserted more than once, the matrix entries
       * will be added upon calling compress() (since Epetra does not track
       * values to the same entry before the final compress() is called), even
       * if VectorOperation::insert is specified as argument to compress(). In
       * the case you cannot make sure that matrix entries are only set once,
       * initialize the matrix with a sparsity pattern to fix the matrix
       * structure before inserting elements.
       */
      void
      set(const size_type  row,
          const size_type  n_cols,
          const size_type *col_indices,
          const Number    *values,
          const bool       elide_zero_values = false);

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
      clear_row(const size_type row, const Number new_diag_value = 0);

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
      clear_rows(const ArrayView<const size_type> &rows,
                 const Number                      new_diag_value = 0);

      /**
       * Release all memory and return to a state just like after having called
       * the default constructor.
       *
       * This is a @ref GlossCollectiveOperation "collective operation" that needs to be called on all
       * processors in order to avoid a dead lock.
       */
      void
      clear();

      /** @} */
      /**
       * @name Entry Access
       */
      /** @{ */

      /**
       * Return the value of the entry (<i>i,j</i>).  This may be an expensive
       * operation and you should always take care where to call this function.
       * As in the deal.II sparse matrix class, we throw an exception if the
       * respective entry doesn't exist in the sparsity pattern of this class,
       * which is requested from Trilinos. Moreover, an exception will be thrown
       * when the requested element is not saved on the calling process.
       */
      Number
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
      Number
      el(const size_type i, const size_type j) const;

      /**
       * Return the main diagonal element in the <i>i</i>th row. This function
       * throws an error if the matrix is not quadratic and it also throws an
       * error if <i>(i,i)</i> is not element of the local matrix.
       */
      Number
      diag_element(const size_type i) const;

      /** @} */
      /**
       * @name Multiplications
       */
      /** @{ */
      /*
       * Matrix-vector multiplication: let <i>dst = M*src</i> with <i>M</i>
       * being this matrix.
       *
       * Source and destination must not be the same vector.
       *
       * The vector @p dst has to be initialized with the same IndexSet that was
       * used for the row indices of the matrix and the vector @p src has to be
       * initialized with the same IndexSet that was used for the column indices
       * of the matrix.
       */
      template <typename InputVectorType>
      void
      vmult(InputVectorType &dst, const InputVectorType &src) const;

      /*
       * Matrix-vector multiplication: let <i>dst = M<sup>T</sup>*src</i> with
       * <i>M</i> being this matrix. This function does the same as vmult() but
       * takes the transposed matrix.
       *
       * Source and destination must not be the same vector.
       */
      template <typename InputVectorType>
      void
      Tvmult(InputVectorType &dst, const InputVectorType &src) const;

      /**
       * Adding matrix-vector multiplication. Add <i>M*src</i> on <i>dst</i>
       * with <i>M</i> being this matrix.
       *
       * Source and destination must not be the same vector.
       */
      template <typename InputVectorType>
      void
      vmult_add(InputVectorType &dst, const InputVectorType &src) const;

      /**
       * Adding matrix-vector multiplication. Add <i>M<sup>T</sup>*src</i> to
       * <i>dst</i> with <i>M</i> being this matrix. This function does the same
       * as vmult_add() but takes the transposed matrix.
       *
       * Source and destination must not be the same vector.
       */
      template <typename InputVectorType>
      void
      Tvmult_add(InputVectorType &dst, const InputVectorType &src) const;

      /**
       * Return the square of the norm of the vector $v$ with respect to the
       * norm induced by this matrix, i.e., $\left(v,Mv\right)$. This is useful,
       * e.g. in the finite element context, where the $L_2$ norm of a function
       * equals the matrix norm with respect to the @ref GlossMassMatrix "mass matrix" of the vector
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
      Number
      matrix_norm_square(const Vector<Number, MemorySpace> &v) const;

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
      Number
      matrix_scalar_product(const Vector<Number, MemorySpace> &u,
                            const Vector<Number, MemorySpace> &v) const;

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
      Number
      residual(Vector<Number, MemorySpace>       &dst,
               const Vector<Number, MemorySpace> &x,
               const Vector<Number, MemorySpace> &b) const;

      /** @} */

      /**
       * @name Matrix norms
       */
      /** @{ */

      /**
       * Return the frobenius norm of the matrix, i.e. the square root of the
       * sum of squares of all entries in the matrix.
       */
      Number
      frobenius_norm() const;

      /** @} */

      /**
       * @name Mixed Stuff
       */
      /** @{ */
      /**
       * Print the matrix to the given stream, using the format (line,col)
       * value, i.e. one nonzero entry of the matrix per line. The optional flag
       * outputs the sparsity pattern in Trilinos style, where the data is
       * sorted according to the processor number when printed to the stream, as
       * well as a summary of the matrix like the global size.
       */
      void
      print(std::ostream &out,
            const bool    print_detailed_trilinos_information = false) const;

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
       *
       * @note The @p operation can be safely omitted, as that parameter is not
       * used at all and is only present to ensure compatibility with other
       * SparseMatrix classes.
       */
      void
      compress(VectorOperation::values operation);

      /**
       * This function must be called to allow for changes to the structure
       * of the matrix again after compress() was called.
       * Once you are done modifying the matrix structure, you must call
       * compress() again.
       */
      void
      resume_fill();

      /**
       * Return a const reference to the underlying Trilinos
       * <a
       * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html">Tpetra::CrsMatrix</a>
       * class.
       */
      const TpetraTypes::MatrixType<Number, MemorySpace> &
      trilinos_matrix() const;

      /**
       * Return a (modifiable) reference to the underlying Trilinos
       * <a
       * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html">Tpetra::CrsMatrix</a>
       * class.
       */
      TpetraTypes::MatrixType<Number, MemorySpace> &
      trilinos_matrix();

      /**
       * Return a const
       * <a
       * href="https://docs.trilinos.org/dev/packages/teuchos/doc/html/classTeuchos_1_1RCP.html">Teuchos::RCP</a>
       * to the underlying Trilinos
       * <a
       * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html">Tpetra::CrsMatrix</a>
       * class.
       */
      Teuchos::RCP<const TpetraTypes::MatrixType<Number, MemorySpace>>
      trilinos_rcp() const;

      /**
       * Return a (modifiable)
       * <a
       * href="https://docs.trilinos.org/dev/packages/teuchos/doc/html/classTeuchos_1_1RCP.html">Teuchos::RCP</a>
       * to the underlying Trilinos
       * <a
       * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html">Tpetra::CrsMatrix</a>
       * class.
       */
      Teuchos::RCP<TpetraTypes::MatrixType<Number, MemorySpace>>
      trilinos_rcp();
      /** @} */

      /**
       * @name Partitioners
       */
      /** @{ */

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

      /** @} */

      /**
       * @name Iterators
       */
      /** @{ */

      /**
       * Return an iterator pointing to the first element of the matrix.
       *
       * The elements accessed by iterators within each row are ordered in the
       * way in which Trilinos stores them, though the implementation guarantees
       * that all elements of one row are accessed before the elements of the
       * next row. If your algorithm relies on visiting elements within one row,
       * you will need to consult with the Trilinos documentation on the order
       * in which it stores data. It is, however, generally not a good and
       * long-term stable idea to rely on the order in which receive elements if
       * you iterate over them.
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
       * in which it stores data. It is, however, generally not a good and
       * long-term stable idea to rely on the order in which receive elements if
       * you iterate over them.
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

      /** @} */
      /**
       * @addtogroup Exceptions
       */
      /** @{ */

      /**
       * Exception
       */
      DeclException0(ExcMatrixNotCompressed);

      /**
       * Exception
       */
      DeclExceptionMsg(
        ExcSourceEqualsDestination,
        "You are attempting an operation on two vectors that "
        "are the same object, but the operation requires that the "
        "two objects are in fact different.");

      /*
       * Exception
       */
      DeclExceptionMsg(ExcColMapMismatch,
                       "The column partitioning of a matrix does not match "
                       "the partitioning of a vector you are trying to "
                       "multiply it with. Are you multiplying the "
                       "matrix with a vector that has ghost elements?");

      /*
       * Exception
       */
      DeclExceptionMsg(ExcDomainMapMismatch,
                       "The row partitioning of a matrix does not match "
                       "the partitioning of a vector you are trying to "
                       "put the result of a matrix-vector product in. "
                       "Are you trying to put the product of the "
                       "matrix with a vector into a vector that has "
                       "ghost elements?");

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
      DeclException4(ExcAccessToNonLocalElement,
                     size_type,
                     size_type,
                     size_type,
                     size_type,
                     << "You tried to access element (" << arg1 << '/' << arg2
                     << ')'
                     << " of a distributed matrix, but only rows in range ["
                     << arg3 << ',' << arg4
                     << "] are stored locally and can be accessed.");

      /** @} */

    private:
      /**
       * Helper function for operator() and el().
       */
      Number
      element(const size_type i, const size_type j, const bool no_error) const;

      /**
       * Pointer to the user-supplied Tpetra Trilinos mapping of the matrix
       * columns that assigns parts of the matrix to the individual processes.
       *
       * @note The Trilinos matrix is row-oriented, and the row_space_map is
       * therefore stored in the Trilinos matrix itself. The additional
       * information from the column space map is used to speed up the
       * assembly process.
       */
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> column_space_map;

      /**
       * A sparse matrix object in Trilinos to be used for finite element based
       * problems which allows for assembling into non-local elements.  The
       * actual type, a sparse matrix, is set in the constructor.
       */
      Teuchos::RCP<TpetraTypes::MatrixType<Number, MemorySpace>> matrix;

      /**
       * A boolean variable to hold information on whether the matrix is
       * fill complete or if the matrix is in compute mode.
       */
      bool compressed;

      /**
       * For some matrix storage formats, in particular for the PETSc
       * distributed blockmatrices, set and add operations on individual
       * elements can not be freely mixed. Rather, one has to synchronize
       * operations when one wants to switch from setting elements to adding to
       * elements.  BlockMatrixBase automatically synchronizes the access by
       * calling this helper function for each block.  This function ensures
       * that the matrix is in a state that allows adding elements; if it
       * previously already was in this state, the function does nothing.
       *
       * This function is called from BlockMatrixBase.
       */
      void
      prepare_add();

      /**
       * Same as prepare_add() but prepare the matrix for setting elements if
       * the representation of elements in this class requires such an
       * operation.
       *
       * This function is called from BlockMatrixBase.
       */
      void
      prepare_set();

      // To allow calling protected prepare_add() and prepare_set().
      friend class BlockMatrixBase<SparseMatrix<Number, MemorySpace>>;
    }; // class SparseMatrix

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
       * Handling of indices for both constant and non constant Accessor objects
       *
       * For a regular dealii::SparseMatrix, we would use an accessor for the
       * sparsity pattern. For Trilinos matrices, this does not seem so simple,
       * therefore, we write a little base class here.
       */
      template <typename Number, typename MemorySpace>
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
        AccessorBase(SparseMatrix<Number, MemorySpace> *matrix,
                     const size_type                    row,
                     const size_type                    index);

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
        mutable SparseMatrix<Number, MemorySpace> *matrix;
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
        std::shared_ptr<std::vector<dealii::types::signed_global_dof_index>>
          colnum_cache;

        /**
         * Cache for the values of this row.
         */
        std::shared_ptr<std::vector<Number>> value_cache;

      private:
        friend class Iterator<Number, MemorySpace, false>;
        friend class Iterator<Number, MemorySpace, true>;
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
      template <typename Number, typename MemorySpace, bool Constness>
      class Accessor : public AccessorBase<Number, MemorySpace>
      {
        /**
         * Value of this matrix entry.
         */
        Number
        value() const;

        /**
         * Value of this matrix entry.
         */
        Number &
        value();
      };

      /**
       * The specialization for a const Accessor.
       */
      template <typename Number, typename MemorySpace>
      class Accessor<Number, MemorySpace, true>
        : public AccessorBase<Number, MemorySpace>
      {
      public:
        /**
         * Typedef for the type (including constness) of the matrix to be used
         * here.
         */
        using MatrixType = const SparseMatrix<Number, MemorySpace>;

        /**
         * Typedef for the size type of the matrix to be used here.
         */
        using size_type = typename AccessorBase<Number, MemorySpace>::size_type;

        /**
         * Constructor. Since we use accessors only for read access, a const
         * matrix pointer is sufficient.
         */
        Accessor(MatrixType     *matrix,
                 const size_type row,
                 const size_type index);

        /**
         * Copy constructor to get from a const or non-const accessor to a const
         * accessor.
         */
        template <bool Other>
        Accessor(const Accessor<Number, MemorySpace, Other> &a);

        /**
         * Value of this matrix entry.
         */
        Number
        value() const;
      };

      /**
       * The specialization for a mutable Accessor.
       */
      template <typename Number, typename MemorySpace>
      class Accessor<Number, MemorySpace, false>
        : public AccessorBase<Number, MemorySpace>
      {
        class Reference
        {
        public:
          /**
           * Constructor.
           */
          Reference(const Accessor<Number, MemorySpace, false> &accessor);

          /**
           * Conversion operator to the data type of the matrix.
           */
          operator Number() const;

          /**
           * Set the element of the matrix we presently point to to @p n.
           */
          const Reference &
          operator=(const Number n) const;

          /**
           * Add @p n to the element of the matrix we presently point to.
           */
          const Reference &
          operator+=(const Number n) const;

          /**
           * Subtract @p n from the element of the matrix we presently point to.
           */
          const Reference &
          operator-=(const Number n) const;

          /**
           * Multiply the element of the matrix we presently point to by @p n.
           */
          const Reference &
          operator*=(const Number n) const;

          /**
           * Divide the element of the matrix we presently point to by @p n.
           */
          const Reference &
          operator/=(const Number n) const;

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
        using MatrixType = SparseMatrix<Number, MemorySpace>;

        /**
         * Typedef for the size type of the matrix to be used here.
         */
        using size_type = typename AccessorBase<Number, MemorySpace>::size_type;

        /**
         * Constructor. Since we use accessors only for read access, a const
         * matrix pointer is sufficient.
         */
        Accessor(MatrixType     *matrix,
                 const size_type row,
                 const size_type index);

        /**
         * Value of this matrix entry.
         */
        Reference
        value() const;

      private:
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
       * @ingroup TpetraWrappers
       */
      template <typename Number, typename MemorySpace, bool Constness>
      class Iterator
      {
      public:
        /**
         * Declare type for container size.
         */
        using size_type = dealii::types::global_dof_index;

        /**
         * A type that denotes what data types is used to express the difference
         * between two iterators.
         */
        using difference_type = dealii::types::global_dof_index;

        /**
         * An alias for the type you get when you dereference an iterator of the
         * current kind.
         */
        using value_type = Number;

        /**
         * Typedef for the matrix type (including constness) we are to operate
         * on.
         */
        using MatrixType =
          typename Accessor<Number, MemorySpace, Constness>::MatrixType;

        /**
         * Constructor. Create an iterator into the matrix @p matrix for the
         * given row and the index within it.
         */
        Iterator(MatrixType     *matrix,
                 const size_type row,
                 const size_type index);

        /**
         * Copy constructor with optional change of constness.
         */
        template <bool Other>
        Iterator(const Iterator<Number, MemorySpace, Other> &other);

        /**
         * Prefix increment.
         */
        Iterator<Number, MemorySpace, Constness> &
        operator++();

        /**
         * Postfix increment.
         */
        Iterator<Number, MemorySpace, Constness>
        operator++(int);

        /**
         * Dereferencing operator.
         */
        const Accessor<Number, MemorySpace, Constness> &
        operator*() const;

        /**
         * Dereferencing operator.
         */
        const Accessor<Number, MemorySpace, Constness> *
        operator->() const;

        /**
         * Comparison. True, if both iterators point to the same matrix
         * position.
         */
        template <bool OtherConstness>
        bool
        operator==(const Iterator<Number, MemorySpace, OtherConstness> &) const;

        /**
         * Inverse of <tt>==</tt>.
         */
        template <bool OtherConstness>
        bool
        operator!=(const Iterator<Number, MemorySpace, OtherConstness> &) const;

        /**
         * Comparison operator. Result is true if either the first row number is
         * smaller or if the row numbers are equal and the first index is
         * smaller.
         */
        template <bool OtherConstness>
        bool
        operator<(const Iterator<Number, MemorySpace, OtherConstness> &) const;

        /**
         * Comparison operator. The opposite of the previous operator
         */
        template <bool OtherConstness>
        bool
        operator>(const Iterator<Number, MemorySpace, OtherConstness> &) const;

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
        Accessor<Number, MemorySpace, Constness> accessor;

        friend class Iterator<Number, MemorySpace, true>;
        friend class Iterator<Number, MemorySpace, false>;
      };
    } // namespace SparseMatrixIterators

  } // namespace TpetraWrappers

} // namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

namespace std
{
  template <typename Number, typename MemorySpace, bool Constness>
  struct iterator_traits<
    dealii::LinearAlgebra::TpetraWrappers::SparseMatrixIterators::
      Iterator<Number, MemorySpace, Constness>>
  {
    using iterator_category = forward_iterator_tag;
    using value_type =
      typename dealii::LinearAlgebra::TpetraWrappers::SparseMatrixIterators::
        Iterator<Number, MemorySpace, Constness>::value_type;
    using difference_type =
      typename dealii::LinearAlgebra::TpetraWrappers::SparseMatrixIterators::
        Iterator<Number, MemorySpace, Constness>::difference_type;
  };
} // namespace std

/* ------------------------- Inline functions ---------------------- */


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{

  namespace TpetraWrappers
  {
    template <typename Number, typename MemorySpace>
    inline void
    SparseMatrix<Number, MemorySpace>::set(const size_type i,
                                           const size_type j,
                                           const Number    value)
    {
      set(i, 1, &j, &value, false);
    }



    template <typename Number, typename MemorySpace>
    inline void
    SparseMatrix<Number, MemorySpace>::add(const size_type i,
                                           const size_type j,
                                           const Number    value)
    {
      add(i, 1, &j, &value, false);
    }



    template <typename Number, typename MemorySpace>
    inline Number
    SparseMatrix<Number, MemorySpace>::residual(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &x,
      const Vector<Number, MemorySpace> &b) const
    {
      vmult(dst, x);
      dst -= b;
      dst *= -1.;

      return dst.l2_norm();
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::size_type
    SparseMatrix<Number, MemorySpace>::m() const
    {
      return matrix->getRowMap()->getGlobalNumElements();
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::size_type
    SparseMatrix<Number, MemorySpace>::n() const
    {
      // If the matrix structure has not been fixed (i.e., we did not have a
      // sparsity pattern), it does not know about the number of columns, so we
      // must always take this from the additional column space map
      Assert(column_space_map.get() != nullptr, ExcInternalError());
      return column_space_map->getGlobalNumElements();
    }



    template <typename Number, typename MemorySpace>
    inline bool
    SparseMatrix<Number, MemorySpace>::is_compressed() const
    {
      return compressed;
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::const_iterator
    SparseMatrix<Number, MemorySpace>::begin() const
    {
      return begin(0);
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::const_iterator
    SparseMatrix<Number, MemorySpace>::end() const
    {
      return const_iterator(this, m(), 0);
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::const_iterator
    SparseMatrix<Number, MemorySpace>::begin(const size_type r) const
    {
      AssertIndexRange(r, m());
      if (in_local_range(r) && (row_length(r) > 0))
        return const_iterator(this, r, 0);
      else
        return end(r);
    }


    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::const_iterator
    SparseMatrix<Number, MemorySpace>::end(const size_type r) const
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



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::iterator
    SparseMatrix<Number, MemorySpace>::begin()
    {
      return begin(0);
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::iterator
    SparseMatrix<Number, MemorySpace>::end()
    {
      return iterator(this, m(), 0);
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::iterator
    SparseMatrix<Number, MemorySpace>::begin(const size_type r)
    {
      AssertIndexRange(r, m());
      if (in_local_range(r) && (row_length(r) > 0))
        return iterator(this, r, 0);
      else
        return end(r);
    }



    template <typename Number, typename MemorySpace>
    inline typename SparseMatrix<Number, MemorySpace>::iterator
    SparseMatrix<Number, MemorySpace>::end(const size_type r)
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



    template <typename Number, typename MemorySpace>
    inline bool
    SparseMatrix<Number, MemorySpace>::in_local_range(
      const size_type index) const
    {
      const size_type begin = matrix->getRowMap()->getMinGlobalIndex();
      const size_type end   = matrix->getRowMap()->getMaxGlobalIndex() + 1;

      return ((index >= begin) && (index < end));
    }



    template <typename Number, typename MemorySpace>
    unsigned int
    SparseMatrix<Number, MemorySpace>::row_length(const size_type row) const
    {
      auto n_entries = matrix->getNumEntriesInGlobalRow(row);
      Assert(n_entries !=
               Teuchos::OrdinalTraits<decltype(n_entries)>::invalid(),
             ExcAccessToNonlocalRow(row));

      return n_entries;
    }



    template <typename Number, typename MemorySpace>
    inline void
    SparseMatrix<Number, MemorySpace>::prepare_add()
    {
      // nothing to do here
    }



    template <typename Number, typename MemorySpace>
    inline void
    SparseMatrix<Number, MemorySpace>::prepare_set()
    {
      // nothing to do here
    }



    template <typename Number, typename MemorySpace>
    inline const TpetraTypes::MatrixType<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::trilinos_matrix() const
    {
      return *matrix;
    }



    template <typename Number, typename MemorySpace>
    inline TpetraTypes::MatrixType<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::trilinos_matrix()
    {
      return *matrix;
    }



    template <typename Number, typename MemorySpace>
    inline Teuchos::RCP<const TpetraTypes::MatrixType<Number, MemorySpace>>
    SparseMatrix<Number, MemorySpace>::trilinos_rcp() const
    {
      return matrix.getConst();
    }



    template <typename Number, typename MemorySpace>
    inline Teuchos::RCP<TpetraTypes::MatrixType<Number, MemorySpace>>
    SparseMatrix<Number, MemorySpace>::trilinos_rcp()
    {
      return matrix;
    }



    template <typename Number, typename MemorySpace>
    inline IndexSet
    SparseMatrix<Number, MemorySpace>::locally_owned_domain_indices() const
    {
      return IndexSet(matrix->getDomainMap());
    }



    template <typename Number, typename MemorySpace>
    inline IndexSet
    SparseMatrix<Number, MemorySpace>::locally_owned_range_indices() const
    {
      return IndexSet(matrix->getRangeMap());
    }


    namespace SparseMatrixIterators
    {
      template <typename Number, typename MemorySpace>
      inline AccessorBase<Number, MemorySpace>::AccessorBase(
        SparseMatrix<Number, MemorySpace> *matrix,
        size_type                          row,
        size_type                          index)
        : matrix(matrix)
        , a_row(row)
        , a_index(index)
      {
        visit_present_row();
      }


      template <typename Number, typename MemorySpace>
      inline typename AccessorBase<Number, MemorySpace>::size_type
      AccessorBase<Number, MemorySpace>::row() const
      {
        Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
        return a_row;
      }


      template <typename Number, typename MemorySpace>
      inline typename AccessorBase<Number, MemorySpace>::size_type
      AccessorBase<Number, MemorySpace>::column() const
      {
        Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
        return (*colnum_cache)[a_index];
      }


      template <typename Number, typename MemorySpace>
      inline typename AccessorBase<Number, MemorySpace>::size_type
      AccessorBase<Number, MemorySpace>::index() const
      {
        Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
        return a_index;
      }



      template <typename Number, typename MemorySpace>
      void
      AccessorBase<Number, MemorySpace>::visit_present_row()
      {
        // if we are asked to visit the past-the-end line, then simply
        // release all our caches and go on with life.
        //
        // do the same if the row we're supposed to visit is not locally
        // owned. this is simply going to make non-locally owned rows
        // look like they're empty
        if ((this->a_row == matrix->m()) ||
            (matrix->in_local_range(this->a_row) == false))
          {
            colnum_cache.reset();
            value_cache.reset();

            return;
          }

        // get a representation of the present row
        size_t                            ncols;
        TrilinosWrappers::types::int_type colnums =
          matrix->row_length(this->a_row);
        if (value_cache.get() == nullptr)
          {
            value_cache  = std::make_shared<std::vector<Number>>(colnums);
            colnum_cache = std::make_shared<
              std::vector<dealii::types::signed_global_dof_index>>(colnums);
          }
        else
          {
            value_cache->resize(colnums);
            colnum_cache->resize(colnums);
          }

        typename TpetraTypes::MatrixType<Number, MemorySpace>::
          nonconst_global_inds_host_view_type col_indices(colnum_cache->data(),
                                                          colnums);
        typename TpetraTypes::MatrixType<Number, MemorySpace>::
          nonconst_values_host_view_type values(value_cache->data(), colnums);

        matrix->trilinos_matrix().getGlobalRowCopy(this->a_row,
                                                   col_indices,
                                                   values,
                                                   ncols);

        AssertDimension(ncols, colnums);

        // copy it into our caches if the
        // line isn't empty. if it is, then
        // we've done something wrong, since
        // we shouldn't have initialized an
        // iterator for an empty line (what
        // would it point to?)
      }



      template <typename Number, typename MemorySpace>
      inline Accessor<Number, MemorySpace, true>::Accessor(
        MatrixType     *matrix,
        const size_type row,
        const size_type index)
        : AccessorBase<Number, MemorySpace>(
            const_cast<SparseMatrix<Number, MemorySpace> *>(matrix),
            row,
            index)
      {}



      template <typename Number, typename MemorySpace>
      template <bool Other>
      inline Accessor<Number, MemorySpace, true>::Accessor(
        const Accessor<Number, MemorySpace, Other> &other)
        : AccessorBase<Number, MemorySpace>(other)
      {}



      template <typename Number, typename MemorySpace>
      inline Number
      Accessor<Number, MemorySpace, true>::value() const
      {
        Assert((AccessorBase<Number, MemorySpace>::a_row <
                AccessorBase<Number, MemorySpace>::matrix->m()),
               ExcBeyondEndOfMatrix());
        return (*AccessorBase<Number, MemorySpace>::value_cache)
          [AccessorBase<Number, MemorySpace>::a_index];
      }



      template <typename Number, typename MemorySpace>
      inline Accessor<Number, MemorySpace, false>::Reference::Reference(
        const Accessor<Number, MemorySpace, false> &acc)
        : accessor(const_cast<Accessor<Number, MemorySpace, false> &>(acc))
      {}



      template <typename Number, typename MemorySpace>
      inline Accessor<Number, MemorySpace, false>::Reference::operator Number()
        const
      {
        return (*accessor.value_cache)[accessor.a_index];
      }



      template <typename Number, typename MemorySpace>
      inline const typename Accessor<Number, MemorySpace, false>::Reference &
      Accessor<Number, MemorySpace, false>::Reference::operator=(
        const Number n) const
      {
        (*accessor.value_cache)[accessor.a_index] = n;
        accessor.matrix->set(accessor.row(),
                             accessor.column(),
                             static_cast<Number>(*this));
        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const typename Accessor<Number, MemorySpace, false>::Reference &
      Accessor<Number, MemorySpace, false>::Reference::operator+=(
        const Number n) const
      {
        (*accessor.value_cache)[accessor.a_index] += n;
        accessor.matrix->set(accessor.row(),
                             accessor.column(),
                             static_cast<Number>(*this));
        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const typename Accessor<Number, MemorySpace, false>::Reference &
      Accessor<Number, MemorySpace, false>::Reference::operator-=(
        const Number n) const
      {
        (*accessor.value_cache)[accessor.a_index] -= n;
        accessor.matrix->set(accessor.row(),
                             accessor.column(),
                             static_cast<Number>(*this));
        return *this;
      }


      template <typename Number, typename MemorySpace>
      inline const typename Accessor<Number, MemorySpace, false>::Reference &
      Accessor<Number, MemorySpace, false>::Reference::operator*=(
        const Number n) const
      {
        (*accessor.value_cache)[accessor.a_index] *= n;
        accessor.matrix->set(accessor.row(),
                             accessor.column(),
                             static_cast<Number>(*this));
        return *this;
      }


      template <typename Number, typename MemorySpace>
      inline const typename Accessor<Number, MemorySpace, false>::Reference &
      Accessor<Number, MemorySpace, false>::Reference::operator/=(
        const Number n) const
      {
        (*accessor.value_cache)[accessor.a_index] /= n;
        accessor.matrix->set(accessor.row(),
                             accessor.column(),
                             static_cast<Number>(*this));
        return *this;
      }


      template <typename Number, typename MemorySpace>
      inline Accessor<Number, MemorySpace, false>::Accessor(
        MatrixType     *matrix,
        const size_type row,
        const size_type index)
        : AccessorBase<Number, MemorySpace>(matrix, row, index)
      {}


      template <typename Number, typename MemorySpace>
      inline typename Accessor<Number, MemorySpace, false>::Reference
      Accessor<Number, MemorySpace, false>::value() const
      {
        Assert((AccessorBase<Number, MemorySpace>::a_row <
                AccessorBase<Number, MemorySpace>::matrix->m()),
               ExcBeyondEndOfMatrix());
        return {*this};
      }



      template <typename Number, typename MemorySpace, bool Constness>
      inline Iterator<Number, MemorySpace, Constness>::Iterator(
        MatrixType     *matrix,
        const size_type row,
        const size_type index)
        : accessor(matrix, row, index)
      {}


      template <typename Number, typename MemorySpace, bool Constness>
      template <bool Other>
      inline Iterator<Number, MemorySpace, Constness>::Iterator(
        const Iterator<Number, MemorySpace, Other> &other)
        : accessor(other.accessor)
      {}



      template <typename Number, typename MemorySpace, bool Constness>
      inline Iterator<Number, MemorySpace, Constness> &
      Iterator<Number, MemorySpace, Constness>::operator++()
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

            while (
              (accessor.a_row < accessor.matrix->m()) &&
              ((accessor.matrix->in_local_range(accessor.a_row) == false) ||
               (accessor.matrix->row_length(accessor.a_row) == 0)))
              ++accessor.a_row;

            accessor.visit_present_row();
          }
        return *this;
      }


      template <typename Number, typename MemorySpace, bool Constness>
      inline Iterator<Number, MemorySpace, Constness>
      Iterator<Number, MemorySpace, Constness>::operator++(int)
      {
        const Iterator<Number, MemorySpace, Constness> old_state = *this;
        ++(*this);
        return old_state;
      }



      template <typename Number, typename MemorySpace, bool Constness>
      inline const Accessor<Number, MemorySpace, Constness> &
      Iterator<Number, MemorySpace, Constness>::operator*() const
      {
        return accessor;
      }



      template <typename Number, typename MemorySpace, bool Constness>
      inline const Accessor<Number, MemorySpace, Constness> *
      Iterator<Number, MemorySpace, Constness>::operator->() const
      {
        return &accessor;
      }



      template <typename Number, typename MemorySpace, bool Constness>
      template <bool OtherConstness>
      inline bool
      Iterator<Number, MemorySpace, Constness>::operator==(
        const Iterator<Number, MemorySpace, OtherConstness> &other) const
      {
        return (accessor.a_row == other.accessor.a_row &&
                accessor.a_index == other.accessor.a_index);
      }



      template <typename Number, typename MemorySpace, bool Constness>
      template <bool OtherConstness>
      inline bool
      Iterator<Number, MemorySpace, Constness>::operator!=(
        const Iterator<Number, MemorySpace, OtherConstness> &other) const
      {
        return !(*this == other);
      }



      template <typename Number, typename MemorySpace, bool Constness>
      template <bool OtherConstness>
      inline bool
      Iterator<Number, MemorySpace, Constness>::operator<(
        const Iterator<Number, MemorySpace, OtherConstness> &other) const
      {
        return (accessor.row() < other.accessor.row() ||
                (accessor.row() == other.accessor.row() &&
                 accessor.index() < other.accessor.index()));
      }


      template <typename Number, typename MemorySpace, bool Constness>
      template <bool OtherConstness>
      inline bool
      Iterator<Number, MemorySpace, Constness>::operator>(
        const Iterator<Number, MemorySpace, OtherConstness> &other) const
      {
        return (other < *this);
      }

    } // namespace SparseMatrixIterators
  }   // namespace TpetraWrappers

} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif // dealii_trilinos_tpetra_sparse_matrix_h
