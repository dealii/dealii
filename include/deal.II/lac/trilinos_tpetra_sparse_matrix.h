// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef dealii_trilinos_tpetra_sparse_matrix_h
#define dealii_trilinos_tpetra_sparse_matrix_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>

// Tpetra includes
#  include <Tpetra_Core.hpp>
#  include <Tpetra_CrsMatrix.hpp>


DEAL_II_NAMESPACE_OPEN

// forward declarations
#  ifndef DOXYGEN
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename MemorySpace>
    class SparsityPattern;
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
    class SparseMatrix : public Subscriptor
    {
    public:
      /**
       * Declare the type for container size.
       */
      using size_type = dealii::types::global_dof_index;

      /**
       * Declare an alias for the type used to store matrix elements, in analogy
       * to all the other container classes.
       */
      using value_type = Number;

      /**
       * Typedef for the NodeType
       */
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 2, 0)
      using NodeType = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<
        typename MemorySpace::kokkos_space::execution_space,
        typename MemorySpace::kokkos_space>;
#  else
      using NodeType = Kokkos::Compat::KokkosDeviceWrapperNode<
        typename MemorySpace::kokkos_space::execution_space,
        typename MemorySpace::kokkos_space>;
#  endif

      /**
       * Typedef for Tpetra::CrsMatrix
       */
      using MatrixType =
        Tpetra::CrsMatrix<Number,
                          int,
                          dealii::types::signed_global_dof_index,
                          NodeType>;

      /**
       * Typedef for Tpetra::Map
       */
      using MapType =
        Tpetra::Map<int, dealii::types::signed_global_dof_index, NodeType>;

      /**
       * Typedef for Tpetra::CrsGraph
       */
      using GraphType =
        Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index, NodeType>;

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
      void
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
      void
      reinit(const IndexSet            &row_parallel_partitioning,
             const IndexSet            &col_parallel_partitioning,
             const SparsityPatternType &sparsity_pattern,
             const MPI_Comm             communicator  = MPI_COMM_WORLD,
             const bool                 exchange_data = false);
      /** @} */

      /**
       * @name Information on the matrix
       */
      /** @{ */
      /**
       * Return the number of rows in this matrix.
       */
      dealii::types::signed_global_dof_index
      m() const;

      /**
       * Return the number of columns in this matrix.
       */
      dealii::types::signed_global_dof_index
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
       * Return the total number of nonzero elements of this matrix (summed
       * over all MPI processes).
       */
      size_t
      n_nonzero_elements() const;

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
       * Add @p value to the element (<i>i,j</i>).
       * Just as the respective call in deal.II SparseMatrix<Number,
       * MemorySpace> class. Moreover, if <tt>value</tt> is not a finite number
       * an exception is thrown.
       *
       * @note When add is called on a compressed matrix, the matrix is set
       * back to an uncompressed state.
       */
      void
      add(const size_type i, const size_type j, const TrilinosScalar value);

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
      add(const size_type       row,
          const size_type       n_cols,
          const size_type      *col_indices,
          const TrilinosScalar *values,
          const bool            elide_zero_values      = true,
          const bool            col_indices_are_sorted = false);
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
      void
      vmult(Vector<Number> &dst, const Vector<Number> &src) const;

      /*
       * Matrix-vector multiplication: let <i>dst = M<sup>T</sup>*src</i> with
       * <i>M</i> being this matrix. This function does the same as vmult() but
       * takes the transposed matrix.
       *
       * Source and destination must not be the same vector.
       */
      void
      Tvmult(Vector<Number> &dst, const Vector<Number> &src) const;

      /**
       * Adding matrix-vector multiplication. Add <i>M*src</i> on <i>dst</i>
       * with <i>M</i> being this matrix.
       *
       * Source and destination must not be the same vector.
       */
      void
      vmult_add(Vector<Number> &dst, const Vector<Number> &src) const;


      /**
       * Adding matrix-vector multiplication. Add <i>M<sup>T</sup>*src</i> to
       * <i>dst</i> with <i>M</i> being this matrix. This function does the same
       * as vmult_add() but takes the transposed matrix.
       *
       * Source and destination must not be the same vector.
       */
      void
      Tvmult_add(Vector<Number> &dst, const Vector<Number> &src) const;
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
      const MatrixType &
      trilinos_matrix() const;

      /**
       * Return a (modifiable) reference to the underlying Trilinos
       * <a
       * href="https://docs.trilinos.org/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html">Tpetra::CrsMatrix</a>
       * class.
       */
      MatrixType &
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
      Teuchos::RCP<const MatrixType>
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
      Teuchos::RCP<MatrixType>
      trilinos_rcp();
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
      DeclExceptionMsg(ExcColMapMissmatch,
                       "The column partitioning of a matrix does not match "
                       "the partitioning of a vector you are trying to "
                       "multiply it with. Are you multiplying the "
                       "matrix with a vector that has ghost elements?");

      /*
       * Exception
       */
      DeclExceptionMsg(ExcDomainMapMissmatch,
                       "The row partitioning of a matrix does not match "
                       "the partitioning of a vector you are trying to "
                       "put the result of a matrix-vector product in. "
                       "Are you trying to put the product of the "
                       "matrix with a vector into a vector that has "
                       "ghost elements?");
      /** @} */

    private:
      /**
       * Pointer to the user-supplied Tpetra Trilinos mapping of the matrix
       * columns that assigns parts of the matrix to the individual processes.
       *
       * @note The Trilinos matrix is row-oriented, and the row_space_map is
       * therefore stored in the Trilinos matrix itself. The additional
       * information from the column space map is used to speed up the
       * assembly process.
       */
      Teuchos::RCP<MapType> column_space_map;

      /**
       * A sparse matrix object in Trilinos to be used for finite element based
       * problems which allows for assembling into non-local elements.  The
       * actual type, a sparse matrix, is set in the constructor.
       */
      Teuchos::RCP<MatrixType> matrix;

      /**
       * A boolean variable to hold information on whether the matrix is
       * fill complete or if the matrix is in compute mode.
       */
      bool compressed;

    }; // class SparseMatrix


    /* ------------------------- Inline functions ---------------------- */

    template <typename Number, typename MemorySpace>
    inline void
    SparseMatrix<Number, MemorySpace>::add(const size_type      i,
                                           const size_type      j,
                                           const TrilinosScalar value)
    {
      add(i, 1, &j, &value, false);
    }



    template <typename Number, typename MemorySpace>
    inline dealii::types::signed_global_dof_index
    SparseMatrix<Number, MemorySpace>::m() const
    {
      return matrix->getRowMap()->getGlobalNumElements();
    }



    template <typename Number, typename MemorySpace>
    inline dealii::types::signed_global_dof_index
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
    inline const Tpetra::CrsMatrix<
      Number,
      int,
      types::signed_global_dof_index,
      typename SparseMatrix<Number, MemorySpace>::NodeType> &
    SparseMatrix<Number, MemorySpace>::trilinos_matrix() const
    {
      return *matrix;
    }



    template <typename Number, typename MemorySpace>
    inline Tpetra::CrsMatrix<
      Number,
      int,
      types::signed_global_dof_index,
      typename SparseMatrix<Number, MemorySpace>::NodeType> &
    SparseMatrix<Number, MemorySpace>::trilinos_matrix()
    {
      return *matrix;
    }



    template <typename Number, typename MemorySpace>
    inline Teuchos::RCP<const Tpetra::CrsMatrix<
      Number,
      int,
      types::signed_global_dof_index,
      typename SparseMatrix<Number, MemorySpace>::NodeType>>
    SparseMatrix<Number, MemorySpace>::trilinos_rcp() const
    {
      return matrix.getConst();
    }



    template <typename Number, typename MemorySpace>
    inline Teuchos::RCP<
      Tpetra::CrsMatrix<Number,
                        int,
                        types::signed_global_dof_index,
                        typename SparseMatrix<Number, MemorySpace>::NodeType>>
    SparseMatrix<Number, MemorySpace>::trilinos_rcp()
    {
      return matrix;
    }

  } // namespace TpetraWrappers

} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif // dealii_trilinos_tpetra_sparse_matrix_h
