// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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

#ifndef __deal2__trilinos_sparse_matrix_h
#define __deal2__trilinos_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/std_cxx11/shared_ptr.h>
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector_view.h>

#  include <vector>
#  include <cmath>
#  include <memory>

#  define TrilinosScalar double
#  include <Epetra_FECrsMatrix.h>
#  include <Epetra_Map.h>
#  include <Epetra_CrsGraph.h>
#  include <Epetra_MultiVector.h>
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#    include "mpi.h"
#  else
#    include "Epetra_SerialComm.h"
#  endif

class Epetra_Export;

DEAL_II_NAMESPACE_OPEN

// forward declarations
template <typename MatrixType> class BlockMatrixBase;

template <typename number> class SparseMatrix;
class SparsityPattern;


namespace TrilinosWrappers
{
  // forward declarations
  class SparseMatrix;
  class SparsityPattern;

  /**
   * Iterators for Trilinos matrices
   */
  namespace SparseMatrixIterators
  {
    // forward declaration
    template <bool Constness> class Iterator;

    /**
     * Exception
     */
    DeclException0 (ExcBeyondEndOfMatrix);

    /**
     * Exception
     */
    DeclException3 (ExcAccessToNonlocalRow,
                    std::size_t, std::size_t, std::size_t,
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
     *
     * @author Guido Kanschat
     * @date 2012
     */
    class AccessorBase
    {
    public:
      /**
       * Declare the type for container size.
       */
      typedef dealii::types::global_dof_index size_type;

      /**
       * Constructor.
       */
      AccessorBase (SparseMatrix *matrix,
                    const size_type  row,
                    const size_type  index);

      /**
       * Row number of the element represented by this object.
       */
      size_type row() const;

      /**
       * Index in row of the element represented by this object.
       */
      size_type index() const;

      /**
       * Column number of the element represented by this object.
       */
      size_type column() const;

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
      void visit_present_row ();

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
      std_cxx11::shared_ptr<std::vector<size_type> > colnum_cache;

      /**
       * Cache for the values of this row.
       */
      std_cxx11::shared_ptr<std::vector<TrilinosScalar> > value_cache;
    };

    /**
     * General template for sparse matrix accessors. The first template
     * argument denotes the underlying numeric type, the second the
     * constness of the matrix.
     *
     * The general template is not implemented, only the specializations
     * for the two possible values of the second template
     * argument. Therefore, the interface listed here only serves as a
     * template provided since doxygen does not link the specializations.
     */
    template <bool Constess>
    class Accessor : public AccessorBase
    {
      /**
      * Value of this matrix entry.
      */
      TrilinosScalar value() const;

      /**
       * Value of this matrix entry.
       */
      TrilinosScalar &value();
    };

    /**
     * The specialization for a const Accessor.
     */
    template<>
    class Accessor<true> : public AccessorBase
    {
    public:
      /**
       * Typedef for the type (including constness) of the matrix to be used
       * here.
       */
      typedef const SparseMatrix MatrixType;

      /**
       * Constructor. Since we use accessors only for read access, a const
       * matrix pointer is sufficient.
       */
      Accessor (MatrixType *matrix,
                const size_type  row,
                const size_type  index);

      /**
       * Copy constructor to get from a const or non-const accessor to a const
       * accessor.
       */
      template <bool Other>
      Accessor (const Accessor<Other> &a);

      /**
       * Value of this matrix entry.
       */
      TrilinosScalar value() const;

    private:
      /**
       * Make iterator class a friend.
       */
      template <bool> friend class Iterator;
    };

    /**
     * The specialization for a mutable Accessor.
     */
    template<>
    class Accessor<false> : public AccessorBase
    {
      class Reference
      {
      public:
        /**
         * Constructor.
         */
        Reference (const Accessor<false> &accessor);

        /**
         * Conversion operator to the data type of the matrix.
         */
        operator TrilinosScalar () const;

        /**
         * Set the element of the matrix we presently point to to @p n.
         */
        const Reference &operator = (const TrilinosScalar n) const;

        /**
         * Add @p n to the element of the matrix we presently point to.
         */
        const Reference &operator += (const TrilinosScalar n) const;

        /**
         * Subtract @p n from the element of the matrix we presently point to.
         */
        const Reference &operator -= (const TrilinosScalar n) const;

        /**
         * Multiply the element of the matrix we presently point to by @p n.
         */
        const Reference &operator *= (const TrilinosScalar n) const;

        /**
         * Divide the element of the matrix we presently point to by @p n.
         */
        const Reference &operator /= (const TrilinosScalar n) const;

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
      typedef SparseMatrix MatrixType;

      /**
       * Constructor. Since we use accessors only for read access, a const
       * matrix pointer is sufficient.
       */
      Accessor (MatrixType *matrix,
                const size_type  row,
                const size_type  index);

      /**
       * Value of this matrix entry.
       */
      Reference value() const;

    private:
      /**
       * Make iterator class a friend.
       */
      template <bool> friend class Iterator;
      /**
       * Make Reference object a friend.
       */
      friend class Reference;
    };

    /**
     * STL conforming iterator. This class acts as an iterator walking
     * over the elements of Trilinos matrices. The implementation of this
     * class is similar to the one for PETSc matrices.
     *
     * Note that Trilinos stores the elements within each row in ascending
     * order. This is opposed to the deal.II sparse matrix style where the
     * diagonal element (if it exists) is stored before all other values, and
     * the PETSc sparse matrices, where one can't guarantee a certain order of
     * the elements.
     *
     * @ingroup TrilinosWrappers
     * @author Martin Kronbichler, Wolfgang Bangerth, 2008
     */
    template <bool Constness>
    class Iterator
    {
    public:
      /**
       * Declare type for container size.
       */
      typedef dealii::types::global_dof_index size_type;

      /**
       * Typedef for the matrix type (including constness) we are to operate
       * on.
       */
      typedef typename Accessor<Constness>::MatrixType MatrixType;

      /**
       * Constructor. Create an iterator into the matrix @p matrix for the
       * given row and the index within it.
       */
      Iterator (MatrixType *matrix,
                const size_type  row,
                const size_type  index);

      /**
       * Copy constructor with optional change of constness.
       */
      template <bool Other>
      Iterator(const Iterator<Other> &other);

      /**
       * Prefix increment.
       */
      Iterator<Constness> &operator++ ();

      /**
       * Postfix increment.
       */
      Iterator<Constness> operator++ (int);

      /**
       * Dereferencing operator.
       */
      const Accessor<Constness> &operator* () const;

      /**
       * Dereferencing operator.
       */
      const Accessor<Constness> *operator-> () const;

      /**
       * Comparison. True, if both iterators point to the same matrix
       * position.
       */
      bool operator == (const Iterator<Constness> &) const;

      /**
       * Inverse of <tt>==</tt>.
       */
      bool operator != (const Iterator<Constness> &) const;

      /**
       * Comparison operator. Result is true if either the first row number is
       * smaller or if the row numbers are equal and the first index is
       * smaller.
       */
      bool operator < (const Iterator<Constness> &) const;

      /**
       * Comparison operator. The opposite of the previous operator
       */
      bool operator > (const Iterator<Constness> &) const;

      /**
       * Exception
       */
      DeclException2 (ExcInvalidIndexWithinRow,
                      size_type, size_type,
                      << "Attempt to access element " << arg2
                      << " of row " << arg1
                      << " which doesn't have that many elements.");

    private:
      /**
       * Store an object of the accessor class.
       */
      Accessor<Constness> accessor;

      template <bool Other> friend class Iterator;
    };

  }


  /**
   * This class implements a wrapper to use the Trilinos distributed
   * sparse matrix class Epetra_FECrsMatrix. This is precisely the kind of
   * matrix we deal with all the time - we most likely get it from some
   * assembly process, where also entries not locally owned might need to
   * be written and hence need to be forwarded to the owner process.  This
   * class is designed to be used in a distributed memory architecture
   * with an MPI compiler on the bottom, but works equally well also for
   * serial processes. The only requirement for this class to work is that
   * Trilinos has been installed with the same compiler as is used for
   * generating deal.II.
   *
   * The interface of this class is modeled after the existing
   * SparseMatrix class in deal.II. It has almost the same member
   * functions, and is often exchangable. However, since Trilinos only
   * supports a single scalar type (double), it is not templated, and only
   * works with doubles.
   *
   * Note that Trilinos only guarantees that operations do what you expect
   * if the functions @p GlobalAssemble has been called after matrix
   * assembly.  Therefore, you need to call SparseMatrix::compress()
   * before you actually use the matrix. This also calls @p FillComplete
   * that compresses the storage format for sparse matrices by discarding
   * unused elements. Trilinos allows to continue with assembling the
   * matrix after calls to these functions, though.
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
   *   <li> The matrix uses only one MPI process.
   *   <li> The matrix has been initialized with the reinit() method
   *   with a CompressedSimpleSparsityPattern (that includes the set of locally
   *   relevant rows, i.e., the rows that an assembly routine will possibly
   *   write into).
   *   <li> The matrix has been initialized from a
   *   TrilinosWrappers::SparsityPattern object that in turn has been
   *   initialized with the reinit function specifying three index sets, one
   *   for the rows, one for the columns and for the larger set of @p
   *   writeable_rows, and the operation is an addition. If Trilinos version
   *   11.10 and greater is used, initializing from a
   *   TrilinosWrappers::SparsityPattern that has been filled by a function
   *   similar to DoFTools::make_sparsity_pattern always results in a matrix
   *   that allows several processes to write into the same matrix row.
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
   * @author Martin Kronbichler, Wolfgang Bangerth, 2008, 2009
   */
  class SparseMatrix : public Subscriptor
  {
  public:
    /**
     * Declare the type for container size.
     */
    typedef dealii::types::global_dof_index size_type;

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
     * Declare a typedef for the iterator class.
     */
    typedef SparseMatrixIterators::Iterator<false> iterator;

    /**
     * Declare a typedef for the const iterator class.
     */
    typedef SparseMatrixIterators::Iterator<true> const_iterator;

    /**
     * Declare a typedef in analogy to all the other container classes.
     */
    typedef TrilinosScalar value_type;

    /**
     * @name Constructors and initialization.
     */
//@{
    /**
     * Default constructor. Generates an empty (zero-size) matrix.
     */
    SparseMatrix ();

    /**
     * Generate a matrix that is completely stored locally, having #m rows and
     * #n columns.
     *
     * The number of columns entries per row is specified as the maximum
     * number of entries argument.
     */
    SparseMatrix (const size_type  m,
                  const size_type  n,
                  const unsigned int  n_max_entries_per_row);

    /**
     * Generate a matrix that is completely stored locally, having #m rows and
     * #n columns.
     *
     * The vector <tt>n_entries_per_row</tt> specifies the number of entries
     * in each row.
     */
    SparseMatrix (const size_type                  m,
                  const size_type                  n,
                  const std::vector<unsigned int> &n_entries_per_row);

    /**
     * Generate a matrix from a Trilinos sparsity pattern object.
     */
    SparseMatrix (const SparsityPattern &InputSparsityPattern);

    /**
     * Copy constructor. Sets the calling matrix to be the same as the input
     * matrix, i.e., using the same sparsity pattern and entries.
     */
    SparseMatrix (const SparseMatrix &InputMatrix);

    /**
     * Destructor. Made virtual so that one can use pointers to this class.
     */
    virtual ~SparseMatrix ();

    /**
     * This function initializes the Trilinos matrix with a deal.II sparsity
     * pattern, i.e. it makes the Trilinos Epetra matrix know the position of
     * nonzero entries according to the sparsity pattern. This function is
     * meant for use in serial programs, where there is no need to specify how
     * the matrix is going to be distributed among different processors. This
     * function works in %parallel, too, but it is recommended to manually
     * specify the %parallel partioning of the matrix using an
     * Epetra_Map. When run in %parallel, it is currently necessary that each
     * processor holds the sparsity_pattern structure because each processor
     * sets its rows.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template<typename SparsityType>
    void reinit (const SparsityType &sparsity_pattern);

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
    void reinit (const SparsityPattern &sparsity_pattern);

    /**
     * This function copies the content in <tt>sparse_matrix</tt> to the
     * calling matrix.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    void reinit (const SparseMatrix &sparse_matrix);

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
     * (i.e., one that differs from the one used in the sparse matrix given
     * in the first argument), then the resulting Trilinos matrix will have
     * the sparsity pattern so given. This of course also means that all
     * entries in the given matrix that are not part of this separate
     * sparsity pattern will in fact be dropped.
     */
    template <typename number>
    void reinit (const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                 const double                          drop_tolerance=1e-13,
                 const bool                            copy_values=true,
                 const ::dealii::SparsityPattern      *use_this_sparsity=0);

    /**
     * This reinit function takes as input a Trilinos Epetra_CrsMatrix and
     * copies its sparsity pattern. If so requested, even the content (values)
     * will be copied.
     */
    void reinit (const Epetra_CrsMatrix &input_matrix,
                 const bool              copy_values = true);
//@}
    /**
     * @name Constructors and initialization using an Epetra_Map description
     */
//@{
    /**
     * Constructor using an Epetra_Map to describe the %parallel
     * partitioning. The parameter @p n_max_entries_per_row sets the number of
     * nonzero entries in each row that will be allocated. Note that this
     * number does not need to be exact, and it is even allowed that the
     * actual matrix structure has more nonzero entries than specified in the
     * constructor. However it is still advantageous to provide good estimates
     * here since this will considerably increase the performance of the
     * matrix setup. However, there is no effect in the performance of
     * matrix-vector products, since Trilinos reorganizes the matrix memory
     * prior to use (in the compress() step).
     */
    SparseMatrix (const Epetra_Map  &parallel_partitioning,
                  const size_type    n_max_entries_per_row = 0);

    /**
     * Same as before, but now set a value of nonzeros for each matrix
     * row. Since we know the number of elements in the matrix exactly in this
     * case, we can already allocate the right amount of memory, which makes
     * the creation process including the insertion of nonzero elements by the
     * respective SparseMatrix::reinit call considerably faster.
     */
    SparseMatrix (const Epetra_Map                &parallel_partitioning,
                  const std::vector<unsigned int> &n_entries_per_row);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different Epetra maps for rows and columns. This interface is meant to
     * be used for generating rectangular matrices, where one map describes
     * the %parallel partitioning of the dofs associated with the matrix rows
     * and the other one the partitioning of dofs in the matrix columns. Note
     * that there is no real parallelism along the columns &ndash; the
     * processor that owns a certain row always owns all the column elements,
     * no matter how far they might be spread out. The second Epetra_Map is
     * only used to specify the number of columns and for internal
     * arrangements when doing matrix-vector products with vectors based on
     * that column map.
     *
     * The integer input @p n_max_entries_per_row defines the number of
     * columns entries per row that will be allocated.
     */
    SparseMatrix (const Epetra_Map &row_parallel_partitioning,
                  const Epetra_Map &col_parallel_partitioning,
                  const size_type   n_max_entries_per_row = 0);

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
    SparseMatrix (const Epetra_Map                &row_parallel_partitioning,
                  const Epetra_Map                &col_parallel_partitioning,
                  const std::vector<unsigned int> &n_entries_per_row);

    /**
     * This function is initializes the Trilinos Epetra matrix according to
     * the specified sparsity_pattern, and also reassigns the matrix rows to
     * different processes according to a user-supplied Epetra map. In
     * programs following the style of the tutorial programs, this function
     * (and the respective call for a rectangular matrix) are the natural way
     * to initialize the matrix size, its distribution among the MPI processes
     * (if run in %parallel) as well as the locatoin of non-zero
     * elements. Trilinos stores the sparsity pattern internally, so it won't
     * be needed any more after this call, in contrast to the deal.II own
     * object. The optional argument @p exchange_data can be used for
     * reinitialization with a sparsity pattern that is not fully
     * constructed. This feature is only implemented for input sparsity
     * patterns of type CompressedSimpleSparsityPattern. If the flag is not
     * set, each processor just sets the elements in the sparsity pattern that
     * belong to its rows.
     *
     * If the sparsity pattern given to this function is of type
     * CompressedSimpleSparsity pattern, then a matrix will be created that
     * allows several threads to write into different rows of the matrix at
     * the same also with MPI, as opposed to most other reinit() methods.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template<typename SparsityType>
    void reinit (const Epetra_Map    &parallel_partitioning,
                 const SparsityType  &sparsity_pattern,
                 const bool          exchange_data = false);

    /**
     * This function is similar to the other initialization function above,
     * but now also reassigns the matrix rows and columns according to two
     * user-supplied Epetra maps.  To be used for rectangular matrices. The
     * optional argument @p exchange_data can be used for reinitialization
     * with a sparsity pattern that is not fully constructed. This feature is
     * only implemented for input sparsity patterns of type
     * CompressedSimpleSparsityPattern.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template<typename SparsityType>
    void reinit (const Epetra_Map    &row_parallel_partitioning,
                 const Epetra_Map    &col_parallel_partitioning,
                 const SparsityType  &sparsity_pattern,
                 const bool          exchange_data = false);

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
    void reinit (const Epetra_Map                     &parallel_partitioning,
                 const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                 const double                          drop_tolerance=1e-13,
                 const bool                            copy_values=true,
                 const ::dealii::SparsityPattern      *use_this_sparsity=0);

    /**
     * This function is similar to the other initialization function with
     * deal.II sparse matrix input above, but now takes Epetra maps for both
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
    void reinit (const Epetra_Map                      &row_parallel_partitioning,
                 const Epetra_Map                      &col_parallel_partitioning,
                 const ::dealii::SparseMatrix<number>  &dealii_sparse_matrix,
                 const double                           drop_tolerance=1e-13,
                 const bool                             copy_values=true,
                 const ::dealii::SparsityPattern      *use_this_sparsity=0);
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
     * matrix setup. However, there is no effect in the performance of
     * matrix-vector products, since Trilinos reorganizes the matrix memory
     * prior to use (in the compress() step).
     */
    SparseMatrix (const IndexSet    &parallel_partitioning,
                  const MPI_Comm    &communicator = MPI_COMM_WORLD,
                  const unsigned int n_max_entries_per_row = 0);

    /**
     * Same as before, but now set the number of nonzeros in each matrix row
     * separately. Since we know the number of elements in the matrix exactly
     * in this case, we can already allocate the right amount of memory, which
     * makes the creation process including the insertion of nonzero elements
     * by the respective SparseMatrix::reinit call considerably faster.
     */
    SparseMatrix (const IndexSet                  &parallel_partitioning,
                  const MPI_Comm                  &communicator,
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
    SparseMatrix (const IndexSet  &row_parallel_partitioning,
                  const IndexSet  &col_parallel_partitioning,
                  const MPI_Comm  &communicator = MPI_COMM_WORLD,
                  const size_type  n_max_entries_per_row = 0);

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
    SparseMatrix (const IndexSet                  &row_parallel_partitioning,
                  const IndexSet                  &col_parallel_partitioning,
                  const MPI_Comm                  &communicator,
                  const std::vector<unsigned int> &n_entries_per_row);

    /**
     * This function is initializes the Trilinos Epetra matrix according to
     * the specified sparsity_pattern, and also reassigns the matrix rows to
     * different processes according to a user-supplied index set and
     * %parallel communicator. In programs following the style of the tutorial
     * programs, this function (and the respective call for a rectangular
     * matrix) are the natural way to initialize the matrix size, its
     * distribution among the MPI processes (if run in %parallel) as well as
     * the locatoin of non-zero elements. Trilinos stores the sparsity pattern
     * internally, so it won't be needed any more after this call, in contrast
     * to the deal.II own object. The optional argument @p exchange_data can
     * be used for reinitialization with a sparsity pattern that is not fully
     * constructed. This feature is only implemented for input sparsity
     * patterns of type CompressedSimpleSparsityPattern. If the flag is not
     * set, each processor just sets the elements in the sparsity pattern that
     * belong to its rows.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template<typename SparsityType>
    void reinit (const IndexSet      &parallel_partitioning,
                 const SparsityType  &sparsity_pattern,
                 const MPI_Comm      &communicator = MPI_COMM_WORLD,
                 const bool           exchange_data = false);

    /**
     * This function is similar to the other initialization function above,
     * but now also reassigns the matrix rows and columns according to two
     * user-supplied index sets.  To be used for rectangular matrices. The
     * optional argument @p exchange_data can be used for reinitialization
     * with a sparsity pattern that is not fully constructed. This feature is
     * only implemented for input sparsity patterns of type
     * CompressedSimpleSparsityPattern.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    template<typename SparsityType>
    void reinit (const IndexSet      &row_parallel_partitioning,
                 const IndexSet      &col_parallel_partitioning,
                 const SparsityType  &sparsity_pattern,
                 const MPI_Comm      &communicator = MPI_COMM_WORLD,
                 const bool           exchange_data = false);

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
    void reinit (const IndexSet                       &parallel_partitioning,
                 const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                 const MPI_Comm                       &communicator = MPI_COMM_WORLD,
                 const double                          drop_tolerance=1e-13,
                 const bool                            copy_values=true,
                 const ::dealii::SparsityPattern      *use_this_sparsity=0);

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
    void reinit (const IndexSet                        &row_parallel_partitioning,
                 const IndexSet                        &col_parallel_partitioning,
                 const ::dealii::SparseMatrix<number>  &dealii_sparse_matrix,
                 const MPI_Comm                        &communicator = MPI_COMM_WORLD,
                 const double                           drop_tolerance=1e-13,
                 const bool                             copy_values=true,
                 const ::dealii::SparsityPattern      *use_this_sparsity=0);
//@}
    /**
     * @name Information on the matrix
     */
//@{

    /**
     * Return the number of rows in this matrix.
     */
    size_type m () const;

    /**
     * Return the number of columns in this matrix.
     */
    size_type n () const;

    /**
     * Return the local dimension of the matrix, i.e. the number of rows
     * stored on the present MPI process. For sequential matrices, this number
     * is the same as m(), but for %parallel matrices it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range().
     */
    unsigned int local_size () const;

    /**
     * Return a pair of indices indicating which rows of this matrix are
     * stored locally. The first number is the index of the first row stored,
     * the second the index of the one past the last one that is stored
     * locally. If this is a sequential matrix, then the result will be the
     * pair (0,m()), otherwise it will be a pair (i,i+n), where
     * <tt>n=local_size()</tt>.
     */
    std::pair<size_type, size_type>
    local_range () const;

    /**
     * Return whether @p index is in the local range or not, see also
     * local_range().
     */
    bool in_local_range (const size_type index) const;

    /**
     * Return the number of nonzero elements of this matrix.
     */
    size_type n_nonzero_elements () const;

    /**
     * Number of entries in a specific row.
     */
    unsigned int row_length (const size_type row) const;

    /**
     * Returns the state of the matrix, i.e., whether compress() needs to be
     * called after an operation requiring data exchange. A call to compress()
     * is also needed when the method set() has been called (even when working
     * in serial).
     */
    bool is_compressed () const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object. Note that only the memory reserved on the current processor is
     * returned in case this is called in an MPI-based program.
     */
    size_type memory_consumption () const;

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
    operator = (const double d);

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    void clear ();

    /**
     * This command does two things:
     * <ul>
     * <li> If the matrix was initialized without a sparsity pattern, elements
     * have been added manually using the set() command. When this process is
     * completed, a call to compress() reorganizes the internal data
     * structures (aparsity pattern) so that a fast access to data is possible
     * in matrix-vector products.
     * <li> If the matrix structure has already been fixed (either by
     * initialization with a sparsity pattern or by calling compress() during
     * the setup phase), this command does the %parallel exchange of
     * data. This is necessary when we perform assembly on more than one (MPI)
     * process, because then some non-local row data will accumulate on nodes
     * that belong to the current's processor element, but are actually held
     * by another. This command is usually called after all elements have been
     * traversed.
     * </ul>
     *
     * In both cases, this function compresses the data structures and allows
     * the resulting matrix to be used in all other operations like
     * matrix-vector products. This is a collective operation, i.e., it needs
     * to be run on all processors when used in %parallel.
     *
     * See @ref GlossCompress "Compressing distributed objects"
     * for more information.
     */
    void compress (::dealii::VectorOperation::values operation);

    /**
     * @deprecated: use compress() with VectorOperation instead.
     */
    void compress () DEAL_II_DEPRECATED;

    /**
     * Set the element (<i>i,j</i>) to @p value.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. When compress() is called for the first time, then this is no
     * longer possible and an insertion of elements at positions which have
     * not been initialized will throw an exception. Note that in case
     * elements need to be inserted, it is mandatory that elements are
     * inserted only once. Otherwise, the elements will actually be added in
     * the end (since it is not possible to efficiently find values to the
     * same entry before compress() has been called). In the case that an
     * element is set more than once, initialize the matrix with a sparsity
     * pattern first.
     */
    void set (const size_type i,
              const size_type j,
              const TrilinosScalar value);

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
     * extended. When compress() is called for the first time, then this is no
     * longer possible and an insertion of elements at positions which have
     * not been initialized will throw an exception.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     */
    void set (const std::vector<size_type>     &indices,
              const FullMatrix<TrilinosScalar> &full_matrix,
              const bool                        elide_zero_values = false);

    /**
     * Same function as before, but now including the possibility to use
     * rectangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     */
    void set (const std::vector<size_type>     &row_indices,
              const std::vector<size_type>     &col_indices,
              const FullMatrix<TrilinosScalar> &full_matrix,
              const bool                        elide_zero_values = false);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by <tt>col_indices</tt> to the respective value.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. When compress() is called for the first time, then this is no
     * longer possible and an insertion of elements at positions which have
     * not been initialized will throw an exception.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     */
    void set (const size_type                    row,
              const std::vector<size_type>      &col_indices,
              const std::vector<TrilinosScalar> &values,
              const bool                         elide_zero_values = false);

    /**
     * Set several elements to values given by <tt>values</tt> in a given row
     * in columns given by col_indices into the sparse matrix.
     *
     * This function is able to insert new elements into the matrix as long as
     * compress() has not been called, so the sparsity pattern will be
     * extended. When compress() is called for the first time, then this is no
     * longer possible and an insertion of elements at positions which have
     * not been initialized will throw an exception.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     */
    void set (const size_type       row,
              const size_type       n_cols,
              const size_type      *col_indices,
              const TrilinosScalar *values,
              const bool            elide_zero_values = false);

    /**
     * Add @p value to the element (<i>i,j</i>).
     *
     * Just as the respective call in deal.II SparseMatrix<Number> class (but
     * in contrast to the situation for PETSc based matrices), this function
     * throws an exception if an entry does not exist in the sparsity
     * pattern. Moreover, if <tt>value</tt> is not a finite number an
     * exception is thrown.
     */
    void add (const size_type      i,
              const size_type      j,
              const TrilinosScalar value);

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
    void add (const std::vector<size_type>  &indices,
              const FullMatrix<TrilinosScalar> &full_matrix,
              const bool                        elide_zero_values = true);

    /**
     * Same function as before, but now including the possibility to use
     * rectangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     */
    void add (const std::vector<size_type>     &row_indices,
              const std::vector<size_type>     &col_indices,
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
    void add (const size_type                    row,
              const std::vector<size_type>      &col_indices,
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
    void add (const size_type       row,
              const size_type       n_cols,
              const size_type      *col_indices,
              const TrilinosScalar *values,
              const bool            elide_zero_values = true,
              const bool            col_indices_are_sorted = false);

    /**
     * Multiply the entire matrix by a fixed factor.
     */
    SparseMatrix &operator *= (const TrilinosScalar factor);

    /**
     * Divide the entire matrix by a fixed factor.
     */
    SparseMatrix &operator /= (const TrilinosScalar factor);

    /**
     * Copy the given (Trilinos) matrix
     * (sparsity pattern and entries).
     */
    void copy_from (const SparseMatrix &source);

    /**
     * Add <tt>matrix</tt> scaled by <tt>factor</tt> to this matrix, i.e. the
     * matrix <tt>factor*matrix</tt> is added to <tt>this</tt>. If the
     * sparsity pattern of the calling matrix does not contain all the
     * elements in the sparsity pattern of the input matrix, this function
     * will throw an exception.
     */
    void add (const TrilinosScalar  factor,
              const SparseMatrix   &matrix);

    /**
     * Remove all elements from this <tt>row</tt> by setting them to zero. The
     * function does not modify the number of allocated nonzero entries, it
     * only sets some entries to zero. It may drop them from the sparsity
     * pattern, though (but retains the allocated memory in case new entries
     * are again added later). Note that this is a global operation, so this
     * needs to be done on all MPI processes.
     *
     * This operation is used in eliminating constraints (e.g. due to hanging
     * nodes) and makes sure that we can write this modification to the matrix
     * without having to read entries (such as the locations of non-zero
     * elements) from it &mdash; without this operation, removing constraints
     * on %parallel matrices is a rather complicated procedure.
     *
     * The second parameter can be used to set the diagonal entry of this row
     * to a value different from zero. The default is to set it to zero.
     */
    void clear_row (const size_type      row,
                    const TrilinosScalar new_diag_value = 0);

    /**
     * Same as clear_row(), except that it works on a number of rows at once.
     *
     * The second parameter can be used to set the diagonal entries of all
     * cleared rows to something different from zero. Note that all of these
     * diagonal entries get the same value -- if you want different values for
     * the diagonal entries, you have to set them by hand.
     */
    void clear_rows (const std::vector<size_type> &rows,
                     const TrilinosScalar          new_diag_value = 0);

    /**
     * Sets an internal flag so that all operations performed by the matrix,
     * i.e., multiplications, are done in transposed order. However, this does
     * not reshape the matrix to transposed form directly, so care should be
     * taken when using this flag.
     */
    void transpose ();

//@}
    /**
     * @name Entry Access
     */
//@{

    /**
     * Return the value of the entry (<i>i,j</i>).  This may be an expensive
     * operation and you should always take care where to call this
     * function. As in the deal.II sparse matrix class, we throw an exception
     * if the respective entry doesn't exist in the sparsity pattern of this
     * class, which is requested from Trilinos. Moreover, an exception will be
     * thrown when the requested element is not saved on the calling process.
     */
    TrilinosScalar operator () (const size_type i,
                                const size_type j) const;

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
    TrilinosScalar el (const size_type i,
                       const size_type j) const;

    /**
     * Return the main diagonal element in the <i>i</i>th row. This function
     * throws an error if the matrix is not quadratic and it also throws an
     * error if <i>(i,i)</i> is not element of the local matrix.  See also the
     * comment in trilinos_sparse_matrix.cc.
     */
    TrilinosScalar diag_element (const size_type i) const;

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
     * This function can be called with several different vector objects,
     * namely TrilinosWrappers::Vector, TrilinosWrappers::MPI::Vector as well
     * as deal.II's own vector classes Vector<double> and
     * parallel::distributed::Vector<double>.
     *
     * Note that both vectors have to be distributed vectors generated using
     * the same Map as was used for the matrix in case you work on a
     * distributed memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class (or one of the two derived classes
     * Vector and MPI::Vector).
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     */
    template<typename VectorType>
    void vmult (VectorType       &dst,
                const VectorType &src) const;

    /**
     * Matrix-vector multiplication: let <i>dst = M<sup>T</sup>*src</i> with
     * <i>M</i> being this matrix. This function does the same as vmult() but
     * takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several different vector objects,
     * namely TrilinosWrappers::Vector, TrilinosWrappers::MPI::Vector as well
     * as deal.II's own vector classes Vector<double> and
     * parallel::distributed::Vector<double>.
     *
     * Note that both vectors have to be distributed vectors generated using
     * the same Map as was used for the matrix in case you work on a
     * distributed memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class (or one of the two derived classes
     * Vector and MPI::Vector).
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     */
    template <typename VectorType>
    void Tvmult (VectorType       &dst,
                 const VectorType &src) const;

    /**
     * Adding matrix-vector multiplication. Add <i>M*src</i> on <i>dst</i>
     * with <i>M</i> being this matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several different vector objects,
     * namely TrilinosWrappers::Vector, TrilinosWrappers::MPI::Vector as well
     * as deal.II's own vector classes Vector<double> and
     * parallel::distributed::Vector<double>.
     *
     * When using a vector of type TrilinosWrappers::MPI::Vector, both vectors
     * have to be distributed vectors generated using the same Map as was used
     * for the matrix rows and columns in case you work on a distributed
     * memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class.
     *
     * In case of a localized Vector (i.e., TrilinosWrappers::Vector or
     * Vector<double>), this function will only work when running on one
     * processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     *
     */
    template<typename VectorType>
    void vmult_add (VectorType       &dst,
                    const VectorType &src) const;

    /**
     * Adding matrix-vector multiplication. Add <i>M<sup>T</sup>*src</i> to
     * <i>dst</i> with <i>M</i> being this matrix. This function does the same
     * as vmult_add() but takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     *
     * This function can be called with several different vector objects,
     * namely TrilinosWrappers::Vector, TrilinosWrappers::MPI::Vector as well
     * as deal.II's own vector classes Vector<double> and
     * parallel::distributed::Vector<double>.
     *
     * When using a vector of type TrilinosWrappers::MPI::Vector, both vectors
     * have to be distributed vectors generated using the same Map as was used
     * for the matrix rows and columns in case you work on a distributed
     * memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class.
     *
     * In case of a localized Vector (i.e., TrilinosWrappers::Vector or
     * Vector<double>), this function will only work when running on one
     * processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     */
    template <typename VectorType>
    void Tvmult_add (VectorType       &dst,
                     const VectorType &src) const;

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
     * Note that both vectors have to be distributed vectors generated using
     * the same Map as was used for the matrix in case you work on a
     * distributed memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class (or one of the two derived classes
     * Vector and MPI::Vector).
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     */
    TrilinosScalar matrix_norm_square (const VectorBase &v) const;

    /**
     * Compute the matrix scalar product $\left(u,Mv\right)$.
     *
     * The implementation of this function is not as efficient as the one in
     * the @p SparseMatrix class used in deal.II (i.e. the original one, not
     * the Trilinos wrapper class) since Trilinos doesn't support this
     * operation and needs a temporary vector.
     *
     * Note that both vectors have to be distributed vectors generated using
     * the same Map as was used for the matrix in case you work on a
     * distributed memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class (or one of the two derived classes
     * Vector and MPI::Vector).
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     */
    TrilinosScalar matrix_scalar_product (const VectorBase &u,
                                          const VectorBase &v) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * Source <i>x</i> and destination <i>dst</i> must not be the same vector.
     *
     * Note that both vectors have to be distributed vectors generated using
     * the same Map as was used for the matrix in case you work on a
     * distributed memory architecture, using the interface in the
     * TrilinosWrappers::VectorBase class (or one of the two derived classes
     * Vector and MPI::Vector).
     *
     * In case of a localized Vector, this function will only work when
     * running on one processor, since the matrix object is inherently
     * distributed. Otherwise, and exception will be thrown.
     */
    TrilinosScalar residual (VectorBase       &dst,
                             const VectorBase &x,
                             const VectorBase &b) const;

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
    void mmult (SparseMatrix       &C,
                const SparseMatrix &B,
                const VectorBase   &V = VectorBase()) const;


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
    void Tmmult (SparseMatrix       &C,
                 const SparseMatrix &B,
                 const VectorBase   &V = VectorBase()) const;

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
    TrilinosScalar l1_norm () const;

    /**
     * Return the linfty-norm of the matrix, that is
     * $|M|_\infty=\max_{\mathrm{all\ rows\ } i}\sum_{\mathrm{all\ columns\ }
     * j} |M_{ij}|$, (max. sum of rows).  This is the natural matrix norm that
     * is compatible to the linfty-norm of vectors, i.e.  $|Mv|_\infty \leq
     * |M|_\infty |v|_\infty$.  (cf. Haemmerlin-Hoffmann: Numerische
     * Mathematik)
     */
    TrilinosScalar linfty_norm () const;

    /**
     * Return the frobenius norm of the matrix, i.e. the square root of the
     * sum of squares of all entries in the matrix.
     */
    TrilinosScalar frobenius_norm () const;

//@}
    /**
     * @name Access to underlying Trilinos data
     */
//@{

    /**
     * Return a const reference to the underlying Trilinos Epetra_CrsMatrix
     * data.
     */
    const Epetra_CrsMatrix &trilinos_matrix () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_CrsGraph
     * data that stores the sparsity pattern of the matrix.
     */
    const Epetra_CrsGraph &trilinos_sparsity_pattern () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the domain space of this matrix, i.e., the
     * partitioning of the vectors this matrix has to be multiplied with.
     */
    const Epetra_Map &domain_partitioner () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the range space of this matrix, i.e., the
     * partitioning of the vectors that are result from matrix-vector
     * products.
     */
    const Epetra_Map &range_partitioner () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the matrix rows. Equal to the partitioning of
     * the range.
     */
    const Epetra_Map &row_partitioner () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the matrix columns. This is in general not
     * equal to the partitioner Epetra_Map for the domain because of overlap
     * in the matrix.
     */
    const Epetra_Map &col_partitioner () const;
//@}
    /**
     * @name Iterators
     */
//@{

    /**
     * STL-like iterator with the first entry.
     */
    const_iterator begin () const;

    /**
     * Final iterator.
     */
    const_iterator end () const;

    /**
     * STL-like iterator with the first entry of row @p r.
     *
     * Note that if the given row is empty, i.e. does not contain any nonzero
     * entries, then the iterator returned by this function equals
     * <tt>end(r)</tt>. Note also that the iterator may not be dereferencable
     * in that case.
     */
    const_iterator begin (const size_type r) const;

    /**
     * Final iterator of row <tt>r</tt>. It points to the first element past
     * the end of line @p r, or past the end of the entire sparsity pattern.
     *
     * Note that the end iterator is not necessarily dereferencable. This is
     * in particular the case if it is the end iterator for the last row of a
     * matrix.
     */
    const_iterator end (const size_type r) const;

    /**
     * STL-like iterator with the first entry.
     */
    iterator begin ();

    /**
     * Final iterator.
     */
    iterator end ();

    /**
     * STL-like iterator with the first entry of row @p r.
     *
     * Note that if the given row is empty, i.e. does not contain any nonzero
     * entries, then the iterator returned by this function equals
     * <tt>end(r)</tt>. Note also that the iterator may not be dereferencable
     * in that case.
     */
    iterator begin (const size_type r);

    /**
     * Final iterator of row <tt>r</tt>. It points to the first element past
     * the end of line @p r, or past the end of the entire sparsity pattern.
     *
     * Note that the end iterator is not necessarily dereferencable. This is
     * in particular the case if it is the end iterator for the last row of a
     * matrix.
     */
    iterator end (const size_type r);

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
    void write_ascii ();

    /**
     * Print the matrix to the given stream, using the format <tt>(line,col)
     * value</tt>, i.e. one nonzero entry of the matrix per line. The optional
     * flag outputs the sparsity pattern in Trilinos style, where the data is
     * sorted according to the processor number when printed to the stream, as
     * well as a summary of the matrix like the global size.
     */
    void print (std::ostream &out,
                const bool    write_extended_trilinos_info = false) const;

//@}
    /** @addtogroup Exceptions
     *
     */
//@{
    /**
     * Exception
     */
    DeclException1 (ExcTrilinosError,
                    int,
                    << "An error with error number " << arg1
                    << " occurred while calling a Trilinos function");

    /**
     * Exception
     */
    DeclException2 (ExcInvalidIndex,
                    size_type, size_type,
                    << "The entry with index <" << arg1 << ',' << arg2
                    << "> does not exist.");

    /**
     * Exception
     */
    DeclException0 (ExcSourceEqualsDestination);

    /**
     * Exception
     */
    DeclException0 (ExcMatrixNotCompressed);

    /**
     * Exception
     */
    DeclException4 (ExcAccessToNonLocalElement,
                    size_type, size_type, size_type, size_type,
                    << "You tried to access element (" << arg1
                    << "/" << arg2 << ")"
                    << " of a distributed matrix, but only rows "
                    << arg3 << " through " << arg4
                    << " are stored locally and can be accessed.");

    /**
     * Exception
     */
    DeclException2 (ExcAccessToNonPresentElement,
                    size_type, size_type,
                    << "You tried to access element (" << arg1
                    << "/" << arg2 << ")"
                    << " of a sparse matrix, but it appears to not"
                    << " exist in the Trilinos sparsity pattern.");
//@}



  protected:

    /**
    * For some matrix storage formats, in particular for the PETSc distributed
    * blockmatrices, set and add operations on individual elements can not be
    * freely mixed. Rather, one has to synchronize operations when one wants
    * to switch from setting elements to adding to elements.  BlockMatrixBase
    * automatically synchronizes the access by calling this helper function
    * for each block.  This function ensures that the matrix is in a state
    * that allows adding elements; if it previously already was in this state,
    * the function does nothing.
    */
    void prepare_add();

    /**
    * Same as prepare_add() but prepare the matrix for setting elements if the
    * representation of elements in this class requires such an operation.
    */
    void prepare_set();



  private:

    /**
     * Pointer to the user-supplied Epetra Trilinos mapping of the matrix
     * columns that assigns parts of the matrix to the individual processes.
     */
    std_cxx11::shared_ptr<Epetra_Map> column_space_map;

    /**
     * A sparse matrix object in Trilinos to be used for finite element based
     * problems which allows for assembling into non-local elements.  The
     * actual type, a sparse matrix, is set in the constructor.
     */
    std_cxx11::shared_ptr<Epetra_FECrsMatrix> matrix;

    /**
     * A sparse matrix object in Trilinos to be used for collecting the
     * non-local elements if the matrix was constructed from a Trilinos
     * sparsity pattern with the respective option.
     */
    std_cxx11::shared_ptr<Epetra_CrsMatrix> nonlocal_matrix;

    /**
     * An export object used to communicate the nonlocal matrix.
     */
    std_cxx11::shared_ptr<Epetra_Export>    nonlocal_matrix_exporter;

    /**
     * Trilinos doesn't allow to mix additions to matrix entries and
     * overwriting them (to make synchronisation of %parallel computations
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

    /**
     *  To allow calling protected prepare_add() and prepare_set().
     */
    friend class BlockMatrixBase<SparseMatrix>;
  };



// -------------------------- inline and template functions ----------------------


#ifndef DOXYGEN

  namespace SparseMatrixIterators
  {
    inline
    AccessorBase::AccessorBase(SparseMatrix *matrix, size_type row, size_type index)
      :
      matrix(matrix),
      a_row(row),
      a_index(index)
    {
      visit_present_row ();
    }


    inline
    AccessorBase::size_type
    AccessorBase::row() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_row;
    }


    inline
    AccessorBase::size_type
    AccessorBase::column() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*colnum_cache)[a_index];
    }


    inline
    AccessorBase::size_type
    AccessorBase::index() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_index;
    }


    inline
    Accessor<true>::Accessor (MatrixType *matrix,
                              const size_type  row,
                              const size_type  index)
      :
      AccessorBase(const_cast<SparseMatrix *>(matrix), row, index)
    {}


    template <bool Other>
    inline
    Accessor<true>::Accessor(const Accessor<Other> &other)
      :
      AccessorBase(other)
    {}


    inline
    TrilinosScalar
    Accessor<true>::value() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*value_cache)[a_index];
    }


    inline
    Accessor<false>::Reference::Reference (
      const Accessor<false> &acc)
      :
      accessor(const_cast<Accessor<false>&>(acc))
    {}


    inline
    Accessor<false>::Reference::operator TrilinosScalar () const
    {
      return (*accessor.value_cache)[accessor.a_index];
    }

    inline
    const Accessor<false>::Reference &
    Accessor<false>::Reference::operator = (const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] = n;
      accessor.matrix->set(accessor.row(), accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline
    const Accessor<false>::Reference &
    Accessor<false>::Reference::operator += (const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] += n;
      accessor.matrix->set(accessor.row(), accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline
    const Accessor<false>::Reference &
    Accessor<false>::Reference::operator -= (const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] -= n;
      accessor.matrix->set(accessor.row(), accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline
    const Accessor<false>::Reference &
    Accessor<false>::Reference::operator *= (const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] *= n;
      accessor.matrix->set(accessor.row(), accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline
    const Accessor<false>::Reference &
    Accessor<false>::Reference::operator /= (const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] /= n;
      accessor.matrix->set(accessor.row(), accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline
    Accessor<false>::Accessor (MatrixType *matrix,
                               const size_type  row,
                               const size_type  index)
      :
      AccessorBase(matrix, row, index)
    {}


    inline
    Accessor<false>::Reference
    Accessor<false>::value() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return Reference(*this);
    }



    template <bool Constness>
    inline
    Iterator<Constness>::Iterator(MatrixType *matrix,
                                  const size_type  row,
                                  const size_type  index)
      :
      accessor(matrix, row, index)
    {}


    template <bool Constness>
    template <bool Other>
    inline
    Iterator<Constness>::Iterator(const Iterator<Other> &other)
      :
      accessor(other.accessor)
    {}


    template <bool Constness>
    inline
    Iterator<Constness> &
    Iterator<Constness>::operator++ ()
    {
      Assert (accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

      ++accessor.a_index;

      // If at end of line: do one
      // step, then cycle until we
      // find a row with a nonzero
      // number of entries.
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while ((accessor.a_row < accessor.matrix->m())
                 &&
                 (accessor.matrix->row_length(accessor.a_row) == 0))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }


    template <bool Constness>
    inline
    Iterator<Constness>
    Iterator<Constness>::operator++ (int)
    {
      const Iterator<Constness> old_state = *this;
      ++(*this);
      return old_state;
    }



    template <bool Constness>
    inline
    const Accessor<Constness> &
    Iterator<Constness>::operator* () const
    {
      return accessor;
    }



    template <bool Constness>
    inline
    const Accessor<Constness> *
    Iterator<Constness>::operator-> () const
    {
      return &accessor;
    }



    template <bool Constness>
    inline
    bool
    Iterator<Constness>::operator == (const Iterator<Constness> &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }



    template <bool Constness>
    inline
    bool
    Iterator<Constness>::operator != (const Iterator<Constness> &other) const
    {
      return ! (*this == other);
    }



    template <bool Constness>
    inline
    bool
    Iterator<Constness>::operator < (const Iterator<Constness> &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }


    template <bool Constness>
    inline
    bool
    Iterator<Constness>::operator > (const Iterator<Constness> &other) const
    {
      return (other < *this);
    }

  }



  inline
  SparseMatrix::const_iterator
  SparseMatrix::begin() const
  {
    return const_iterator(this, 0, 0);
  }



  inline
  SparseMatrix::const_iterator
  SparseMatrix::end() const
  {
    return const_iterator(this, m(), 0);
  }



  inline
  SparseMatrix::const_iterator
  SparseMatrix::begin(const size_type r) const
  {
    Assert (r < m(), ExcIndexRange(r, 0, m()));
    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end (r);
  }



  inline
  SparseMatrix::const_iterator
  SparseMatrix::end(const size_type r) const
  {
    Assert (r < m(), ExcIndexRange(r, 0, m()));

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i=r+1; i<m(); ++i)
      if (row_length(i) > 0)
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline
  SparseMatrix::iterator
  SparseMatrix::begin()
  {
    return iterator(this, 0, 0);
  }



  inline
  SparseMatrix::iterator
  SparseMatrix::end()
  {
    return iterator(this, m(), 0);
  }



  inline
  SparseMatrix::iterator
  SparseMatrix::begin(const size_type r)
  {
    Assert (r < m(), ExcIndexRange(r, 0, m()));
    if (row_length(r) > 0)
      return iterator(this, r, 0);
    else
      return end (r);
  }



  inline
  SparseMatrix::iterator
  SparseMatrix::end(const size_type r)
  {
    Assert (r < m(), ExcIndexRange(r, 0, m()));

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i=r+1; i<m(); ++i)
      if (row_length(i) > 0)
        return iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline
  bool
  SparseMatrix::in_local_range (const size_type index) const
  {
    TrilinosWrappers::types::int_type begin, end;
#ifndef DEAL_II_WITH_64BIT_INDICES
    begin = matrix->RowMap().MinMyGID();
    end = matrix->RowMap().MaxMyGID()+1;
#else
    begin = matrix->RowMap().MinMyGID();
    end = matrix->RowMap().MaxMyGID()+1;
#endif

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline
  bool
  SparseMatrix::is_compressed () const
  {
    return compressed;
  }



  inline
  void
  SparseMatrix::compress ()
  {
    compress(::dealii::VectorOperation::unknown);
  }



  inline
  SparseMatrix &
  SparseMatrix::operator = (const double d)
  {
    Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());
    compress (::dealii::VectorOperation::unknown); // TODO: why do we do this? Should we not check for is_compressed?

    const int ierr = matrix->PutScalar(d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    if (nonlocal_matrix.get() != 0)
      nonlocal_matrix->PutScalar(d);

    return *this;
  }



  // Inline the set() and add() functions, since they will be called
  // frequently, and the compiler can optimize away some unnecessary loops
  // when the sizes are given at compile time.
  inline
  void
  SparseMatrix::set (const size_type      i,
                     const size_type      j,
                     const TrilinosScalar value)
  {

    Assert (numbers::is_finite(value), ExcNumberNotFinite());

    set (i, 1, &j, &value, false);
  }



  inline
  void
  SparseMatrix::set (const std::vector<size_type>  &indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (indices.size() == values.m(),
            ExcDimensionMismatch(indices.size(), values.m()));
    Assert (values.m() == values.n(), ExcNotQuadratic());

    for (size_type i=0; i<indices.size(); ++i)
      set (indices[i], indices.size(), &indices[0], &values(i,0),
           elide_zero_values);
  }



  inline
  void
  SparseMatrix::set (const std::vector<size_type>     &row_indices,
                     const std::vector<size_type>     &col_indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (row_indices.size() == values.m(),
            ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert (col_indices.size() == values.n(),
            ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i=0; i<row_indices.size(); ++i)
      set (row_indices[i], col_indices.size(), &col_indices[0], &values(i,0),
           elide_zero_values);
  }



  inline
  void
  SparseMatrix::set (const size_type                    row,
                     const std::vector<size_type>      &col_indices,
                     const std::vector<TrilinosScalar> &values,
                     const bool                         elide_zero_values)
  {
    Assert (col_indices.size() == values.size(),
            ExcDimensionMismatch(col_indices.size(), values.size()));

    set (row, col_indices.size(), &col_indices[0], &values[0],
         elide_zero_values);
  }



  inline
  void
  SparseMatrix::set (const size_type       row,
                     const size_type       n_cols,
                     const size_type      *col_indices,
                     const TrilinosScalar *values,
                     const bool            elide_zero_values)
  {
    AssertIndexRange(row, this->m());

    int ierr;
    if (last_action == Add)
      {
        ierr = matrix->GlobalAssemble (*column_space_map, matrix->RowMap(),
                                       true);

        Assert (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Insert;

    TrilinosWrappers::types::int_type *col_index_ptr;
    TrilinosScalar *col_value_ptr;
    TrilinosWrappers::types::int_type n_columns;

    TrilinosScalar short_val_array[100];
    TrilinosWrappers::types::int_type short_index_array[100];
    std::vector<TrilinosScalar> long_val_array;
    std::vector<TrilinosWrappers::types::int_type> long_index_array;


    // If we don't elide zeros, the pointers are already available... need to
    // cast to non-const pointers as that is the format taken by Trilinos (but
    // we will not modify const data)
    if (elide_zero_values == false)
      {
        col_index_ptr = (TrilinosWrappers::types::int_type *)col_indices;
        col_value_ptr = const_cast<TrilinosScalar *>(values);
        n_columns = n_cols;
      }
    else
      {
        // Otherwise, extract nonzero values in each row and get the
        // respective indices.
        if (n_cols > 100)
          {
            long_val_array.resize(n_cols);
            long_index_array.resize(n_cols);
            col_index_ptr = &long_index_array[0];
            col_value_ptr = &long_val_array[0];
          }
        else
          {
            col_index_ptr = &short_index_array[0];
            col_value_ptr = &short_val_array[0];
          }

        n_columns = 0;
        for (size_type j=0; j<n_cols; ++j)
          {
            const double value = values[j];
            Assert (numbers::is_finite(value), ExcNumberNotFinite());
            if (value != 0)
              {
                col_index_ptr[n_columns] = col_indices[j];
                col_value_ptr[n_columns] = value;
                n_columns++;
              }
          }

        Assert(n_columns <= (TrilinosWrappers::types::int_type)n_cols, ExcInternalError());
      }


    // If the calling matrix owns the row to which we want to insert values,
    // we can directly call the Epetra_CrsMatrix input function, which is much
    // faster than the Epetra_FECrsMatrix function. We distinguish between two
    // cases: the first one is when the matrix is not filled (i.e., it is
    // possible to add new elements to the sparsity pattern), and the second
    // one is when the pattern is already fixed. In the former case, we add
    // the possibility to insert new values, and in the second we just replace
    // data.
    if (row_partitioner().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == true)
      {
        if (matrix->Filled() == false)
          {
            ierr = matrix->Epetra_CrsMatrix::InsertGlobalValues(
                     static_cast<TrilinosWrappers::types::int_type>(row),
                     static_cast<int>(n_columns),const_cast<double *>(col_value_ptr),
                     col_index_ptr);

            // When inserting elements, we do not want to create exceptions in
            // the case when inserting non-local data (since that's what we
            // want to do right now).
            if (ierr > 0)
              ierr = 0;
          }
        else
          ierr = matrix->Epetra_CrsMatrix::ReplaceGlobalValues(row, n_columns,
                                                               col_value_ptr,
                                                               col_index_ptr);
      }
    else
      {
        // When we're at off-processor data, we have to stick with the
        // standard Insert/ReplaceGlobalValues function. Nevertheless, the way
        // we call it is the fastest one (any other will lead to repeated
        // allocation and deallocation of memory in order to call the function
        // we already use, which is very inefficient if writing one element at
        // a time).
        compressed = false;

        if (matrix->Filled() == false)
          {
            ierr = matrix->InsertGlobalValues (1,
                                               (TrilinosWrappers::types::int_type *)&row,
                                               n_columns, col_index_ptr,
                                               &col_value_ptr,
                                               Epetra_FECrsMatrix::ROW_MAJOR);
            if (ierr > 0)
              ierr = 0;
          }
        else
          ierr = matrix->ReplaceGlobalValues (1,
                                              (TrilinosWrappers::types::int_type *)&row,
                                              n_columns, col_index_ptr,
                                              &col_value_ptr,
                                              Epetra_FECrsMatrix::ROW_MAJOR);
        // use the FECrsMatrix facilities for set even in the case when we
        // have explicitly set the off-processor rows because that only works
        // properly when adding elements, not when setting them (since we want
        // to only touch elements that have been set explicitly, and there is
        // no way on the receiving processor to identify them otherwise)
      }

    Assert (ierr <= 0, ExcAccessToNonPresentElement(row, col_index_ptr[0]));
    AssertThrow (ierr >= 0, ExcTrilinosError(ierr));
  }



  inline
  void
  SparseMatrix::add (const size_type      i,
                     const size_type      j,
                     const TrilinosScalar value)
  {
    Assert (numbers::is_finite(value), ExcNumberNotFinite());

    if (value == 0)
      {
        // we have to do checkings on Insert/Add in any case to be consistent
        // with the MPI communication model (see the comments in the
        // documentation of TrilinosWrappers::Vector), but we can save some
        // work if the addend is zero. However, these actions are done in case
        // we pass on to the other function.

        // TODO: fix this (do not run compress here, but fail)
        if (last_action == Insert)
          {
            int ierr;
            ierr = matrix->GlobalAssemble(*column_space_map,
                                          row_partitioner(), false);

            Assert (ierr == 0, ExcTrilinosError(ierr));
            (void)ierr; // removes -Wunused-but-set-variable in optimized mode
          }

        last_action = Add;

        return;
      }
    else
      add (i, 1, &j, &value, false);
  }



  inline
  void
  SparseMatrix::add (const std::vector<size_type>     &indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (indices.size() == values.m(),
            ExcDimensionMismatch(indices.size(), values.m()));
    Assert (values.m() == values.n(), ExcNotQuadratic());

    for (size_type i=0; i<indices.size(); ++i)
      add (indices[i], indices.size(), &indices[0], &values(i,0),
           elide_zero_values);
  }



  inline
  void
  SparseMatrix::add (const std::vector<size_type>  &row_indices,
                     const std::vector<size_type>  &col_indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (row_indices.size() == values.m(),
            ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert (col_indices.size() == values.n(),
            ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i=0; i<row_indices.size(); ++i)
      add (row_indices[i], col_indices.size(), &col_indices[0],
           &values(i,0), elide_zero_values);
  }



  inline
  void
  SparseMatrix::add (const size_type                    row,
                     const std::vector<size_type>      &col_indices,
                     const std::vector<TrilinosScalar> &values,
                     const bool                         elide_zero_values)
  {
    Assert (col_indices.size() == values.size(),
            ExcDimensionMismatch(col_indices.size(), values.size()));

    add (row, col_indices.size(), &col_indices[0], &values[0],
         elide_zero_values);
  }



  inline
  void
  SparseMatrix::add (const size_type       row,
                     const size_type       n_cols,
                     const size_type      *col_indices,
                     const TrilinosScalar *values,
                     const bool            elide_zero_values,
                     const bool            /*col_indices_are_sorted*/)
  {
    AssertIndexRange(row, this->m());
    int ierr;
    if (last_action == Insert)
      {
        // TODO: this could lead to a dead lock when only one processor
        // calls GlobalAssemble.
        ierr = matrix->GlobalAssemble(*column_space_map,
                                      row_partitioner(), false);

        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Add;

    TrilinosWrappers::types::int_type *col_index_ptr;
    TrilinosScalar *col_value_ptr;
    TrilinosWrappers::types::int_type n_columns;

    double short_val_array[100];
    TrilinosWrappers::types::int_type short_index_array[100];
    std::vector<TrilinosScalar> long_val_array;
    std::vector<TrilinosWrappers::types::int_type> long_index_array;

    // If we don't elide zeros, the pointers are already available... need to
    // cast to non-const pointers as that is the format taken by Trilinos (but
    // we will not modify const data)
    if (elide_zero_values == false)
      {
        col_index_ptr = (TrilinosWrappers::types::int_type *)col_indices;
        col_value_ptr = const_cast<TrilinosScalar *>(values);
        n_columns = n_cols;
#ifdef DEBUG
        for (size_type j=0; j<n_cols; ++j)
          Assert (numbers::is_finite(values[j]), ExcNumberNotFinite());
#endif
      }
    else
      {
        // Otherwise, extract nonzero values in each row and the corresponding
        // index.
        if (n_cols > 100)
          {
            long_val_array.resize(n_cols);
            long_index_array.resize(n_cols);
            col_index_ptr = &long_index_array[0];
            col_value_ptr = &long_val_array[0];
          }
        else
          {
            col_index_ptr = &short_index_array[0];
            col_value_ptr = &short_val_array[0];
          }

        n_columns = 0;
        for (size_type j=0; j<n_cols; ++j)
          {
            const double value = values[j];

            Assert (numbers::is_finite(value), ExcNumberNotFinite());
            if (value != 0)
              {
                col_index_ptr[n_columns] = col_indices[j];
                col_value_ptr[n_columns] = value;
                n_columns++;
              }
          }

        Assert(n_columns <= (TrilinosWrappers::types::int_type)n_cols, ExcInternalError());

      }

    // If the calling processor owns the row to which we want to add values, we
    // can directly call the Epetra_CrsMatrix input function, which is much
    // faster than the Epetra_FECrsMatrix function.
    if (row_partitioner().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == true)
      {
        ierr = matrix->Epetra_CrsMatrix::SumIntoGlobalValues(row, n_columns,
                                                             col_value_ptr,
                                                             col_index_ptr);
      }
    else if (nonlocal_matrix.get() != 0)
      {
        compressed = false;
        // this is the case when we have explicitly set the off-processor rows
        // and want to create a separate matrix object for them (to retain
        // thread-safety)
        Assert (nonlocal_matrix->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row)) != -1,
                ExcMessage("Attempted to write into off-processor matrix row "
                           "that has not be specified as being writable upon "
                           "initialization"));
        ierr = nonlocal_matrix->SumIntoGlobalValues(row, n_columns,
                                                    col_value_ptr,
                                                    col_index_ptr);
      }
    else
      {
        // When we're at off-processor data, we have to stick with the
        // standard SumIntoGlobalValues function. Nevertheless, the way we
        // call it is the fastest one (any other will lead to repeated
        // allocation and deallocation of memory in order to call the function
        // we already use, which is very inefficient if writing one element at
        // a time).
        compressed = false;

        ierr = matrix->SumIntoGlobalValues (1,
                                            (TrilinosWrappers::types::int_type *)&row, n_columns,
                                            col_index_ptr,
                                            &col_value_ptr,
                                            Epetra_FECrsMatrix::ROW_MAJOR);
      }

#ifdef DEBUG
    if (ierr > 0)
      {
        std::cout << "------------------------------------------"
                  << std::endl;
        std::cout << "Got error " << ierr << " in row " << row
                  << " of proc " << row_partitioner().Comm().MyPID()
                  << " when trying to add the columns:" << std::endl;
        for (TrilinosWrappers::types::int_type i=0; i<n_columns; ++i)
          std::cout << col_index_ptr[i] << " ";
        std::cout << std::endl << std::endl;
        std::cout << "Matrix row "
                  << (row_partitioner().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == false ? "(nonlocal part)" : "")
                  << " has the following indices:" << std::endl;
        std::vector<TrilinosWrappers::types::int_type> indices;
        const Epetra_CrsGraph *graph =
          (nonlocal_matrix.get() != 0 &&
           row_partitioner().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == false) ?
          &nonlocal_matrix->Graph() : &matrix->Graph();

        indices.resize(graph->NumGlobalIndices(static_cast<TrilinosWrappers::types::int_type>(row)));
        int n_indices = 0;
        graph->ExtractGlobalRowCopy(static_cast<TrilinosWrappers::types::int_type>(row),
                                    indices.size(), n_indices, &indices[0]);
        AssertDimension(static_cast<unsigned int>(n_indices), indices.size());

        for (TrilinosWrappers::types::int_type i=0; i<n_indices; ++i)
          std::cout << indices[i] << " ";
        std::cout << std::endl << std::endl;
        Assert (ierr <= 0,
                ExcAccessToNonPresentElement(row, col_index_ptr[0]));
      }
#endif
    Assert (ierr >= 0, ExcTrilinosError(ierr));
  }



  // inline "simple" functions that are called frequently and do only involve
  // a call to some Trilinos function.
  inline
  SparseMatrix::size_type
  SparseMatrix::m () const
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalRows();
#else
    return matrix->NumGlobalRows64();
#endif
  }



  inline
  SparseMatrix::size_type
  SparseMatrix::n () const
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalCols();
#else
    return matrix->NumGlobalCols64();
#endif
  }



  inline
  unsigned int
  SparseMatrix::local_size () const
  {
    return matrix -> NumMyRows();
  }



  inline
  std::pair<SparseMatrix::size_type, SparseMatrix::size_type>
  SparseMatrix::local_range () const
  {
    size_type begin, end;
#ifndef DEAL_II_WITH_64BIT_INDICES
    begin = matrix->RowMap().MinMyGID();
    end = matrix->RowMap().MaxMyGID()+1;
#else
    begin = matrix->RowMap().MinMyGID64();
    end = matrix->RowMap().MaxMyGID64()+1;
#endif

    return std::make_pair (begin, end);
  }



  inline
  SparseMatrix::size_type
  SparseMatrix::n_nonzero_elements () const
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalNonzeros();
#else
    return matrix->NumGlobalNonzeros64();
#endif
  }



  template <typename SparsityType>
  inline
  void SparseMatrix::reinit (const IndexSet      &parallel_partitioning,
                             const SparsityType  &sparsity_pattern,
                             const MPI_Comm      &communicator,
                             const bool           exchange_data)
  {
    Epetra_Map map = parallel_partitioning.make_trilinos_map (communicator, false);
    reinit (map, map, sparsity_pattern, exchange_data);
  }



  template <typename SparsityType>
  inline
  void SparseMatrix::reinit (const IndexSet      &row_parallel_partitioning,
                             const IndexSet      &col_parallel_partitioning,
                             const SparsityType  &sparsity_pattern,
                             const MPI_Comm      &communicator,
                             const bool           exchange_data)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit (row_map, col_map, sparsity_pattern, exchange_data);
  }



  template <typename number>
  inline
  void SparseMatrix::reinit (const IndexSet      &parallel_partitioning,
                             const ::dealii::SparseMatrix<number> &sparse_matrix,
                             const MPI_Comm      &communicator,
                             const double         drop_tolerance,
                             const bool           copy_values,
                             const ::dealii::SparsityPattern *use_this_sparsity)
  {
    Epetra_Map map = parallel_partitioning.make_trilinos_map (communicator, false);
    reinit (map, map, sparse_matrix, drop_tolerance, copy_values,
            use_this_sparsity);
  }



  template <typename number>
  inline
  void SparseMatrix::reinit (const IndexSet      &row_parallel_partitioning,
                             const IndexSet      &col_parallel_partitioning,
                             const ::dealii::SparseMatrix<number> &sparse_matrix,
                             const MPI_Comm      &communicator,
                             const double         drop_tolerance,
                             const bool           copy_values,
                             const ::dealii::SparsityPattern *use_this_sparsity)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit (row_map, col_map, sparse_matrix, drop_tolerance, copy_values,
            use_this_sparsity);
  }



  inline
  TrilinosScalar
  SparseMatrix::l1_norm () const
  {
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormOne();
  }



  inline
  TrilinosScalar
  SparseMatrix::linfty_norm () const
  {
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormInf();
  }



  inline
  TrilinosScalar
  SparseMatrix::frobenius_norm () const
  {
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormFrobenius();
  }



  inline
  SparseMatrix &
  SparseMatrix::operator *= (const TrilinosScalar a)
  {
    const int ierr = matrix->Scale (a);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode

    return *this;
  }



  inline
  SparseMatrix &
  SparseMatrix::operator /= (const TrilinosScalar a)
  {
    Assert (a !=0, ExcDivideByZero());

    const TrilinosScalar factor = 1./a;

    const int ierr = matrix->Scale (factor);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode

    return *this;
  }



  namespace internal
  {
    namespace SparseMatrix
    {
      template <typename VectorType>
      inline
      void check_vector_map_equality(const Epetra_CrsMatrix &,
                                     const VectorType &,
                                     const VectorType &)
      {
      }

      inline
      void check_vector_map_equality(const Epetra_CrsMatrix              &m,
                                     const TrilinosWrappers::MPI::Vector &in,
                                     const TrilinosWrappers::MPI::Vector &out)
      {
        Assert (in.vector_partitioner().SameAs(m.DomainMap()) == true,
                ExcMessage ("Column map of matrix does not fit with vector map!"));
        Assert (out.vector_partitioner().SameAs(m.RangeMap()) == true,
                ExcMessage ("Row map of matrix does not fit with vector map!"));
        (void)m;
        (void)in;
        (void)out;
      }
    }
  }


  template <typename VectorType>
  inline
  void
  SparseMatrix::vmult (VectorType       &dst,
                       const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    (void)src;
    (void)dst;

    internal::SparseMatrix::check_vector_map_equality(*matrix, src, dst);
    const size_type dst_local_size = dst.end() - dst.begin();
    AssertDimension (dst_local_size, static_cast<size_type>(matrix->RangeMap().NumMyElements()));
    (void)dst_local_size;
    const size_type src_local_size = src.end() - src.begin();
    AssertDimension (src_local_size, static_cast<size_type>(matrix->DomainMap().NumMyElements()));
    (void)src_local_size;

    Epetra_MultiVector tril_dst (View, matrix->RangeMap(), dst.begin(),
                                 matrix->DomainMap().NumMyPoints(), 1);
    Epetra_MultiVector tril_src (View, matrix->DomainMap(),
                                 const_cast<TrilinosScalar *>(src.begin()),
                                 matrix->DomainMap().NumMyPoints(), 1);

    const int ierr = matrix->Multiply (false, tril_src, tril_dst);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode
  }



  template <typename VectorType>
  inline
  void
  SparseMatrix::Tvmult (VectorType       &dst,
                        const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());
    Assert (matrix->Filled(), ExcMatrixNotCompressed());

    internal::SparseMatrix::check_vector_map_equality(*matrix, dst, src);
    const size_type dst_local_size = dst.end() - dst.begin();
    AssertDimension (dst_local_size, static_cast<size_type>(matrix->DomainMap().NumMyElements()));
    const size_type src_local_size = src.end() - src.begin();
    AssertDimension (src_local_size, static_cast<size_type>(matrix->RangeMap().NumMyElements()));

    Epetra_MultiVector tril_dst (View, matrix->DomainMap(), dst.begin(),
                                 matrix->DomainMap().NumMyPoints(), 1);
    Epetra_MultiVector tril_src (View, matrix->RangeMap(),
                                 const_cast<double *>(src.begin()),
                                 matrix->DomainMap().NumMyPoints(), 1);

    const int ierr = matrix->Multiply (true, tril_src, tril_dst);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode
  }



  template <typename VectorType>
  inline
  void
  SparseMatrix::vmult_add (VectorType       &dst,
                           const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    // Reinit a temporary vector with fast argument set, which does not
    // overwrite the content (to save time). However, the
    // TrilinosWrappers::Vector classes do not support this, so create a
    // deal.II local vector that has this fast setting. It will be accepted in
    // vmult because it only checks the local size.
    dealii::Vector<TrilinosScalar> temp_vector;
    temp_vector.reinit(dst.end()-dst.begin(), true);
    dealii::VectorView<TrilinosScalar> src_view(src.end()-src.begin(), src.begin());
    dealii::VectorView<TrilinosScalar> dst_view(dst.end()-dst.begin(), dst.begin());
    vmult (temp_vector, static_cast<const dealii::Vector<TrilinosScalar>&>(src_view));
    if (dst_view.size() > 0)
      dst_view += temp_vector;
  }



  template <typename VectorType>
  inline
  void
  SparseMatrix::Tvmult_add (VectorType       &dst,
                            const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    // Reinit a temporary vector with fast argument set, which does not
    // overwrite the content (to save time). However, the
    // TrilinosWrappers::Vector classes do not support this, so create a
    // deal.II local vector that has this fast setting. It will be accepted in
    // vmult because it only checks the local size.
    dealii::Vector<TrilinosScalar> temp_vector;
    temp_vector.reinit(dst.end()-dst.begin(), true);
    dealii::VectorView<TrilinosScalar> src_view(src.end()-src.begin(), src.begin());
    dealii::VectorView<TrilinosScalar> dst_view(dst.end()-dst.begin(), dst.begin());
    Tvmult (temp_vector, static_cast<const dealii::Vector<TrilinosScalar>&>(src_view));
    if (dst_view.size() > 0)
      dst_view += temp_vector;
  }



  inline
  TrilinosScalar
  SparseMatrix::matrix_norm_square (const VectorBase &v) const
  {
    Assert (row_partitioner().SameAs(domain_partitioner()),
            ExcNotQuadratic());

    VectorBase temp_vector;
    temp_vector.reinit(v, true);

    vmult (temp_vector, v);
    return temp_vector*v;
  }



  inline
  TrilinosScalar
  SparseMatrix::matrix_scalar_product (const VectorBase &u,
                                       const VectorBase &v) const
  {
    Assert (row_partitioner().SameAs(domain_partitioner()),
            ExcNotQuadratic());

    VectorBase temp_vector;
    temp_vector.reinit(v, true);

    vmult (temp_vector, v);
    return u*temp_vector;
  }



  inline
  TrilinosScalar
  SparseMatrix::residual (VectorBase       &dst,
                          const VectorBase &x,
                          const VectorBase &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }


  inline
  const Epetra_CrsMatrix &
  SparseMatrix::trilinos_matrix () const
  {
    return static_cast<const Epetra_CrsMatrix &>(*matrix);
  }



  inline
  const Epetra_CrsGraph &
  SparseMatrix::trilinos_sparsity_pattern () const
  {
    return matrix->Graph();
  }



  inline
  const Epetra_Map &
  SparseMatrix::domain_partitioner () const
  {
    return matrix->DomainMap();
  }



  inline
  const Epetra_Map &
  SparseMatrix::range_partitioner () const
  {
    return matrix->RangeMap();
  }



  inline
  const Epetra_Map &
  SparseMatrix::row_partitioner () const
  {
    return matrix->RowMap();
  }



  inline
  const Epetra_Map &
  SparseMatrix::col_partitioner () const
  {
    return matrix->ColMap();
  }



  inline
  void
  SparseMatrix::prepare_add()
  {
    //nothing to do here
  }



  inline
  void
  SparseMatrix::prepare_set()
  {
    //nothing to do here
  }



#endif // DOXYGEN

}


DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_TRILINOS


/*-----------------------   trilinos_sparse_matrix.h     --------------------*/

#endif
/*-----------------------   trilinos_sparse_matrix.h     --------------------*/
