// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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

#ifndef dealii_petsc_matrix_base_h
#  define dealii_petsc_matrix_base_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/full_matrix.h>
#    include <deal.II/lac/petsc_compatibility.h>
#    include <deal.II/lac/petsc_vector_base.h>
#    include <deal.II/lac/vector_operation.h>

#    include <petscmat.h>

#    include <cmath>
#    include <memory>
#    include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#    ifndef DOXYGEN
template <typename Matrix>
class BlockMatrixBase;
#    endif


namespace PETScWrappers
{
  // forward declarations
  class MatrixBase;

  namespace MatrixIterators
  {
    /**
     * This class acts as an iterator walking over the elements of PETSc
     * matrices. Since PETSc offers a uniform interface for all types of
     * matrices, this iterator can be used to access both sparse and full
     * matrices.
     *
     * Note that PETSc does not give any guarantees as to the order of
     * elements within each row. Note also that accessing the elements of a
     * full matrix surprisingly only shows the nonzero elements of the matrix,
     * not all elements.
     *
     * @ingroup PETScWrappers
     * @author Guido Kanschat, Roy Stogner, Wolfgang Bangerth, 2004
     */
    class const_iterator
    {
    private:
      /**
       * Accessor class for iterators
       */
      class Accessor
      {
      public:
        /**
         * Declare type for container size.
         */
        using size_type = types::global_dof_index;

        /**
         * Constructor. Since we use accessors only for read access, a const
         * matrix pointer is sufficient.
         */
        Accessor(const MatrixBase *matrix,
                 const size_type   row,
                 const size_type   index);

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

        /**
         * Value of this matrix entry.
         */
        PetscScalar
        value() const;

        /**
         * Exception
         */
        DeclException0(ExcBeyondEndOfMatrix);
        /**
         * Exception
         */
        DeclException3(ExcAccessToNonlocalRow,
                       int,
                       int,
                       int,
                       << "You tried to access row " << arg1
                       << " of a distributed matrix, but only rows " << arg2
                       << " through " << arg3
                       << " are stored locally and can be accessed.");

      private:
        /**
         * The matrix accessed.
         */
        mutable MatrixBase *matrix;

        /**
         * Current row number.
         */
        size_type a_row;

        /**
         * Current index in row.
         */
        size_type a_index;

        /**
         * Cache where we store the column indices of the present row. This is
         * necessary, since PETSc makes access to the elements of its matrices
         * rather hard, and it is much more efficient to copy all column
         * entries of a row once when we enter it than repeatedly asking PETSc
         * for individual ones. This also makes some sense since it is likely
         * that we will access them sequentially anyway.
         *
         * In order to make copying of iterators/accessor of acceptable
         * performance, we keep a shared pointer to these entries so that more
         * than one accessor can access this data if necessary.
         */
        std::shared_ptr<const std::vector<size_type>> colnum_cache;

        /**
         * Similar cache for the values of this row.
         */
        std::shared_ptr<const std::vector<PetscScalar>> value_cache;

        /**
         * Discard the old row caches (they may still be used by other
         * accessors) and generate new ones for the row pointed to presently
         * by this accessor.
         */
        void
        visit_present_row();

        // Make enclosing class a friend.
        friend class const_iterator;
      };

    public:
      /**
       * Declare type for container size.
       */
      using size_type = types::global_dof_index;

      /**
       * Constructor. Create an iterator into the matrix @p matrix for the
       * given row and the index within it.
       */
      const_iterator(const MatrixBase *matrix,
                     const size_type   row,
                     const size_type   index);

      /**
       * Prefix increment.
       */
      const_iterator &
      operator++();

      /**
       * Postfix increment.
       */
      const_iterator
      operator++(int);

      /**
       * Dereferencing operator.
       */
      const Accessor &operator*() const;

      /**
       * Dereferencing operator.
       */
      const Accessor *operator->() const;

      /**
       * Comparison. True, if both iterators point to the same matrix
       * position.
       */
      bool
      operator==(const const_iterator &) const;
      /**
       * Inverse of <tt>==</tt>.
       */
      bool
      operator!=(const const_iterator &) const;

      /**
       * Comparison operator. Result is true if either the first row number is
       * smaller or if the row numbers are equal and the first index is
       * smaller.
       */
      bool
      operator<(const const_iterator &) const;

      /**
       * Exception
       */
      DeclException2(ExcInvalidIndexWithinRow,
                     int,
                     int,
                     << "Attempt to access element " << arg2 << " of row "
                     << arg1 << " which doesn't have that many elements.");

    private:
      /**
       * Store an object of the accessor class.
       */
      Accessor accessor;
    };

  } // namespace MatrixIterators


  /**
   * Base class for all matrix classes that are implemented on top of the
   * PETSc matrix types. Since in PETSc all matrix types (i.e. sequential and
   * parallel, sparse, blocked, etc.)  are built by filling the contents of an
   * abstract object that is only referenced through a pointer of a type that
   * is independent of the actual matrix type, we can implement almost all
   * functionality of matrices in this base class. Derived classes will then
   * only have to provide the functionality to create one or the other kind of
   * matrix.
   *
   * The interface of this class is modeled after the existing SparseMatrix
   * class in deal.II. It has almost the same member functions, and is often
   * exchangeable. However, since PETSc only supports a single scalar type
   * (either double, float, or a complex data type), it is not templated, and
   * only works with whatever your PETSc installation has defined the data
   * type PetscScalar to.
   *
   * Note that PETSc only guarantees that operations do what you expect if the
   * functions @p MatAssemblyBegin and @p MatAssemblyEnd have been called
   * after matrix assembly. Therefore, you need to call
   * SparseMatrix::compress() before you actually use the matrix. This also
   * calls @p MatCompress that compresses the storage format for sparse
   * matrices by discarding unused elements. PETSc allows to continue with
   * assembling the matrix after calls to these functions, but since there are
   * no more free entries available after that any more, it is better to only
   * call SparseMatrix::compress() once at the end of the assembly stage and
   * before the matrix is actively used.
   *
   * @ingroup PETScWrappers
   * @ingroup Matrix1
   * @author Wolfgang Bangerth, 2004
   */
  class MatrixBase : public Subscriptor
  {
  public:
    /**
     * Declare an alias for the iterator class.
     */
    using const_iterator = MatrixIterators::const_iterator;

    /**
     * Declare type for container size.
     */
    using size_type = types::global_dof_index;

    /**
     * Declare an alias in analogy to all the other container classes.
     */
    using value_type = PetscScalar;

    /**
     * Default constructor.
     */
    MatrixBase();

    /**
     * Copy constructor. It is deleted as copying this base class
     * without knowing the concrete kind of matrix stored may both
     * miss important details and be expensive if the matrix is large.
     */
    MatrixBase(const MatrixBase &) = delete;

    /**
     * Copy operator. It is deleted as copying this base class
     * without knowing the concrete kind of matrix stored may both
     * miss important details and be expensive if the matrix is large.
     */
    MatrixBase &
    operator=(const MatrixBase &) = delete;

    /**
     * Destructor. Made virtual so that one can use pointers to this class.
     */
    virtual ~MatrixBase() override;

    /**
     * This operator assigns a scalar to a matrix. Since this does usually not
     * make much sense (should we set all matrix entries to this value? Only
     * the nonzero entries of the sparsity pattern?), this operation is only
     * allowed if the actual value to be assigned is zero. This operator only
     * exists to allow for the obvious notation <tt>matrix=0</tt>, which sets
     * all elements of the matrix to zero, but keeps the sparsity pattern
     * previously used.
     */
    MatrixBase &
    operator=(const value_type d);
    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     */
    void
    clear();

    /**
     * Set the element (<i>i,j</i>) to @p value.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds a new entry to the matrix if
     * it didn't exist before, very much in contrast to the SparseMatrix class
     * which throws an error if the entry does not exist. If <tt>value</tt> is
     * not a finite number an exception is thrown.
     */
    void
    set(const size_type i, const size_type j, const PetscScalar value);

    /**
     * Set all elements given in a FullMatrix<double> into the sparse matrix
     * locations given by <tt>indices</tt>. In other words, this function
     * writes the elements in <tt>full_matrix</tt> into the calling matrix,
     * using the local-to-global indexing specified by <tt>indices</tt> for
     * both the rows and the columns of the matrix. This function assumes a
     * quadratic sparse matrix and a quadratic full_matrix, the usual
     * situation in FE calculations.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds some new entries to the matrix
     * if they didn't exist before, very much in contrast to the SparseMatrix
     * class which throws an error if the entry does not exist.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     */
    void
    set(const std::vector<size_type> & indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = false);

    /**
     * Same function as before, but now including the possibility to use
     * rectangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     */
    void
    set(const std::vector<size_type> & row_indices,
        const std::vector<size_type> & col_indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = false);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by <tt>col_indices</tt> to the respective value.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds some new entries to the matrix
     * if they didn't exist before, very much in contrast to the SparseMatrix
     * class which throws an error if the entry does not exist.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     */
    void
    set(const size_type                 row,
        const std::vector<size_type> &  col_indices,
        const std::vector<PetscScalar> &values,
        const bool                      elide_zero_values = false);

    /**
     * Set several elements to values given by <tt>values</tt> in a given row
     * in columns given by col_indices into the sparse matrix.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds some new entries to the matrix
     * if they didn't exist before, very much in contrast to the SparseMatrix
     * class which throws an error if the entry does not exist.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be inserted anyway or they should be
     * filtered away. The default value is <tt>false</tt>, i.e., even zero
     * values are inserted/replaced.
     */
    void
    set(const size_type    row,
        const size_type    n_cols,
        const size_type *  col_indices,
        const PetscScalar *values,
        const bool         elide_zero_values = false);

    /**
     * Add @p value to the element (<i>i,j</i>).
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds a new entry to the matrix if
     * it didn't exist before, very much in contrast to the SparseMatrix class
     * which throws an error if the entry does not exist. If <tt>value</tt> is
     * not a finite number an exception is thrown.
     */
    void
    add(const size_type i, const size_type j, const PetscScalar value);

    /**
     * Add all elements given in a FullMatrix<double> into sparse matrix
     * locations given by <tt>indices</tt>. In other words, this function adds
     * the elements in <tt>full_matrix</tt> to the respective entries in
     * calling matrix, using the local-to-global indexing specified by
     * <tt>indices</tt> for both the rows and the columns of the matrix. This
     * function assumes a quadratic sparse matrix and a quadratic full_matrix,
     * the usual situation in FE calculations.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds some new entries to the matrix
     * if they didn't exist before, very much in contrast to the SparseMatrix
     * class which throws an error if the entry does not exist.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     */
    void
    add(const std::vector<size_type> & indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = true);

    /**
     * Same function as before, but now including the possibility to use
     * rectangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     */
    void
    add(const std::vector<size_type> & row_indices,
        const std::vector<size_type> & col_indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = true);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by <tt>col_indices</tt> to the respective value.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds some new entries to the matrix
     * if they didn't exist before, very much in contrast to the SparseMatrix
     * class which throws an error if the entry does not exist.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     */
    void
    add(const size_type                 row,
        const std::vector<size_type> &  col_indices,
        const std::vector<PetscScalar> &values,
        const bool                      elide_zero_values = true);

    /**
     * Add an array of values given by <tt>values</tt> in the given global
     * matrix row at columns specified by col_indices in the sparse matrix.
     *
     * If the present object (from a derived class of this one) happens to be
     * a sparse matrix, then this function adds some new entries to the matrix
     * if they didn't exist before, very much in contrast to the SparseMatrix
     * class which throws an error if the entry does not exist.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     */
    void
    add(const size_type    row,
        const size_type    n_cols,
        const size_type *  col_indices,
        const PetscScalar *values,
        const bool         elide_zero_values      = true,
        const bool         col_indices_are_sorted = false);

    /**
     * Remove all elements from this <tt>row</tt> by setting them to zero. The
     * function does not modify the number of allocated nonzero entries, it
     * only sets some entries to zero. It may drop them from the sparsity
     * pattern, though (but retains the allocated memory in case new entries
     * are again added later).
     *
     * This operation is used in eliminating constraints (e.g. due to hanging
     * nodes) and makes sure that we can write this modification to the matrix
     * without having to read entries (such as the locations of non-zero
     * elements) from it -- without this operation, removing constraints on
     * parallel matrices is a rather complicated procedure.
     *
     * The second parameter can be used to set the diagonal entry of this row
     * to a value different from zero. The default is to set it to zero.
     */
    void
    clear_row(const size_type row, const PetscScalar new_diag_value = 0);

    /**
     * Same as clear_row(), except that it works on a number of rows at once.
     *
     * The second parameter can be used to set the diagonal entries of all
     * cleared rows to something different from zero. Note that all of these
     * diagonal entries get the same value -- if you want different values for
     * the diagonal entries, you have to set them by hand.
     */
    void
    clear_rows(const std::vector<size_type> &rows,
               const PetscScalar             new_diag_value = 0);

    /**
     * PETSc matrices store their own sparsity patterns. So, in analogy to our
     * own SparsityPattern class, this function compresses the sparsity
     * pattern and allows the resulting matrix to be used in all other
     * operations where before only assembly functions were allowed. This
     * function must therefore be called once you have assembled the matrix.
     *
     * See
     * @ref GlossCompress "Compressing distributed objects"
     * for more information.
     */
    void
    compress(const VectorOperation::values operation);

    /**
     * Return the value of the entry (<i>i,j</i>).  This may be an expensive
     * operation and you should always take care where to call this function.
     * In contrast to the respective function in the @p MatrixBase class, we
     * don't throw an exception if the respective entry doesn't exist in the
     * sparsity pattern of this class, since PETSc does not transmit this
     * information.
     *
     * This function is therefore exactly equivalent to the <tt>el()</tt>
     * function.
     */
    PetscScalar
    operator()(const size_type i, const size_type j) const;

    /**
     * Return the value of the matrix entry (<i>i,j</i>). If this entry does
     * not exist in the sparsity pattern, then zero is returned. While this
     * may be convenient in some cases, note that it is simple to write
     * algorithms that are slow compared to an optimal solution, since the
     * sparsity of the matrix is not used.
     */
    PetscScalar
    el(const size_type i, const size_type j) const;

    /**
     * Return the main diagonal element in the <i>i</i>th row. This function
     * throws an error if the matrix is not quadratic.
     *
     * Since we do not have direct access to the underlying data structure,
     * this function is no faster than the elementwise access using the el()
     * function. However, we provide this function for compatibility with the
     * SparseMatrix class.
     */
    PetscScalar
    diag_element(const size_type i) const;

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
     * is the same as m(), but for parallel matrices it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range().
     */
    size_type
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
     * Return a reference to the MPI communicator object in use with this
     * matrix. This function has to be implemented in derived classes.
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const = 0;

    /**
     * Return the number of nonzero elements of this matrix. Actually, it
     * returns the number of entries in the sparsity pattern; if any of the
     * entries should happen to be zero, it is counted anyway.
     */
    size_type
    n_nonzero_elements() const;

    /**
     * Number of entries in a specific row.
     */
    size_type
    row_length(const size_type row) const;

    /**
     * Return the l1-norm of the matrix, that is $|M|_1=max_{all columns
     * j}\sum_{all rows i} |M_ij|$, (max. sum of columns). This is the natural
     * matrix norm that is compatible to the l1-norm for vectors, i.e.
     * $|Mv|_1\leq |M|_1 |v|_1$. (cf. Haemmerlin-Hoffmann: Numerische
     * Mathematik)
     */
    PetscReal
    l1_norm() const;

    /**
     * Return the linfty-norm of the matrix, that is $|M|_infty=max_{all rows
     * i}\sum_{all columns j} |M_ij|$, (max. sum of rows). This is the natural
     * matrix norm that is compatible to the linfty-norm of vectors, i.e.
     * $|Mv|_infty \leq |M|_infty |v|_infty$. (cf. Haemmerlin-Hoffmann:
     * Numerische Mathematik)
     */
    PetscReal
    linfty_norm() const;

    /**
     * Return the frobenius norm of the matrix, i.e. the square root of the
     * sum of squares of all entries in the matrix.
     */
    PetscReal
    frobenius_norm() const;


    /**
     * Return the square of the norm of the vector $v$ with respect to the
     * norm induced by this matrix, i.e. $\left(v,Mv\right)$. This is useful,
     * e.g. in the finite element context, where the $L_2$ norm of a function
     * equals the matrix norm with respect to the mass matrix of the vector
     * representing the nodal values of the finite element function.
     *
     * Obviously, the matrix needs to be quadratic for this operation.
     *
     * The implementation of this function is not as efficient as the one in
     * the @p MatrixBase class used in deal.II (i.e. the original one, not the
     * PETSc wrapper class) since PETSc doesn't support this operation and
     * needs a temporary vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then the given
     * vector has to be a distributed vector as well. Conversely, if the
     * matrix is not distributed, then neither may the vector be.
     */
    PetscScalar
    matrix_norm_square(const VectorBase &v) const;


    /**
     * Compute the matrix scalar product $\left(u,Mv\right)$.
     *
     * The implementation of this function is not as efficient as the one in
     * the @p MatrixBase class used in deal.II (i.e. the original one, not the
     * PETSc wrapper class) since PETSc doesn't support this operation and
     * needs a temporary vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then both vectors
     * have to be distributed vectors as well. Conversely, if the matrix is
     * not distributed, then neither of the vectors may be.
     */
    PetscScalar
    matrix_scalar_product(const VectorBase &u, const VectorBase &v) const;

    /**
     * Return the trace of the matrix, i.e. the sum of all diagonal entries in
     * the matrix.
     */
    PetscScalar
    trace() const;

    /**
     * Multiply the entire matrix by a fixed factor.
     */
    MatrixBase &
    operator*=(const PetscScalar factor);

    /**
     * Divide the entire matrix by a fixed factor.
     */
    MatrixBase &
    operator/=(const PetscScalar factor);


    /**
     * Add the matrix @p other scaled by the factor @p factor to the current
     * matrix.
     */
    MatrixBase &
    add(const PetscScalar factor, const MatrixBase &other);


    /**
     * Add the matrix @p other scaled by the factor @p factor to the current
     * matrix.
     * @deprecated Use the function with order of arguments reversed instead.
     */
    DEAL_II_DEPRECATED
    MatrixBase &
    add(const MatrixBase &other, const PetscScalar factor);

    /**
     * Matrix-vector multiplication: let <i>dst = M*src</i> with <i>M</i>
     * being this matrix.
     *
     * Source and destination must not be the same vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then both vectors
     * have to be distributed vectors as well. Conversely, if the matrix is
     * not distributed, then neither of the vectors may be.
     */
    void
    vmult(VectorBase &dst, const VectorBase &src) const;

    /**
     * Matrix-vector multiplication: let <i>dst = M<sup>T</sup>*src</i> with
     * <i>M</i> being this matrix. This function does the same as vmult() but
     * takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then both vectors
     * have to be distributed vectors as well. Conversely, if the matrix is
     * not distributed, then neither of the vectors may be.
     */
    void
    Tvmult(VectorBase &dst, const VectorBase &src) const;

    /**
     * Adding Matrix-vector multiplication. Add <i>M*src</i> on <i>dst</i>
     * with <i>M</i> being this matrix.
     *
     * Source and destination must not be the same vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then both vectors
     * have to be distributed vectors as well. Conversely, if the matrix is
     * not distributed, then neither of the vectors may be.
     */
    void
    vmult_add(VectorBase &dst, const VectorBase &src) const;

    /**
     * Adding Matrix-vector multiplication. Add <i>M<sup>T</sup>*src</i> to
     * <i>dst</i> with <i>M</i> being this matrix. This function does the same
     * as vmult_add() but takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then both vectors
     * have to be distributed vectors as well. Conversely, if the matrix is
     * not distributed, then neither of the vectors may be.
     */
    void
    Tvmult_add(VectorBase &dst, const VectorBase &src) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * Source <i>x</i> and destination <i>dst</i> must not be the same vector.
     *
     * Note that if the current object represents a parallel distributed
     * matrix (of type PETScWrappers::MPI::SparseMatrix), then all vectors
     * have to be distributed vectors as well. Conversely, if the matrix is
     * not distributed, then neither of the vectors may be.
     */
    PetscScalar
    residual(VectorBase &dst, const VectorBase &x, const VectorBase &b) const;

    /**
     * Iterator starting at the first entry. This can only be called on a
     * processor owning the entire matrix. In all other cases refer to the
     * version of begin() taking a row number as an argument.
     */
    const_iterator
    begin() const;

    /**
     * Final iterator. This can only be called on a processor owning the entire
     * matrix. In all other cases refer to the version of end() taking a row
     * number as an argument.
     */
    const_iterator
    end() const;

    /**
     * Iterator starting at the first entry of row @p r.
     *
     * Note that if the given row is empty, i.e. does not contain any nonzero
     * entries, then the iterator returned by this function equals
     * <tt>end(r)</tt>. Note also that the iterator may not be dereferenceable
     * in that case.
     */
    const_iterator
    begin(const size_type r) const;

    /**
     * Final iterator of row <tt>r</tt>. It points to the first element past
     * the end of line @p r, or past the end of the entire sparsity pattern.
     *
     * Note that the end iterator is not necessarily dereferenceable. This is
     * in particular the case if it is the end iterator for the last row of a
     * matrix.
     */
    const_iterator
    end(const size_type r) const;

    /**
     * Conversion operator to gain access to the underlying PETSc type. If you
     * do this, you cut this class off some information it may need, so this
     * conversion operator should only be used if you know what you do. In
     * particular, it should only be used for read-only operations into the
     * matrix.
     */
    operator Mat() const;

    /**
     * Return a reference to the underlying PETSc type. It can be used to
     * modify the underlying data, so use it only when you know what you
     * are doing.
     */
    Mat &
    petsc_matrix();

    /**
     * Make an in-place transpose of a matrix.
     */
    void
    transpose();

    /**
     * Test whether a matrix is symmetric.  Default tolerance is
     * $1000\times32$-bit machine precision.
     */
    PetscBool
    is_symmetric(const double tolerance = 1.e-12);

    /**
     * Test whether a matrix is Hermitian, i.e. it is the complex conjugate of
     * its transpose. Default tolerance is $1000\times32$-bit machine
     * precision.
     */
    PetscBool
    is_hermitian(const double tolerance = 1.e-12);

    /**
     * Print the PETSc matrix object values using PETSc internal matrix viewer
     * function <tt>MatView</tt>. The default format prints the non- zero
     * matrix elements. For other valid view formats, consult
     * http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatView.html
     */
    void
    write_ascii(const PetscViewerFormat format = PETSC_VIEWER_DEFAULT);

    /**
     * Print the elements of a matrix to the given output stream.
     *
     * @param[in,out] out The output stream to which to write.
     * @param[in] alternative_output This argument is ignored. It exists for
     * compatibility with similar functions in other matrix classes.
     */
    void
    print(std::ostream &out, const bool alternative_output = false) const;

    /**
     * Return the number bytes consumed by this matrix on this CPU.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Exception
     */
    DeclExceptionMsg(ExcSourceEqualsDestination,
                     "You are attempting an operation on two matrices that "
                     "are the same object, but the operation requires that the "
                     "two objects are in fact different.");

    /**
     * Exception.
     */
    DeclException2(ExcWrongMode,
                   int,
                   int,
                   << "You tried to do a "
                   << (arg1 == 1 ? "'set'" : (arg1 == 2 ? "'add'" : "???"))
                   << " operation but the matrix is currently in "
                   << (arg2 == 1 ? "'set'" : (arg2 == 2 ? "'add'" : "???"))
                   << " mode. You first have to call 'compress()'.");

  protected:
    /**
     * A generic matrix object in PETSc. The actual type, a sparse matrix, is
     * set in the constructor.
     */
    Mat matrix;

    /**
     * Store whether the last action was a write or add operation.
     */
    VectorOperation::values last_action;

    /**
     * Ensure that the add/set mode that is required for actions following
     * this call is compatible with the current mode. Should be called from
     * all internal functions accessing matrix elements.
     */
    void
    prepare_action(const VectorOperation::values new_action);

    /**
     * Internal function that checks that there are no pending insert/add
     * operations. Throws an exception otherwise. Useful before calling any
     * PETSc internal functions modifying the matrix.
     */
    void
    assert_is_compressed();

    /**
     * For some matrix storage formats, in particular for the PETSc
     * distributed blockmatrices, set and add operations on individual
     * elements can not be freely mixed. Rather, one has to synchronize
     * operations when one wants to switch from setting elements to adding to
     * elements. BlockMatrixBase automatically synchronizes the access by
     * calling this helper function for each block. This function ensures that
     * the matrix is in a state that allows adding elements; if it previously
     * already was in this state, the function does nothing.
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

    /**
     * Base function to perform the matrix-matrix multiplication $C = AB$,
     * or, if a vector $V$ whose size is compatible with B is given,
     * $C = A \text{diag}(V) B$, where $\text{diag}(V)$ defines a
     * diagonal matrix with the vector entries.
     *
     * This function assumes that the calling matrix $A$ and $B$
     * have compatible sizes. The size of $C$ will be set within this
     * function.
     *
     * The content as well as the sparsity pattern of the matrix $C$ will be
     * reset by this function, so make sure that the sparsity pattern is not
     * used somewhere else in your program. This is an expensive operation, so
     * think twice before you use this function.
     */
    void
    mmult(MatrixBase &C, const MatrixBase &B, const VectorBase &V) const;

    /**
     * Base function to perform the matrix-matrix multiplication with
     * the transpose of <tt>this</tt>, i.e., $C = A^T B$, or,
     * if an optional vector $V$ whose size is compatible with $B$ is given,
     * $C = A^T \text{diag}(V) B$, where $\text{diag}(V)$ defines a
     * diagonal matrix with the vector entries.
     *
     * This function assumes that the calling matrix $A$ and $B$
     * have compatible sizes. The size of $C$ will be set within this
     * function.
     *
     * The content as well as the sparsity pattern of the matrix $C$ will be
     * changed by this function, so make sure that the sparsity pattern is not
     * used somewhere else in your program. This is an expensive operation, so
     * think twice before you use this function.
     */
    void
    Tmmult(MatrixBase &C, const MatrixBase &B, const VectorBase &V) const;

  private:
    /**
     * An internal array of integer values that is used to store the column
     * indices when adding/inserting local data into the (large) sparse
     * matrix.
     *
     * This variable does not store any "state" of the matrix
     * object. Rather, it is only used as a temporary buffer by some
     * of the member functions of this class. As with all @p mutable
     * member variables, the use of this variable is not thread-safe
     * unless guarded by a mutex. However, since PETSc matrix
     * operations are not thread-safe anyway, there is no need to
     * attempt to make things thread-safe, and so there is no mutex
     * associated with this variable.
     */
    mutable std::vector<PetscInt> column_indices;

    /**
     * An internal array of double values that is used to store the column
     * indices when adding/inserting local data into the (large) sparse
     * matrix.
     *
     * The same comment as for the @p column_indices variable above
     * applies.
     */
    mutable std::vector<PetscScalar> column_values;


    // To allow calling protected prepare_add() and prepare_set().
    template <class>
    friend class dealii::BlockMatrixBase;
  };



#    ifndef DOXYGEN
  // ---------------------- inline and template functions ---------------------


  namespace MatrixIterators
  {
    inline const_iterator::Accessor::Accessor(const MatrixBase *matrix,
                                              const size_type   row,
                                              const size_type   index)
      : matrix(const_cast<MatrixBase *>(matrix))
      , a_row(row)
      , a_index(index)
    {
      visit_present_row();
    }



    inline const_iterator::Accessor::size_type
    const_iterator::Accessor::row() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_row;
    }


    inline const_iterator::Accessor::size_type
    const_iterator::Accessor::column() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*colnum_cache)[a_index];
    }


    inline const_iterator::Accessor::size_type
    const_iterator::Accessor::index() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_index;
    }


    inline PetscScalar
    const_iterator::Accessor::value() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*value_cache)[a_index];
    }


    inline const_iterator::const_iterator(const MatrixBase *matrix,
                                          const size_type   row,
                                          const size_type   index)
      : accessor(matrix, row, index)
    {}



    inline const_iterator &
    const_iterator::operator++()
    {
      Assert(accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

      ++accessor.a_index;

      // if at end of line: do one step, then cycle until we find a
      // row with a nonzero number of entries
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while ((accessor.a_row < accessor.matrix->m()) &&
                 (accessor.a_row < accessor.matrix->local_range().second) &&
                 (accessor.matrix->row_length(accessor.a_row) == 0))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }


    inline const_iterator
    const_iterator::operator++(int)
    {
      const const_iterator old_state = *this;
      ++(*this);
      return old_state;
    }


    inline const const_iterator::Accessor &const_iterator::operator*() const
    {
      return accessor;
    }


    inline const const_iterator::Accessor *const_iterator::operator->() const
    {
      return &accessor;
    }


    inline bool
    const_iterator::operator==(const const_iterator &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }


    inline bool
    const_iterator::operator!=(const const_iterator &other) const
    {
      return !(*this == other);
    }


    inline bool
    const_iterator::operator<(const const_iterator &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }

  } // namespace MatrixIterators



  // Inline the set() and add()
  // functions, since they will be
  // called frequently, and the
  // compiler can optimize away
  // some unnecessary loops when
  // the sizes are given at
  // compile time.
  inline void
  MatrixBase::set(const size_type i, const size_type j, const PetscScalar value)
  {
    AssertIsFinite(value);

    set(i, 1, &j, &value, false);
  }



  inline void
  MatrixBase::set(const std::vector<size_type> & indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
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
  MatrixBase::set(const std::vector<size_type> & row_indices,
                  const std::vector<size_type> & col_indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
  {
    Assert(row_indices.size() == values.m(),
           ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert(col_indices.size() == values.n(),
           ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i = 0; i < row_indices.size(); ++i)
      set(row_indices[i],
          col_indices.size(),
          col_indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  MatrixBase::set(const size_type                 row,
                  const std::vector<size_type> &  col_indices,
                  const std::vector<PetscScalar> &values,
                  const bool                      elide_zero_values)
  {
    Assert(col_indices.size() == values.size(),
           ExcDimensionMismatch(col_indices.size(), values.size()));

    set(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }



  inline void
  MatrixBase::set(const size_type    row,
                  const size_type    n_cols,
                  const size_type *  col_indices,
                  const PetscScalar *values,
                  const bool         elide_zero_values)
  {
    prepare_action(VectorOperation::insert);

    const PetscInt  petsc_i = row;
    PetscInt const *col_index_ptr;

    PetscScalar const *col_value_ptr;
    int                n_columns;

    // If we don't elide zeros, the pointers are already available...
    if (elide_zero_values == false)
      {
        col_index_ptr = reinterpret_cast<const PetscInt *>(col_indices);
        col_value_ptr = values;
        n_columns     = n_cols;
      }
    else
      {
        // Otherwise, extract nonzero values in each row and get the
        // respective index.
        if (column_indices.size() < n_cols)
          {
            column_indices.resize(n_cols);
            column_values.resize(n_cols);
          }

        n_columns = 0;
        for (size_type j = 0; j < n_cols; ++j)
          {
            const PetscScalar value = values[j];
            AssertIsFinite(value);
            if (value != PetscScalar())
              {
                column_indices[n_columns] = col_indices[j];
                column_values[n_columns]  = value;
                n_columns++;
              }
          }
        AssertIndexRange(n_columns, n_cols + 1);

        col_index_ptr = column_indices.data();
        col_value_ptr = column_values.data();
      }

    const PetscErrorCode ierr = MatSetValues(matrix,
                                             1,
                                             &petsc_i,
                                             n_columns,
                                             col_index_ptr,
                                             col_value_ptr,
                                             INSERT_VALUES);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  inline void
  MatrixBase::add(const size_type i, const size_type j, const PetscScalar value)
  {
    AssertIsFinite(value);

    if (value == PetscScalar())
      {
        // we have to check after using Insert/Add in any case to be
        // consistent with the MPI communication model, but we can save
        // some work if the addend is zero. However, these actions are done
        // in case we pass on to the other function.
        prepare_action(VectorOperation::add);

        return;
      }
    else
      add(i, 1, &j, &value, false);
  }



  inline void
  MatrixBase::add(const std::vector<size_type> & indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
  {
    Assert(indices.size() == values.m(),
           ExcDimensionMismatch(indices.size(), values.m()));
    Assert(values.m() == values.n(), ExcNotQuadratic());

    for (size_type i = 0; i < indices.size(); ++i)
      add(indices[i],
          indices.size(),
          indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  MatrixBase::add(const std::vector<size_type> & row_indices,
                  const std::vector<size_type> & col_indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
  {
    Assert(row_indices.size() == values.m(),
           ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert(col_indices.size() == values.n(),
           ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i = 0; i < row_indices.size(); ++i)
      add(row_indices[i],
          col_indices.size(),
          col_indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  MatrixBase::add(const size_type                 row,
                  const std::vector<size_type> &  col_indices,
                  const std::vector<PetscScalar> &values,
                  const bool                      elide_zero_values)
  {
    Assert(col_indices.size() == values.size(),
           ExcDimensionMismatch(col_indices.size(), values.size()));

    add(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }



  inline void
  MatrixBase::add(const size_type    row,
                  const size_type    n_cols,
                  const size_type *  col_indices,
                  const PetscScalar *values,
                  const bool         elide_zero_values,
                  const bool /*col_indices_are_sorted*/)
  {
    (void)elide_zero_values;

    prepare_action(VectorOperation::add);

    const PetscInt  petsc_i = row;
    PetscInt const *col_index_ptr;

    PetscScalar const *col_value_ptr;
    int                n_columns;

    // If we don't elide zeros, the pointers are already available...
    if (elide_zero_values == false)
      {
        col_index_ptr = reinterpret_cast<const PetscInt *>(col_indices);
        col_value_ptr = values;
        n_columns     = n_cols;
      }
    else
      {
        // Otherwise, extract nonzero values in each row and get the
        // respective index.
        if (column_indices.size() < n_cols)
          {
            column_indices.resize(n_cols);
            column_values.resize(n_cols);
          }

        n_columns = 0;
        for (size_type j = 0; j < n_cols; ++j)
          {
            const PetscScalar value = values[j];
            AssertIsFinite(value);
            if (value != PetscScalar())
              {
                column_indices[n_columns] = col_indices[j];
                column_values[n_columns]  = value;
                n_columns++;
              }
          }
        AssertIndexRange(n_columns, n_cols + 1);

        col_index_ptr = column_indices.data();
        col_value_ptr = column_values.data();
      }

    const PetscErrorCode ierr = MatSetValues(
      matrix, 1, &petsc_i, n_columns, col_index_ptr, col_value_ptr, ADD_VALUES);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  inline PetscScalar
  MatrixBase::operator()(const size_type i, const size_type j) const
  {
    return el(i, j);
  }



  inline MatrixBase::const_iterator
  MatrixBase::begin() const
  {
    Assert(
      (in_local_range(0) && in_local_range(m() - 1)),
      ExcMessage(
        "begin() and end() can only be called on a processor owning the entire matrix. If this is a distributed matrix, use begin(row) and end(row) instead."));

    // find the first non-empty row in order to make sure that the returned
    // iterator points to something useful
    size_type first_nonempty_row = 0;
    while ((first_nonempty_row < m()) && (row_length(first_nonempty_row) == 0))
      ++first_nonempty_row;

    return const_iterator(this, first_nonempty_row, 0);
  }


  inline MatrixBase::const_iterator
  MatrixBase::end() const
  {
    Assert(
      (in_local_range(0) && in_local_range(m() - 1)),
      ExcMessage(
        "begin() and end() can only be called on a processor owning the entire matrix. If this is a distributed matrix, use begin(row) and end(row) instead."));

    return const_iterator(this, m(), 0);
  }


  inline MatrixBase::const_iterator
  MatrixBase::begin(const size_type r) const
  {
    Assert(in_local_range(r),
           ExcIndexRange(r, local_range().first, local_range().second));

    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end(r);
  }


  inline MatrixBase::const_iterator
  MatrixBase::end(const size_type r) const
  {
    Assert(in_local_range(r),
           ExcIndexRange(r, local_range().first, local_range().second));

    // place the iterator on the first entry past this line, or at the
    // end of the matrix
    //
    // in the parallel case, we need to put it on the first entry of
    // the first row after the locally owned range. this of course
    // doesn't exist, but we can nevertheless create such an
    // iterator. we need to check whether 'i' is past the locally
    // owned range of rows first, before we ask for the length of the
    // row since the latter query leads to an exception in case we ask
    // for a row that is not locally owned
    for (size_type i = r + 1; i < m(); ++i)
      if (i == local_range().second || (row_length(i) > 0))
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    // we don't allow calling end() directly for distributed matrices so we need
    // to copy the code without the assertion.
    return {this, m(), 0};
  }



  inline bool
  MatrixBase::in_local_range(const size_type index) const
  {
    PetscInt begin, end;

    const PetscErrorCode ierr =
      MatGetOwnershipRange(static_cast<const Mat &>(matrix), &begin, &end);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline void
  MatrixBase::prepare_action(const VectorOperation::values new_action)
  {
    if (last_action == VectorOperation::unknown)
      last_action = new_action;

    Assert(last_action == new_action, ExcWrongMode(last_action, new_action));
  }



  inline void
  MatrixBase::assert_is_compressed()
  {
    // compress() sets the last action to none, which allows us to check if
    // there are pending add/insert operations:
    AssertThrow(last_action == VectorOperation::unknown,
                ExcMessage("Error: missing compress() call."));
  }



  inline void
  MatrixBase::prepare_add()
  {
    prepare_action(VectorOperation::add);
  }



  inline void
  MatrixBase::prepare_set()
  {
    prepare_action(VectorOperation::insert);
  }

#    endif // DOXYGEN
} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_PETSC

#endif
/*---------------------------- petsc_matrix_base.h --------------------------*/
