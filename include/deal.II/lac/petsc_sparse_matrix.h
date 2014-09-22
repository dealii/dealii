// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__petsc_sparse_matrix_h
#define __deal2__petsc_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_matrix_base.h>
#  include <vector>

DEAL_II_NAMESPACE_OPEN
// forward declaration
template <typename MatrixType> class BlockMatrixBase;


namespace PETScWrappers
{
  /**
   * Implementation of a sequential sparse matrix class based on PETSC. All the
   * functionality is actually in the base class, except for the calls to
   * generate a sequential sparse matrix. This is possible since PETSc only works
   * on an abstract matrix type and internally distributes to functions that do
   * the actual work depending on the actual matrix type (much like using
   * virtual functions). Only the functions creating a matrix of specific type
   * differ, and are implemented in this particular class.
   *
   * @ingroup PETScWrappers
   * @ingroup Matrix1
   * @author Wolfgang Bangerth, 2004
   */
  class SparseMatrix : public MatrixBase
  {
  public:

    /**
     * A structure that describes some of
     * the traits of this class in terms of
     * its run-time behavior. Some other
     * classes (such as the block matrix
     * classes) that take one or other of
     * the matrix classes as its template
     * parameters can tune their behavior
     * based on the variables in this
     * class.
     */
    struct Traits
    {
      /**
       * It is safe to elide additions of
       * zeros to individual elements of
       * this matrix.
       */
      static const bool zero_addition_can_be_elided = true;
    };

    /**
     * Default constructor. Create an empty
     * matrix.
     */
    SparseMatrix ();

    /**
     * Create a sparse matrix of dimensions
     * @p m times @p n, with an
     * initial guess of
     * @p n_nonzero_per_row nonzero
     * elements per row. PETSc is able to
     * cope with the situation that more
     * than this number of elements is
     * later allocated for a row, but this
     * involves copying data, and is thus
     * expensive.
     *
     * The @p is_symmetric flag determines
     * whether we should tell PETSc that
     * the matrix is going to be symmetric
     * (as indicated by the call
     * <tt>MatSetOption(mat,
     * MAT_SYMMETRIC)</tt>. Note that the
     * PETSc documentation states that one
     * cannot form an ILU decomposition of
     * a matrix for which this flag has
     * been set to @p true, only an
     * ICC. The default value of this flag
     * is @p false.
     */
    SparseMatrix (const size_type m,
                  const size_type n,
                  const size_type n_nonzero_per_row,
                  const bool      is_symmetric = false);

    /**
     * Initialize a rectangular matrix with
     * @p m rows and @p n
     * columns.  The maximal number of
     * nonzero entries for each row
     * separately is given by the
     * @p row_lengths array.
     *
     * Just as for the other constructors:
     * PETSc is able to cope with the
     * situation that more than this number
     * of elements is later allocated for a
     * row, but this involves copying data,
     * and is thus expensive.
     *
     * The @p is_symmetric flag determines
     * whether we should tell PETSc that
     * the matrix is going to be symmetric
     * (as indicated by the call
     * <tt>MatSetOption(mat,
     * MAT_SYMMETRIC)</tt>. Note that the
     * PETSc documentation states that one
     * cannot form an ILU decomposition of
     * a matrix for which this flag has
     * been set to @p true, only an
     * ICC. The default value of this flag
     * is @p false.
     */
    SparseMatrix (const size_type               m,
                  const size_type               n,
                  const std::vector<size_type> &row_lengths,
                  const bool                    is_symmetric = false);

    /**
     * Initialize a sparse matrix using the
     * given sparsity pattern.
     *
     * Note that PETSc can be very slow
     * if you do not provide it with a
     * good estimate of the lengths of
     * rows. Using the present function
     * is a very efficient way to do
     * this, as it uses the exact number
     * of nonzero entries for each row of
     * the matrix by using the given
     * sparsity pattern argument. If the
     * @p preset_nonzero_locations flag
     * is @p true, this function in
     * addition not only sets the correct
     * row sizes up front, but also
     * pre-allocated the correct nonzero
     * entries in the matrix.
     *
     * PETsc allows to later add
     * additional nonzero entries to a
     * matrix, by simply writing to these
     * elements. However, this will then
     * lead to additional memory
     * allocations which are very
     * inefficient and will greatly slow
     * down your program. It is therefore
     * significantly more efficient to
     * get memory allocation right from
     * the start.
     */
    template <typename SparsityType>
    explicit SparseMatrix (const SparsityType &sparsity_pattern,
                           const bool          preset_nonzero_locations = true);

    /**
     * This operator assigns a scalar to
     * a matrix. Since this does usually
     * not make much sense (should we set
     * all matrix entries to this value?
     * Only the nonzero entries of the
     * sparsity pattern?), this operation
     * is only allowed if the actual
     * value to be assigned is zero. This
     * operator only exists to allow for
     * the obvious notation
     * <tt>matrix=0</tt>, which sets all
     * elements of the matrix to zero,
     * but keep the sparsity pattern
     * previously used.
     */
    SparseMatrix &operator = (const double d);

    /**
     * Throw away the present matrix and
     * generate one that has the same
     * properties as if it were created by
     * the constructor of this class with
     * the same argument list as the
     * present function.
     */
    void reinit (const size_type m,
                 const size_type n,
                 const size_type n_nonzero_per_row,
                 const bool      is_symmetric = false);

    /**
     * Throw away the present matrix and
     * generate one that has the same
     * properties as if it were created by
     * the constructor of this class with
     * the same argument list as the
     * present function.
     */
    void reinit (const size_type               m,
                 const size_type               n,
                 const std::vector<size_type> &row_lengths,
                 const bool                    is_symmetric = false);

    /**
     * Initialize a sparse matrix using the
     * given sparsity pattern.
     *
     * Note that PETSc can be very slow
     * if you do not provide it with a
     * good estimate of the lengths of
     * rows. Using the present function
     * is a very efficient way to do
     * this, as it uses the exact number
     * of nonzero entries for each row of
     * the matrix by using the given
     * sparsity pattern argument. If the
     * @p preset_nonzero_locations flag
     * is @p true, this function in
     * addition not only sets the correct
     * row sizes up front, but also
     * pre-allocated the correct nonzero
     * entries in the matrix.
     *
     * PETsc allows to later add
     * additional nonzero entries to a
     * matrix, by simply writing to these
     * elements. However, this will then
     * lead to additional memory
     * allocations which are very
     * inefficient and will greatly slow
     * down your program. It is therefore
     * significantly more efficient to
     * get memory allocation right from
     * the start.
     *
     * Despite the fact that it would
     * seem to be an obvious win, setting
     * the @p preset_nonzero_locations
     * flag to @p true doesn't seem to
     * accelerate program. Rather on the
     * contrary, it seems to be able to
     * slow down entire programs
     * somewhat. This is suprising, since
     * we can use efficient function
     * calls into PETSc that allow to
     * create multiple entries at once;
     * nevertheless, given the fact that
     * it is inefficient, the respective
     * flag has a default value equal to
     * @p false.
     */
    template <typename SparsityType>
    void reinit (const SparsityType &sparsity_pattern,
                 const bool          preset_nonzero_locations = true);

    /**
     * Return a reference to the MPI
     * communicator object in use with this
     * matrix. Since this is a sequential
     * matrix, it returns the MPI_COMM_SELF
     * communicator.
     */
    virtual const MPI_Comm &get_mpi_communicator () const;

    /**
     * Return the square of the norm
     * of the vector $v$ with respect
     * to the norm induced by this
     * matrix,
     * i.e. $\left(v,Mv\right)$. This
     * is useful, e.g. in the finite
     * element context, where the
     * $L_2$ norm of a function
     * equals the matrix norm with
     * respect to the mass matrix of
     * the vector representing the
     * nodal values of the finite
     * element function.
     *
     * Obviously, the matrix needs to
     * be quadratic for this operation.
     *
     * The implementation of this function
     * is not as efficient as the one in
     * the @p MatrixBase class used in
     * deal.II (i.e. the original one, not
     * the PETSc wrapper class) since PETSc
     * doesn't support this operation and
     * needs a temporary vector.
     */
    PetscScalar matrix_norm_square (const VectorBase &v) const;

    /**
     * Compute the matrix scalar
     * product $\left(u,Mv\right)$.
     *
     * The implementation of this function
     * is not as efficient as the one in
     * the @p MatrixBase class used in
     * deal.II (i.e. the original one, not
     * the PETSc wrapper class) since PETSc
     * doesn't support this operation and
     * needs a temporary vector.
     */
    PetscScalar matrix_scalar_product (const VectorBase &u,
                                       const VectorBase &v) const;

  private:

    /**
     * Purposefully not implemented
     */
    SparseMatrix(const SparseMatrix &);
    /**
     * Purposefully not implemented
     */
    SparseMatrix &operator= (const SparseMatrix &);

    /**
     * Do the actual work for the
     * respective reinit() function and the
     * matching constructor, i.e. create a
     * matrix. Getting rid of the previous
     * matrix is left to the caller.
     */
    void do_reinit (const size_type m,
                    const size_type n,
                    const size_type n_nonzero_per_row,
                    const bool      is_symmetric = false);

    /**
     * Same as previous function.
     */
    void do_reinit (const size_type               m,
                    const size_type               n,
                    const std::vector<size_type> &row_lengths,
                    const bool                    is_symmetric = false);

    /**
     * Same as previous function.
     */
    template <typename SparsityType>
    void do_reinit (const SparsityType &sparsity_pattern,
                    const bool          preset_nonzero_locations);

    /**
     *  To allow calling protected
     *  prepare_add() and
     *  prepare_set().
     */
    friend class BlockMatrixBase<SparseMatrix>;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

/*----------------------------   petsc_sparse_matrix.h     ---------------------------*/

#endif
/*----------------------------   petsc_sparse_matrix.h     ---------------------------*/
