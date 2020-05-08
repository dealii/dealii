// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_sparse_direct_h
#define dealii_sparse_direct_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_UMFPACK
#  include <umfpack.h>
#endif

DEAL_II_NAMESPACE_OPEN

namespace types
{
  /**
   * Index type for UMFPACK. SuiteSparse_long has to be used here for the
   * Windows 64 build.
   */
#ifdef SuiteSparse_long
  using suitesparse_index = SuiteSparse_long;
#else
  using suitesparse_index = long int;
#endif
} // namespace types

/**
 * This class provides an interface to the sparse direct solver UMFPACK, which
 * is part of the SuiteSparse library (see <a
 * href="http://faculty.cse.tamu.edu/davis/suitesparse.html">this link</a>).
 * UMFPACK is a set of routines for solving non-symmetric sparse linear
 * systems, Ax=b, using the Unsymmetric-pattern MultiFrontal method and direct
 * sparse LU factorization. Matrices may have symmetric or unsymmetric
 * sparsity patterns, and may have unsymmetric entries. The use of this class
 * is explained in the step-22 and step-29
 * tutorial programs.
 *
 * This matrix class implements the usual interface of preconditioners, that
 * is a function initialize(const SparseMatrix<double>&matrix,const
 * AdditionalData) for initializing and the whole set of vmult() functions
 * common to all matrices. Implemented here are only vmult() and vmult_add(),
 * which perform multiplication with the inverse matrix. Furthermore, this
 * class provides an older interface, consisting of the functions factorize()
 * and solve(). Both interfaces are interchangeable.
 *
 * @note This class exists if the <a
 * href="http://faculty.cse.tamu.edu/davis/suitesparse.html">UMFPACK</a>
 * interface was not explicitly disabled during configuration.
 *
 * @note UMFPACK has its own license, independent of that of deal.II. If you
 * want to use the UMFPACK you have to accept that license. It is linked to
 * from the deal.II ReadMe file. UMFPACK is included courtesy of its author,
 * <a href="http://faculty.cse.tamu.edu/davis/welcome.html">Timothy
 * A. Davis</a>.
 *
 *
 * <h4>Instantiations</h4>
 *
 * There are instantiations of this class for SparseMatrix<double>,
 * SparseMatrix<float>, SparseMatrixEZ<float>, SparseMatrixEZ<double>,
 * BlockSparseMatrix<double>, and BlockSparseMatrix<float>.
 *
 * @ingroup Solvers Preconditioners
 *
 * @author Wolfgang Bangerth, 2004; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
class SparseDirectUMFPACK : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Dummy class needed for the usual initialization interface of
   * preconditioners.
   */
  class AdditionalData
  {};


  /**
   * Constructor. See the documentation of this class for the meaning of the
   * parameters to this function.
   */
  SparseDirectUMFPACK();

  /**
   * Destructor.
   */
  ~SparseDirectUMFPACK() override;

  /**
   * @name Setting up a sparse factorization
   */
  /**
   * @{
   */

  /**
   * This function does nothing. It is only here to provide a interface
   * consistent with other sparse direct solvers.
   */
  void
  initialize(const SparsityPattern &sparsity_pattern);

  /**
   * Factorize the matrix. This function may be called multiple times for
   * different matrices, after the object of this class has been initialized
   * for a certain sparsity pattern. You may therefore save some computing
   * time if you want to invert several matrices with the same sparsity
   * pattern. However, note that the bulk of the computing time is actually
   * spent in the factorization, so this functionality may not always be of
   * large benefit.
   *
   * In contrast to the other direct solver classes, the initialization method
   * does nothing. Therefore initialize is not automatically called by this
   * method, when the initialization step has not been performed yet.
   *
   * This function copies the contents of the matrix into its own storage; the
   * matrix can therefore be deleted after this operation, even if subsequent
   * solves are required.
   */
  template <class Matrix>
  void
  factorize(const Matrix &matrix);

  /**
   * Initialize memory and call SparseDirectUMFPACK::factorize.
   */
  template <class Matrix>
  void
  initialize(const Matrix &       matrix,
             const AdditionalData additional_data = AdditionalData());

  /**
   * @}
   */

  /**
   * @name Functions that represent the inverse of a matrix
   */
  /**
   * @{
   */

  /**
   * Preconditioner interface function. Usually, given the source vector, this
   * method returns an approximate solution of <i>Ax = b</i>. As this class
   * provides a wrapper to a direct solver, here it is actually the exact
   * solution (exact within the range of numerical accuracy of course).
   *
   * In other words, this function actually multiplies with the exact inverse
   * of the matrix, $A^{-1}$.
   */
  void
  vmult(Vector<double> &dst, const Vector<double> &src) const;

  /**
   * Same as before, but for block vectors.
   */
  void
  vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;

  /**
   * Same as before, but uses the transpose of the matrix, i.e. this function
   * multiplies with $A^{-T}$.
   */
  void
  Tvmult(Vector<double> &dst, const Vector<double> &src) const;

  /**
   * Same as before, but for block vectors
   */
  void
  Tvmult(BlockVector<double> &dst, const BlockVector<double> &src) const;

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   */
  size_type
  n() const;

  /**
   * @}
   */

  /**
   * @name Functions that solve linear systems
   */
  /**
   * @{
   */

  /**
   * Solve for a certain right hand side vector. This function may be called
   * multiple times for different right hand side vectors after the matrix has
   * been factorized. This yields substantial savings in computing time, since
   * the actual solution is fast, compared to the factorization of the matrix.
   *
   * The solution will be returned in place of the right hand side vector.
   *
   * @param[in,out] rhs_and_solution A vector that contains the right hand side
   *   $b$ of a linear system $Ax=b$ upon calling this function, and that
   *   contains the solution $x$ of the linear system after calling this
   *   function.
   * @param[in] transpose If set to true, this function solves the linear
   *   $A^T x = b$ instead of $Ax=b$.
   *
   * @pre You need to call factorize() before this function can be called.
   */
  void
  solve(Vector<double> &rhs_and_solution, const bool transpose = false) const;

  /**
   * Like the previous function, but for a complex-valued right hand side
   * and solution vector.
   *
   * If the matrix that was previously factorized had complex-valued entries,
   * then the `rhs_and_solution` vector will, upon return from this function,
   * simply contain the solution of the linear system $Ax=b$. If the matrix
   * was real-valued, then this is also true, but the solution will simply
   * be computed by applying the factorized $A^{-1}$ to both the
   * real and imaginary parts of the right hand side vector.
   */
  void
  solve(Vector<std::complex<double>> &rhs_and_solution,
        const bool                    transpose = false) const;

  /**
   * Same as before, but for block vectors.
   */
  void
  solve(BlockVector<double> &rhs_and_solution,
        const bool           transpose = false) const;

  /**
   * Same as before, but for complex-valued block vectors.
   */
  void
  solve(BlockVector<std::complex<double>> &rhs_and_solution,
        const bool                         transpose = false) const;

  /**
   * Call the two functions factorize() and solve() in that order, i.e.
   * perform the whole solution process for the given right hand side vector.
   *
   * The solution will be returned in place of the right hand side vector.
   */
  template <class Matrix>
  void
  solve(const Matrix &  matrix,
        Vector<double> &rhs_and_solution,
        const bool      transpose = false);

  /**
   * Same as before, but for complex-valued solution vectors.
   */
  template <class Matrix>
  void
  solve(const Matrix &                matrix,
        Vector<std::complex<double>> &rhs_and_solution,
        const bool                    transpose = false);

  /**
   * Same as before, but for block vectors.
   */
  template <class Matrix>
  void
  solve(const Matrix &       matrix,
        BlockVector<double> &rhs_and_solution,
        const bool           transpose = false);

  /**
   * Same as before, but for complex-valued block vectors.
   */
  template <class Matrix>
  void
  solve(const Matrix &                     matrix,
        BlockVector<std::complex<double>> &rhs_and_solution,
        const bool                         transpose = false);

  /**
   * @}
   */

  /**
   * One of the UMFPack routines threw an error. The error code is included in
   * the output and can be looked up in the UMFPack user manual. The name of
   * the routine is included for reference.
   */
  DeclException2(
    ExcUMFPACKError,
    std::string,
    int,
    << "UMFPACK routine " << arg1 << " returned error status " << arg2 << "."
    << "\n\n"
    << ("A complete list of error codes can be found in the file "
        "<bundled/umfpack/UMFPACK/Include/umfpack.h>."
        "\n\n"
        "That said, the two most common errors that can happen are "
        "that your matrix cannot be factorized because it is "
        "rank deficient, and that UMFPACK runs out of memory "
        "because your problem is too large."
        "\n\n"
        "The first of these cases most often happens if you "
        "forget terms in your bilinear form necessary to ensure "
        "that the matrix has full rank, or if your equation has a "
        "spatially variable coefficient (or nonlinearity) that is "
        "supposed to be strictly positive but, for whatever "
        "reasons, is negative or zero. In either case, you probably "
        "want to check your assembly procedure. Similarly, a "
        "matrix can be rank deficient if you forgot to apply the "
        "appropriate boundary conditions. For example, the "
        "Laplace equation for a problem where only Neumann boundary "
        "conditions are posed (or where you forget to apply Dirichlet "
        "boundary conditions) has exactly one eigenvalue equal to zero "
        "and its rank is therefore deficient by one. Finally, the matrix "
        "may be rank deficient because you are using a quadrature "
        "formula with too few quadrature points."
        "\n\n"
        "The other common situation is that you run out of memory. "
        "On a typical laptop or desktop, it should easily be possible "
        "to solve problems with 100,000 unknowns in 2d. If you are "
        "solving problems with many more unknowns than that, in "
        "particular if you are in 3d, then you may be running out "
        "of memory and you will need to consider iterative "
        "solvers instead of the direct solver employed by "
        "UMFPACK."));

private:
  /**
   * The dimension of the range space, i.e., the number of rows of the matrix.
   */
  size_type n_rows;

  /**
   * The dimension of the domain space, i.e., the number of columns of the
   * matrix.
   */
  size_type n_cols;

  /**
   * The UMFPACK routines allocate objects in which they store information
   * about symbolic and numeric values of the decomposition. The actual data
   * type of these objects is opaque, and only passed around as void pointers.
   */
  void *symbolic_decomposition;
  void *numeric_decomposition;

  /**
   * Free all memory that hasn't been freed yet.
   */
  void
  clear();

  /**
   * Make sure that the arrays Ai and Ap are sorted in each row. UMFPACK wants
   * it this way. We need to have three versions of this function, one for the
   * usual SparseMatrix, one for the SparseMatrixEZ, and one for the
   * BlockSparseMatrix classes
   */
  template <typename number>
  void
  sort_arrays(const SparseMatrixEZ<number> &);

  template <typename number>
  void
  sort_arrays(const SparseMatrix<number> &);

  template <typename number>
  void
  sort_arrays(const BlockSparseMatrix<number> &);

  /**
   * The arrays in which we store the data for the solver. These are documented
   * in the descriptions of the umfpack_*_symbolic() and umfpack_*_numeric()
   * functions, but in short:
   * - `Ap` is the array saying which row starts where in `Ai`
   * - `Ai` is the array that stores the column indices of nonzero entries
   * - `Ax` is the array that stores the values of nonzero entries; if the
   *   matrix is complex-valued, then it stores the real parts
   * - `Az` is the array that stores the imaginary parts of nonzero entries,
   *   and is used only if the matrix is complex-valued.
   */
  std::vector<types::suitesparse_index> Ap;
  std::vector<types::suitesparse_index> Ai;
  std::vector<double>                   Ax;
  std::vector<double>                   Az;

  /**
   * Control and work arrays for the solver routines.
   */
  std::vector<double> control;
};

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_sparse_direct_h
