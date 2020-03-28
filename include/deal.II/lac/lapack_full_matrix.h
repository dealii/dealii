// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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

#ifndef dealii_lapack_full_matrix_h
#define dealii_lapack_full_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector_memory.h>

#include <complex>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class BlockVector;
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;
#endif

/**
 * A variant of FullMatrix using LAPACK functions wherever possible. In order
 * to do this, the matrix is stored in transposed order. The element access
 * functions hide this fact by reverting the transposition.
 *
 * @note In order to perform LAPACK functions, the class contains a lot of
 * auxiliary data in the private section. The names of these data vectors are
 * usually the names chosen for the arguments in the LAPACK documentation.
 *
 * @ingroup Matrix1
 * @author Guido Kanschat, 2005, Denis Davydov, 2017, 2018
 */
template <typename number>
class LAPACKFullMatrix : public TransposeTable<number>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = std::make_unsigned<types::blas_int>::type;

  /**
   * Constructor. Initialize the matrix as a square matrix with dimension
   * @p size.
   *
   * In order to avoid the implicit conversion of integers and other types to
   * a matrix, this constructor is declared <tt>explicit</tt>.
   *
   * By default, no memory is allocated.
   */
  explicit LAPACKFullMatrix(const size_type size = 0);


  /**
   * Constructor. Initialize the matrix as a rectangular matrix $\rm{rows}
   * \times \rm{cols}$.
   */
  LAPACKFullMatrix(const size_type rows, const size_type cols);


  /**
   * Copy constructor. This constructor does a deep copy of the matrix.
   * Therefore, it poses a possible efficiency problem, if for example,
   * function arguments are passed by value rather than by reference.
   * Unfortunately, we can't mark this copy constructor <tt>explicit</tt>,
   * since that prevents the use of this class in containers, such as
   * <code>std::vector</code>. The responsibility to check performance of
   * programs must therefore remain with the user of this class.
   */
  LAPACKFullMatrix(const LAPACKFullMatrix &);

  /**
   * Assignment operator.
   */
  LAPACKFullMatrix<number> &
  operator=(const LAPACKFullMatrix<number> &);

  /**
   * Assignment operator from a regular FullMatrix.
   *
   * @note Since LAPACK expects matrices in transposed order, this
   * transposition is included here.
   */
  template <typename number2>
  LAPACKFullMatrix<number> &
  operator=(const FullMatrix<number2> &);

  /**
   * Assignment operator from a regular SparseMatrix.
   *
   * @note Since LAPACK expects matrices in transposed order, this
   * transposition is included here.
   */
  template <typename number2>
  LAPACKFullMatrix<number> &
  operator=(const SparseMatrix<number2> &);

  /**
   * This operator assigns a scalar to a matrix. To avoid confusion with
   * constructors, zero (when cast to the @p number type) is the only
   * value allowed for @p d.
   */
  LAPACKFullMatrix<number> &
  operator=(const number d);

  /**
   * This operator multiplies all entries by a fixed @p factor.
   */
  LAPACKFullMatrix<number> &
  operator*=(const number factor);

  /**
   * This operator divides all entries by a fixed @p factor.
   */
  LAPACKFullMatrix<number> &
  operator/=(const number factor);

  /**
   * Set a particular entry of the matrix to a @p value.
   * Thus, calling <code>A.set(1,2,3.141);</code> is entirely equivalent to the
   * operation <code>A(1,2) = 3.141;</code>. This function exists for
   * compatibility with the various sparse matrix objects.
   *
   * @param i The row index of the element to be set.
   * @param j The column index of the element to be set.
   * @param value The value to be written into the element.
   */
  void
  set(const size_type i, const size_type j, const number value);

  /**
   * Simple addition of a scaled matrix, i.e. $\mathbf A \mathrel{+}= a \,
   * \mathbf B$.
   */
  void
  add(const number a, const LAPACKFullMatrix<number> &B);

  /**
   * Perform a rank-1 update of a symmetric matrix
   * $ \mathbf A \leftarrow \mathbf A + a \, \mathbf v \mathbf v^T $.
   *
   * This function also works for Cholesky factorization.
   * In that case, updating ($a>0$) is
   * performed via Givens rotations, whereas downdating ($a<0$) via hyperbolic
   * rotations. Note that the latter case might lead to a negative definite
   * matrix in which case the error will be thrown (because Cholesky
   * factorizations are only valid for symmetric and positive definite
   * matrices).
   */
  void
  rank1_update(const number a, const Vector<number> &v);

  /**
   * Apply <a href="https://en.wikipedia.org/wiki/Givens_rotation">Givens
   * rotation</a>
   * @p csr (a triplet of cosine, sine and radius, see
   * Utilities::LinearAlgebra::givens_rotation() for the definition of the
   * rotation matrix $\mathbf G$)
   * to this matrix in the plane spanned by the @p i'th and @p k'th unit vectors.
   * If @p left is <code>true</code>, the rotation is applied from left
   * $\mathbf A \leftarrow \mathbf G \mathbf A$
   * and only rows @p i and @p k are affected.
   * Otherwise, transpose of the rotation matrix is applied from right
   * $\mathbf A \leftarrow \mathbf A \mathbf G^T$
   * and only columns @p i and @p k are affected.
   */
  void
  apply_givens_rotation(const std::array<number, 3> &csr,
                        const size_type              i,
                        const size_type              k,
                        const bool                   left = true);

  /**
   * Assignment from different matrix classes, performing the usual conversion
   * to the transposed format expected by LAPACK. This assignment operator
   * uses iterators of the typename MatrixType. Therefore, sparse matrices are
   * possible sources.
   */
  template <typename MatrixType>
  void
  copy_from(const MatrixType &);

  /**
   * Regenerate the current matrix by one that has the same properties as if
   * it were created by the constructor of this class with the same argument
   * list as this present function.
   */
  void
  reinit(const size_type size);

  /**
   * Same as above but will preserve the values of matrix upon resizing.
   * The original values
   * of the matrix are kept on increasing the size
   * \f[
   * \mathbf A \rightarrow
   * \left(
   * \begin{array}{cc}
   * \mathbf A & \mathbf 0 \\
   * \mathbf 0 & \mathbf 0
   * \end{array}
   * \right)
   * \f]
   * Whereas if the new size is smaller, the matrix will contain the upper left
   * block of the original one \f[ \left( \begin{array}{cc}
   * \mathbf A_{11} & \mathbf A_{12} \\
   * \mathbf A_{21} & \mathbf A_{22}
   * \end{array}
   * \right)
   * \rightarrow \mathbf A_{11}
   * \f]
   */
  void
  grow_or_shrink(const size_type size);

  /**
   * Remove row @p row and column @p col from the matrix.
   * \f[
   * \left(
   * \begin{array}{ccc}
   * \mathbf A_{11}   & \mathbf a_{12} & \mathbf A_{13}   \\
   * \mathbf a_{21}^T & a_{22}         & \mathbf a_{23}^T \\
   * \mathbf A_{31}   & \mathbf a_{32} & \mathbf A_{33}
   * \end{array}
   * \right)
   * \rightarrow
   * \left(
   * \begin{array}{cc}
   * \mathbf A_{11} & \mathbf A_{13} \\
   * \mathbf A_{31} & \mathbf A_{33}
   * \end{array}
   * \right)
   * \f]
   */
  void
  remove_row_and_column(const size_type row, const size_type col);

  /**
   * Regenerate the current matrix by one that has the same properties as if
   * it were created by the constructor of this class with the same argument
   * list as this present function.
   */
  void
  reinit(const size_type rows, const size_type cols);

  /**
   * Assign @p property to this matrix.
   */
  void
  set_property(const LAPACKSupport::Property property);

  /**
   * Return the dimension of the codomain (or range) space.
   *
   * @note The matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space.
   *
   * @note The matrix is of dimension $m \times n$.
   */
  size_type
  n() const;

  /**
   * Fill rectangular block.
   *
   * A rectangular block of the matrix <tt>src</tt> is copied into
   * <tt>this</tt>. The upper left corner of the block being copied is
   * <tt>(src_offset_i,src_offset_j)</tt>.  The upper left corner of the
   * copied block is <tt>(dst_offset_i,dst_offset_j)</tt>.  The size of the
   * rectangular block being copied is the maximum size possible, determined
   * either by the size of <tt>this</tt> or <tt>src</tt>.
   *
   * The final two arguments allow to enter a multiple of the source or its
   * transpose.
   */
  template <typename MatrixType>
  void
  fill(const MatrixType &src,
       const size_type   dst_offset_i = 0,
       const size_type   dst_offset_j = 0,
       const size_type   src_offset_i = 0,
       const size_type   src_offset_j = 0,
       const number      factor       = 1.,
       const bool        transpose    = false);


  /**
   * Matrix-vector-multiplication.
   *
   * Depending on previous transformations recorded in #state, the result of
   * this function is one of
   * <ul>
   * <li> If #state is LAPACKSupport::matrix or LAPACKSupport::inverse_matrix,
   * this will be a regular matrix vector product using LAPACK gemv().
   * <li> If #state is LAPACKSupport::svd or LAPACKSupport::inverse_svd, this
   * function first multiplies with the right transformation matrix, then with
   * the diagonal matrix of singular values or their reciprocal values, and
   * finally with the left transformation matrix.
   * </ul>
   *
   * The optional parameter @p adding determines, whether the result is
   * stored in the vector
   * $\mathbf w = \mathbf A \cdot \mathbf v$
   * or added to it
   * $\mathbf w \mathrel{+}= \mathbf A \cdot \mathbf v$.
   *
   * @note Source and destination must not be the same vector.
   *
   * @note The template with @p number2 only exists for compile-time
   * compatibility with FullMatrix. Only the case @p number2 = @p number is
   * implemented due to limitations in the underlying LAPACK interface. All
   * other variants throw an error upon invocation.
   */
  template <typename number2>
  void
  vmult(Vector<number2> &      w,
        const Vector<number2> &v,
        const bool             adding = false) const;

  /**
   * Specialization of above function for compatible Vector::value_type.
   */
  void
  vmult(Vector<number> &      w,
        const Vector<number> &v,
        const bool            adding = false) const;

  /**
   * Adding Matrix-vector-multiplication $\mathbf w \mathrel{+}= \mathbf A \cdot
   * \mathbf v$.
   *
   * See the documentation of vmult() for details on the implementation.
   */
  template <typename number2>
  void
  vmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * Specialization of above function for compatible Vector::value_type.
   */
  void
  vmult_add(Vector<number> &w, const Vector<number> &v) const;

  /**
   * Transpose matrix-vector-multiplication.
   *
   * The optional parameter @p adding determines, whether the result is
   * stored in the vector
   * $\mathbf w = \mathbf A^T \cdot \mathbf v$
   * or added to it
   * $\mathbf w \mathrel{+}= \mathbf A^T \cdot \mathbf v$.
   *
   * See the documentation of vmult() for details on the implementation.
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &      w,
         const Vector<number2> &v,
         const bool             adding = false) const;

  /**
   * Specialization of above function for compatible Vector::value_type.
   */
  void
  Tvmult(Vector<number> &      w,
         const Vector<number> &v,
         const bool            adding = false) const;

  /**
   * Adding transpose matrix-vector-multiplication $\mathbf w \mathrel{+}=
   * \mathbf A^T \cdot \mathbf v$.
   *
   * See the documentation of vmult() for details on the implementation.
   */
  template <typename number2>
  void
  Tvmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * Specialization of above function for compatible Vector::value_type.
   */
  void
  Tvmult_add(Vector<number> &w, const Vector<number> &v) const;


  /**
   * Matrix-matrix-multiplication.
   *
   * The optional parameter @p adding determines, whether the result is
   * stored in the matrix
   * $\mathbf C            = \mathbf A \cdot \mathbf B$
   * or added to it
   * $\mathbf C \mathrel{+}= \mathbf A \cdot \mathbf B$.
   *
   * @note It is assumed that @p A and @p B have compatible sizes and that
   * @p C already has the right size.
   *
   * @p This function uses the BLAS function Xgemm.
   */
  void
  mmult(LAPACKFullMatrix<number> &      C,
        const LAPACKFullMatrix<number> &B,
        const bool                      adding = false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void
  mmult(FullMatrix<number> &            C,
        const LAPACKFullMatrix<number> &B,
        const bool                      adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt>.
   *
   * The optional parameter @p adding determines, whether the result is
   * stored in the matrix
   * $\mathbf C = \mathbf A^T \cdot \mathbf B$
   * or added to it
   * $\mathbf C \mathrel{+}= \mathbf A^T \cdot \mathbf B$.
   *
   * @note It is assumed that @p A and @p B have compatible sizes and that
   * @p C already has the right size.
   *
   * @note This function uses the BLAS function Xgemm.
   */
  void
  Tmmult(LAPACKFullMatrix<number> &      C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void
  Tmmult(FullMatrix<number> &            C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt> and a
   * diagonal vector @p V.
   *
   * If the <code>adding=false</code> then the result is stored in the matrix
   * $\mathbf C = \mathbf A^T \cdot \rm{diag}(\mathbf V) \cdot \mathbf B$
   * otherwise it is added $\mathbf C \mathrel{+}= \mathbf A^T \cdot
   * \rm{diag}(\mathbf V) \cdot \mathbf B$.
   *
   * @note It is assumed that @p A, @p B and @p V have compatible sizes and that
   * @p C already has the right size.
   *
   * @note This function is not provided by LAPACK. The function first forms $\rm{diag}(\mathbf V) \cdot \mathbf B$ product and
   * then uses Xgemm function.
   */
  void
  Tmmult(LAPACKFullMatrix<number> &      C,
         const LAPACKFullMatrix<number> &B,
         const Vector<number> &          V,
         const bool                      adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of @p B.
   *
   * The optional parameter @p adding determines, whether the result is
   * stored in the matrix
   * $\mathbf C = \mathbf A \cdot \mathbf B^T$
   * or added to it
   * $\mathbf C \mathrel{+}= \mathbf A \cdot \mathbf B^T$.
   *
   * @note It is assumed that @p A and @p B have compatible sizes and that
   * @p C already has the right size.
   *
   * @note This function uses the BLAS function Xgemm.
   */
  void
  mTmult(LAPACKFullMatrix<number> &      C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void
  mTmult(FullMatrix<number> &            C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt> and
   * @p B.
   *
   * The optional parameter @p adding determines, whether the result is
   * stored in the matrix
   * $\mathbf C = \mathbf A^T \cdot \mathbf B^T$
   * or added to it
   * $\mathbf C \mathrel{+}= \mathbf A^T \cdot \mathbf B^T$.
   *
   * @note It is assumed that @p A and @p B have compatible sizes and that
   * @p C already has the right size.
   *
   * @note This function uses the BLAS function Xgemm.
   */
  void
  TmTmult(LAPACKFullMatrix<number> &      C,
          const LAPACKFullMatrix<number> &B,
          const bool                      adding = false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void
  TmTmult(FullMatrix<number> &            C,
          const LAPACKFullMatrix<number> &B,
          const bool                      adding = false) const;

  /**
   * Performs out-place transposition.
   * Matrix @p B should be appropriately sized.
   *
   * @note for complex number types, conjugate transpose will be performed.
   *
   * @note If deal.II is configured with Intel-MKL, `mkl_?omatcopy` will be used,
   * otherwise transposition is done element by element.
   */
  void
  transpose(LAPACKFullMatrix<number> &B) const;

  /**
   * Scale rows of this matrix by @p V . This is equivalent to premultiplication
   * with a diagonal matrix $\mathbf A\leftarrow {\rm diag}(\mathbf V)\mathbf
   * A$.
   */
  void
  scale_rows(const Vector<number> &V);

  /**
   * Compute the LU factorization of the matrix using LAPACK function Xgetrf.
   */
  void
  compute_lu_factorization();

  /**
   * Compute the Cholesky factorization of the matrix using LAPACK function
   * Xpotrf.
   *
   * @note The factorization is stored in the lower-triangular part of the matrix.
   */
  void
  compute_cholesky_factorization();

  /**
   * Estimate the reciprocal of the condition number $1/k(\mathbf A)$ in $L_1$
   * norm ($1/(||\mathbf A||_1 \, ||\mathbf A^{-1}||_1)$) of a symmetric
   * positive definite matrix using Cholesky factorization. This function can
   * only be called if the matrix is already factorized.
   *
   * @note The condition number $k(\mathbf A)$ can be used to estimate the numerical
   * error related to the matrix inversion or the solution of the
   * system of linear algebraic equations as
   * <code>error = std::numeric_limits<Number>::epsilon * k</code>.
   * Alternatively one can get the number of accurate digits
   * <code>std::floor(std::log10(k))</code>.
   *
   * @note The function computes reciprocal of the condition number to
   * avoid possible overflow if the matrix is nearly singular.
   *
   * @param[in] l1_norm Is the $L_1$ norm of the matrix before calling Cholesky
   * factorization. It can be obtained by calling l1_norm().
   */
  number
  reciprocal_condition_number(const number l1_norm) const;

  /**
   * Estimate the reciprocal of the condition number $1/k(\mathbf A)$ in $L_1$
   * norm for triangular matrices. The matrix has to have the
   * LAPACKSupport::Property set to either
   * LAPACKSupport::Property::upper_triangular or
   * LAPACKSupport::Property::lower_triangular, see set_property().
   */
  number
  reciprocal_condition_number() const;

  /**
   * Compute the determinant of a matrix. As it requires the LU factorization of
   * the matrix, this function can only be called after
   * compute_lu_factorization() has been called.
   */
  number
  determinant() const;

  /**
   * Compute $L_1$ norm.
   */
  number
  l1_norm() const;

  /**
   * Compute $L_\infty$ norm.
   */
  number
  linfty_norm() const;

  /**
   * Compute Frobenius norm
   */
  number
  frobenius_norm() const;

  /**
   * Compute trace of the matrix, i.e. the sum of the diagonal values.
   * Obviously, the matrix needs to be quadratic for this function.
   */
  number
  trace() const;

  /**
   * Invert the matrix by first computing an LU/Cholesky factorization with the
   * LAPACK function Xgetrf/Xpotrf and then building the actual inverse using
   * Xgetri/Xpotri.
   */
  void
  invert();

  /**
   * Solve the linear system with right hand side @p v and put the solution
   * back to @p v. The matrix should be either triangular or LU/Cholesky
   * factorization should be previously computed.
   *
   * The flag transposed indicates whether the solution of the transposed
   * system is to be performed.
   */
  void
  solve(Vector<number> &v, const bool transposed = false) const;

  /**
   * Same as above but for multiple right hand sides (as many as there
   * are columns in the matrix @p B).
   */
  void
  solve(LAPACKFullMatrix<number> &B, const bool transposed = false) const;

  /**
   * Compute eigenvalues of the matrix. After this routine has been called,
   * eigenvalues can be retrieved using the eigenvalue() function. The matrix
   * itself will be LAPACKSupport::unusable after this operation.
   *
   * The optional arguments allow to compute left and right eigenvectors as
   * well.
   *
   * Note that the function does not return the computed eigenvalues right
   * away since that involves copying data around between the output arrays of
   * the LAPACK functions and any return array. This is often unnecessary
   * since one may not be interested in all eigenvalues at once, but for
   * example only the extreme ones. In that case, it is cheaper to just have
   * this function compute the eigenvalues and have a separate function that
   * returns whatever eigenvalue is requested.
   *
   * @note Calls the LAPACK function Xgeev.
   */
  void
  compute_eigenvalues(const bool right_eigenvectors = false,
                      const bool left_eigenvectors  = false);

  /**
   * Compute eigenvalues and eigenvectors of a real symmetric matrix. Only
   * eigenvalues in the interval $(\rm{lower\_bound}, \rm{upper\_bound}]$ are
   * computed with the absolute tolerance $\rm abs\_accuracy$. An approximate
   * eigenvalue is accepted as converged when it is determined to lie in an
   * interval $[a,b]$ of width less than or equal to $\rm{abs\_accuracy} + eps *
   * \rm{max}(|a|,|b|)$, where $eps$ is the machine precision.  If
   * $\rm{abs\_accuracy}$ is less than or equal to zero, then
   * $eps\,|\mathbf{T}|_1$ will be used in its place, where $|\mathbf{T}|_1$ is
   * the 1-norm of the tridiagonal matrix obtained by reducing $\mathbf A$ to
   * tridiagonal form. Eigenvalues will be computed most accurately when
   * $\rm{abs\_accuracy}$ is set to twice the underflow threshold, not zero.
   * After this routine has been called, all eigenvalues in $(\rm{lower\_bound},
   * \rm{upper\_bound}]$ will be stored in eigenvalues and the corresponding
   * eigenvectors will be stored in the columns of eigenvectors, whose dimension
   * is set accordingly.
   *
   * @note Calls the LAPACK function Xsyevx.
   */
  void
  compute_eigenvalues_symmetric(const number        lower_bound,
                                const number        upper_bound,
                                const number        abs_accuracy,
                                Vector<number> &    eigenvalues,
                                FullMatrix<number> &eigenvectors);

  /**
   * Compute generalized eigenvalues and eigenvectors of a real generalized
   * symmetric eigenproblem of the form
   * - itype = 1: $\mathbf A \cdot \mathbf x=\lambda \mathbf B \cdot \mathbf x$
   * - itype = 2: $\mathbf A \cdot \mathbf B \cdot \mathbf x=\lambda \mathbf x$
   * - itype = 3: $\mathbf B \cdot \mathbf A \cdot \mathbf x=\lambda \mathbf x$
   *
   * where $\mathbf A$ is this matrix. $\mathbf A$
   * and $\mathbf B$ are assumed to be symmetric, and $\mathbf B$ has to be
   * positive definite. Only eigenvalues in the interval $(\rm{lower\_bound},
   * \rm{upper\_bound}]$ are computed with the absolute tolerance
   * $\rm{abs\_accuracy}$.  An approximate eigenvalue is accepted as converged
   * when it is determined to lie in an interval $[a,b]$ of width less than or
   * equal to $\rm{abs\_accuracy} + eps * \rm{max}( |a|,|b| )$, where $eps$ is
   * the machine precision. If $\rm{abs\_accuracy}$ is less than or equal to
   * zero, then $eps \, |\mathbf{T}|_1$ will be used in its place, where
   * $|\mathbf{T}|_1$ is the 1-norm of the tridiagonal matrix obtained by
   * reducing $\mathbf A$ to tridiagonal form. Eigenvalues will be computed most
   * accurately when $\rm{abs\_accuracy}$ is set to twice the underflow
   * threshold, not zero. After this routine has been called, all eigenvalues in
   * $(\rm{lower\_bound}, \rm{upper\_bound}]$ will be stored in eigenvalues and
   * the corresponding eigenvectors will be stored in eigenvectors, whose
   * dimension is set accordingly.
   *
   * @note Calls the LAPACK function Xsygvx.
   */
  void
  compute_generalized_eigenvalues_symmetric(
    LAPACKFullMatrix<number> &   B,
    const number                 lower_bound,
    const number                 upper_bound,
    const number                 abs_accuracy,
    Vector<number> &             eigenvalues,
    std::vector<Vector<number>> &eigenvectors,
    const types::blas_int        itype = 1);

  /**
   * Same as the other compute_generalized_eigenvalues_symmetric function
   * except that all eigenvalues are computed and the tolerance is set
   * automatically.  Note that this function does not return the computed
   * eigenvalues right away since that involves copying data around between
   * the output arrays of the LAPACK functions and any return array. This is
   * often unnecessary since one may not be interested in all eigenvalues at
   * once, but for example only the extreme ones. In that case, it is cheaper
   * to just have this function compute the eigenvalues and have a separate
   * function that returns whatever eigenvalue is requested. Eigenvalues can
   * be retrieved using the eigenvalue() function.  The number of computed
   * eigenvectors is equal to eigenvectors.size()
   *
   * @note Calls the LAPACK function Xsygv.
   */
  void
  compute_generalized_eigenvalues_symmetric(
    LAPACKFullMatrix<number> &   B,
    std::vector<Vector<number>> &eigenvectors,
    const types::blas_int        itype = 1);

  /**
   * Compute the singular value decomposition of the matrix using LAPACK
   * function Xgesdd.
   *
   * Requires that the #state is LAPACKSupport::matrix, fills the data members
   * #wr, #svd_u, and #svd_vt, and leaves the object in the #state
   * LAPACKSupport::svd.
   *
   * The singular value decomposition factorizes the provided matrix (A) into
   * three parts: U, sigma, and the transpose of V (V^T), such that A = U sigma
   * V^T. Sigma is a MxN matrix which contains the singular values of A on
   * the diagonal while all the other elements are zero. U is a MxM orthogonal
   * matrix containing the left singular vectors corresponding to the singular
   * values of A. V is a NxN orthonal matrix containing the right singular
   * vectors corresponding the the singular values of A.
   *
   * Note that the variable #svd_vt contains the tranpose of V and can be
   * accessed by get_svd_vt(), while U is accessed with get_svd_u().
   */
  void
  compute_svd();

  /**
   * Compute the inverse of the matrix by singular value decomposition.
   *
   * Requires that #state is either LAPACKSupport::matrix or
   * LAPACKSupport::svd. In the first case, this function calls compute_svd().
   * After this function, the object will have the #state
   * LAPACKSupport::inverse_svd.
   *
   * For a singular value decomposition, the inverse is simply computed by
   * replacing all singular values by their reciprocal values. If the matrix
   * does not have maximal rank, singular values 0 are not touched, thus
   * computing the minimal norm right inverse of the matrix.
   *
   * The parameter @p threshold determines, when a singular value should
   * be considered zero. It is the ratio of the smallest to the largest
   * nonzero singular value $s_{max}$. Thus, the inverses of all
   * singular values less than  $s_{max}/\rm{threshold}$ will
   * be set to zero.
   */
  void
  compute_inverse_svd(const double threshold = 0.);

  /**
   * Same as above but provide the size of the kernel instead of a threshold,
   * i.e. the @p kernel_size smallest eigenvalues.
   */
  void
  compute_inverse_svd_with_kernel(const unsigned int kernel_size);

  /**
   * Retrieve eigenvalue after compute_eigenvalues() was called.
   */
  std::complex<number>
  eigenvalue(const size_type i) const;

  /**
   * Retrieve singular values after compute_svd() or compute_inverse_svd() was
   * called.
   */
  number
  singular_value(const size_type i) const;

  /**
   * Retrieve the matrix #svd_u after compute_svd() or compute_inverse_svd() was
   * called.
   */
  inline const LAPACKFullMatrix<number> &
  get_svd_u() const;

  /**
   * Retrieve the matrix #svd_vt after compute_svd() or compute_inverse_svd()
   * was called.
   */
  inline const LAPACKFullMatrix<number> &
  get_svd_vt() const;

  /**
   * Print the matrix and allow formatting of entries.
   *
   * The parameters allow for a flexible setting of the output format:
   *
   * @param out This specifies the stream to write to.
   *
   * @param precision denotes the number of trailing digits.
   *
   * @param scientific is used to determine the number format, where
   * <code>scientific = false</code> means fixed point notation.
   *
   * @param width denotes the with of each column. A zero entry for
   * @p width makes the function compute a width, but it may be changed
   * to a positive value, if output is crude.
   *
   * @param zero_string specifies a string printed for zero entries.
   *
   * @param denominator Multiply the whole matrix by this common
   * denominator to get nicer numbers.
   *
   * @param threshold all entries with absolute value smaller than
   * this are considered zero.
   *
   * @note The entries stored resemble a matrix only if the state is either
   * LAPACKSupport::matrix or LAPACK::inverse_matrix. Otherwise, calling this
   * function is not allowed.
   */
  void
  print_formatted(std::ostream &     out,
                  const unsigned int precision   = 3,
                  const bool         scientific  = true,
                  const unsigned int width       = 0,
                  const char *       zero_string = " ",
                  const double       denominator = 1.,
                  const double       threshold   = 0.) const;

private:
  /**
   * Internal function to compute various norms.
   */
  number
  norm(const char type) const;

  /**
   * Since LAPACK operations notoriously change the meaning of the matrix
   * entries, we record the current state after the last operation here.
   */
  LAPACKSupport::State state;

  /**
   * Additional property of the matrix which may help to select more
   * efficient LAPACK functions.
   */
  LAPACKSupport::Property property;

  /**
   * The working array used for some LAPACK functions.
   */
  mutable std::vector<number> work;

  /**
   * Integer working array used for some LAPACK functions.
   */
  mutable std::vector<types::blas_int> iwork;

  /**
   * The vector storing the permutations applied for pivoting in the LU-
   * factorization.
   *
   * Also used as the scratch array IWORK for LAPACK functions needing it.
   */
  std::vector<types::blas_int> ipiv;

  /**
   * Workspace for calculating the inverse matrix from an LU factorization.
   */
  std::vector<number> inv_work;

  /**
   * Real parts of eigenvalues or the singular values. Filled by
   * compute_eigenvalues() or compute_svd().
   */
  std::vector<typename numbers::NumberTraits<number>::real_type> wr;

  /**
   * Imaginary parts of eigenvalues, or, in the complex scalar case, the
   * eigenvalues themselves. Filled by compute_eigenvalues.
   */
  std::vector<number> wi;

  /**
   * Space where left eigenvectors can be stored.
   */
  std::vector<number> vl;

  /**
   * Space where right eigenvectors can be stored.
   */
  std::vector<number> vr;

  /**
   * The matrix $\mathbf U$ in the singular value decomposition
   * $\mathbf U \cdot \mathbf S \cdot \mathbf V^T$.
   */
  std::unique_ptr<LAPACKFullMatrix<number>> svd_u;

  /**
   * The matrix $\mathbf V^T$  in the singular value decomposition
   * $\mathbf U \cdot \mathbf S \cdot \mathbf V^T$.
   */
  std::unique_ptr<LAPACKFullMatrix<number>> svd_vt;

  /**
   * Thread mutex.
   */
  mutable std::mutex mutex;
};



/**
 * A preconditioner based on the LU-factorization of LAPACKFullMatrix.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat, 2006
 */
template <typename number>
class PreconditionLU : public Subscriptor
{
public:
  void
  initialize(const LAPACKFullMatrix<number> &);
  void
  initialize(const LAPACKFullMatrix<number> &, VectorMemory<Vector<number>> &);
  void
  vmult(Vector<number> &, const Vector<number> &) const;
  void
  Tvmult(Vector<number> &, const Vector<number> &) const;
  void
  vmult(BlockVector<number> &, const BlockVector<number> &) const;
  void
  Tvmult(BlockVector<number> &, const BlockVector<number> &) const;

private:
  SmartPointer<const LAPACKFullMatrix<number>, PreconditionLU<number>> matrix;
  SmartPointer<VectorMemory<Vector<number>>, PreconditionLU<number>>   mem;
};

/*---------------------- Inline functions -----------------------------------*/

template <typename number>
inline void
LAPACKFullMatrix<number>::set(const size_type i,
                              const size_type j,
                              const number    value)
{
  (*this)(i, j) = value;
}


template <typename number>
inline typename LAPACKFullMatrix<number>::size_type
LAPACKFullMatrix<number>::m() const
{
  return static_cast<size_type>(this->n_rows());
}

template <typename number>
inline typename LAPACKFullMatrix<number>::size_type
LAPACKFullMatrix<number>::n() const
{
  return static_cast<size_type>(this->n_cols());
}

template <typename number>
template <typename MatrixType>
inline void
LAPACKFullMatrix<number>::copy_from(const MatrixType &M)
{
  this->reinit(M.m(), M.n());

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // copy them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        this->el(row, entry->column()) = entry->value();
    }

  state = LAPACKSupport::matrix;
}



template <typename number>
template <typename MatrixType>
inline void
LAPACKFullMatrix<number>::fill(const MatrixType &M,
                               const size_type   dst_offset_i,
                               const size_type   dst_offset_j,
                               const size_type   src_offset_i,
                               const size_type   src_offset_j,
                               const number      factor,
                               const bool        transpose)
{
  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class
  for (size_type row = src_offset_i; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        {
          const size_type i = transpose ? entry->column() : row;
          const size_type j = transpose ? row : entry->column();

          const size_type dst_i = dst_offset_i + i - src_offset_i;
          const size_type dst_j = dst_offset_j + j - src_offset_j;
          if (dst_i < this->n_rows() && dst_j < this->n_cols())
            (*this)(dst_i, dst_j) = factor * entry->value();
        }
    }

  state = LAPACKSupport::matrix;
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::vmult(Vector<number2> &,
                                const Vector<number2> &,
                                const bool) const
{
  Assert(false,
         ExcMessage("LAPACKFullMatrix<number>::vmult must be called with a "
                    "matching Vector<double> vector type."));
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::vmult_add(Vector<number2> &,
                                    const Vector<number2> &) const
{
  Assert(false,
         ExcMessage("LAPACKFullMatrix<number>::vmult_add must be called with a "
                    "matching Vector<double> vector type."));
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::Tvmult(Vector<number2> &,
                                 const Vector<number2> &,
                                 const bool) const
{
  Assert(false,
         ExcMessage("LAPACKFullMatrix<number>::Tvmult must be called with a "
                    "matching Vector<double> vector type."));
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::Tvmult_add(Vector<number2> &,
                                     const Vector<number2> &) const
{
  Assert(false,
         ExcMessage(
           "LAPACKFullMatrix<number>::Tvmult_add must be called with a "
           "matching Vector<double> vector type."));
}


template <typename number>
inline std::complex<number>
LAPACKFullMatrix<number>::eigenvalue(const size_type i) const
{
  Assert(state & LAPACKSupport::eigenvalues, ExcInvalidState());
  Assert(wr.size() == this->n_rows(), ExcInternalError());
  Assert(wi.size() == this->n_rows(), ExcInternalError());
  AssertIndexRange(i, this->n_rows());

  if (numbers::NumberTraits<number>::is_complex)
    return std::complex<number>(wi[i]);
  else
    return std::complex<number>(wr[i], wi[i]);
}


template <typename number>
inline number
LAPACKFullMatrix<number>::singular_value(const size_type i) const
{
  Assert(state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd,
         LAPACKSupport::ExcState(state));
  AssertIndexRange(i, wr.size());

  return wr[i];
}


template <typename number>
inline const LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::get_svd_u() const
{
  Assert(state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd,
         LAPACKSupport::ExcState(state));

  return *svd_u;
}


template <typename number>
inline const LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::get_svd_vt() const
{
  Assert(state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd,
         LAPACKSupport::ExcState(state));

  return *svd_vt;
}



DEAL_II_NAMESPACE_CLOSE

#endif
