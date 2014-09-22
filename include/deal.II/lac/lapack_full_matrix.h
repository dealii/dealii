// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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

#ifndef __deal2__lapack_full_matrix_h
#define __deal2__lapack_full_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/table.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <vector>
#include <complex>

DEAL_II_NAMESPACE_OPEN

// forward declarations
template<typename number> class Vector;
template<typename number> class BlockVector;
template<typename number> class FullMatrix;
template<typename number> class SparseMatrix;


/**
 * A variant of FullMatrix using LAPACK functions wherever
 * possible. In order to do this, the matrix is stored in transposed
 * order. The element access functions hide this fact by reverting the
 * transposition.
 *
 * @note In order to perform LAPACK functions, the class contains a lot of
 * auxiliary data in the private section. The names of these data
 * vectors are usually the names chosen for the arguments in the
 * LAPACK documentation.
 *
 * @ingroup Matrix1
 * @author Guido Kanschat, 2005
 */
template <typename number>
class LAPACKFullMatrix : public TransposeTable<number>
{
public:

  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor. Initialize the matrix as a square matrix with dimension
   * <tt>n</tt>.
   *
   * In order to avoid the implicit conversion of integers and other types to
   * a matrix, this constructor is declared <tt>explicit</tt>.
   *
   * By default, no memory is allocated.
   */
  explicit LAPACKFullMatrix (const size_type size = 0);


  /**
   * Constructor. Initialize the matrix as a rectangular matrix.
   */
  LAPACKFullMatrix (const size_type rows,
                    const size_type cols);


  /**
   * Copy constructor. This constructor does a deep copy of the
   * matrix. Therefore, it poses a possible efficiency problem, if for
   * example, function arguments are passed by value rather than by
   * reference. Unfortunately, we can't mark this copy constructor
   * <tt>explicit</tt>, since that prevents the use of this class in
   * containers, such as <tt>std::vector</tt>. The responsibility to check
   * performance of programs must therefore remain with the user of this
   * class.
   */
  LAPACKFullMatrix (const LAPACKFullMatrix &);

  /**
   * Assignment operator.
   */
  LAPACKFullMatrix<number> &
  operator = (const LAPACKFullMatrix<number> &);

  /**
   * Assignment operator from a regular FullMatrix. @note Since LAPACK
   * expects matrices in transposed order, this transposition is
   * included here.
   */
  template <typename number2>
  LAPACKFullMatrix<number> &
  operator = (const FullMatrix<number2> &);

  /**
   * Assignment operator from a regular SparseMatrix. @note Since
   * LAPACK expects matrices in transposed order, this transposition
   * is included here.
   */
  template <typename number2>
  LAPACKFullMatrix<number> &
  operator = (const SparseMatrix<number2> &);

  /**
   * This operator assigns a scalar to a matrix. To avoid confusion with
   * constructors, zero is the only value allowed for <tt>d</tt>
   */
  LAPACKFullMatrix<number> &
  operator = (const double d);

  /**
   * Assignment from different matrix classes, performing the usual conversion
   * to the transposed format expected by LAPACK. This assignment operator
   * uses iterators of the class MATRIX. Therefore, sparse matrices are
   * possible sources.
   */
  template <class MATRIX>
  void copy_from (const MATRIX &);

  /**
   * Regenerate the current matrix by one that has the same properties
   * as if it were created by the constructor of this class with the
   * same argument list as this present function.
   */
  void reinit (const size_type size);

  /**
   * Regenerate the current matrix by one that has the same properties
   * as if it were created by the constructor of this class with the
   * same argument list as this present function.
   */
  void reinit (const size_type rows,
               const size_type cols);

  /**
   * Return the dimension of the range space. @note The matrix is of
   * dimension $m \times n$.
   */
  unsigned int m () const;

  /**
   * Return the number of the range space. @note The matrix is of
   * dimension $m \times n$.
   */
  unsigned int n () const;

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
  template<class MATRIX>
  void fill (const MATRIX &src,
             const size_type dst_offset_i = 0,
             const size_type dst_offset_j = 0,
             const size_type src_offset_i = 0,
             const size_type src_offset_j = 0,
             const number factor = 1.,
             const bool transpose = false);


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
   * finally with the left trandformation matrix.
   * </ul>
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>w</tt> or added to <tt>w</tt>.
   *
   * if (adding)
   *  <i>w += A*v</i>
   *
   * if (!adding)
   *  <i>w = A*v</i>
   *
   * @note Source and destination must not be the same vector.
   *
   * @note This template only exists for compile-time compatibility with
   * FullMatrix. Implementation is only available for
   * <tt>VECTOR=Vector&lt;number&gt;</tt>
   */
  template <class VECTOR>
  void vmult(VECTOR &dst, const VECTOR &src, const bool adding = false) const;

  /**
   * Adding Matrix-vector-multiplication.  <i>w += A*v</i>
   *
   * See the documentation of vmult() for details on the implementation.
   */
  template <class VECTOR>
  void vmult_add (VECTOR &w, const VECTOR &v) const;

  /**
   * Transpose matrix-vector-multiplication.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>w</tt> or added to <tt>w</tt>.
   *
   * if (adding)
   *  <i>w += A<sup>T</sup>*v</i>
   *
   * if (!adding)
   *  <i>w = A<sup>T</sup>*v</i>
   *
  * See the documentation of vmult() for details on the implementation.
   */
  template <class VECTOR>
  void Tvmult (VECTOR &w, const VECTOR &v,
               const bool            adding=false) const;

  /**
   * Adding transpose matrix-vector-multiplication.  <i>w +=
   * A<sup>T</sup>*v</i>
   *
   * See the documentation of vmult() for details on the implementation.
   */
  template <class VECTOR>
  void Tvmult_add (VECTOR &w, const VECTOR &v) const;

  void vmult (Vector<number>   &w,
              const Vector<number> &v,
              const bool            adding=false) const;
  void vmult_add (Vector<number>       &w,
                  const Vector<number> &v) const;
  void Tvmult (Vector<number>       &w,
               const Vector<number> &v,
               const bool            adding=false) const;
  void Tvmult_add (Vector<number>       &w,
                   const Vector<number> &v) const;


  /**
   * Matrix-matrix-multiplication.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding)
   *  <i>C += A*B</i>
   *
   * if (!adding)
   *  <i>C = A*B</i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm.
   */
  void mmult (LAPACKFullMatrix<number>       &C,
              const LAPACKFullMatrix<number> &B,
              const bool                      adding=false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void mmult (FullMatrix<number>             &C,
              const LAPACKFullMatrix<number> &B,
              const bool                      adding=false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt>.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding)
   *  <i>C += A<sup>T</sup>*B</i>
   *
   * if (!adding)
   *  <i>C = A<sup>T</sup>*B</i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm.
   */
  void Tmmult (LAPACKFullMatrix<number>       &C,
               const LAPACKFullMatrix<number> &B,
               const bool                      adding=false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void Tmmult (FullMatrix<number>             &C,
               const LAPACKFullMatrix<number> &B,
               const bool                      adding=false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>B</tt>.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding)
   *  <i>C += A*B<sup>T</sup></i>
   *
   * if (!adding)
   *  <i>C = A*B<sup>T</sup></i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm.
   */
  void mTmult (LAPACKFullMatrix<number>       &C,
               const LAPACKFullMatrix<number> &B,
               const bool                      adding=false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void mTmult (FullMatrix<number>             &C,
               const LAPACKFullMatrix<number> &B,
               const bool                      adding=false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt> and
   * <tt>B</tt>.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding)
   *  <i>C += A<sup>T</sup>*B<sup>T</sup></i>
   *
   * if (!adding)
   *  <i>C = A<sup>T</sup>*B<sup>T</sup></i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm.
   */
  void TmTmult (LAPACKFullMatrix<number>       &C,
                const LAPACKFullMatrix<number> &B,
                const bool                      adding=false) const;

  /**
   * Same as before, but stores the result in a FullMatrix, not in a
   * LAPACKFullMatrix.
   */
  void TmTmult (FullMatrix<number>             &C,
                const LAPACKFullMatrix<number> &B,
                const bool                      adding=false) const;

  /**
   * Compute the LU factorization of the matrix using LAPACK function Xgetrf.
   */
  void compute_lu_factorization ();

  /**
   * Invert the matrix by first computing an LU factorization with the LAPACK
   * function Xgetrf and then building the actual inverse using Xgetri.
   */
  void invert ();

  /**
   * Solve the linear system with right hand side given by applying
   * forward/backward substitution to the previously computed LU
   * factorization. Uses LAPACK function Xgetrs.
   *
   * The flag transposed indicates whether the solution of the transposed
   * system is to be performed.
   */
  void apply_lu_factorization (Vector<number> &v,
                               const bool      transposed) const;

  /**
   * Solve the linear system with multiple right hand sides (as many as there
   * are columns in the matrix b) given by applying forward/backward
   * substitution to the previously computed LU factorization. Uses LAPACK
   * function Xgetrs.
   *
   * The flag transposed indicates whether the solution of the transposed
   * system is to be performed.
   */
  void apply_lu_factorization (LAPACKFullMatrix<number> &B,
                               const bool                transposed) const;

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
  void compute_eigenvalues (const bool right_eigenvectors = false,
                            const bool left_eigenvectors  = false);

  /**
   * Compute eigenvalues and eigenvectors of a real symmetric matrix. Only
   * eigenvalues in the interval (lower_bound, upper_bound] are computed with
   * the absolute tolerance abs_accuracy. An approximate eigenvalue is
   * accepted as converged when it is determined to lie in an interval [a,b]
   * of width less than or equal to abs_accuracy + eps * max( |a|,|b| ), where
   * eps is the machine precision.  If abs_accuracy is less than or equal to
   * zero, then eps*|t| will be used in its place, where |t| is the 1-norm of
   * the tridiagonal matrix obtained by reducing A to tridiagonal form.
   * Eigenvalues will be computed most accurately when abs_accuracy is set to
   * twice the underflow threshold, not zero.  After this routine has been
   * called, all eigenvalues in (lower_bound, upper_bound] will be stored in
   * eigenvalues and the corresponding eigenvectors will be stored in the
   * columns of eigenvectors, whose dimension is set accordingly.
   *
   * @note Calls the LAPACK function Xsyevx. For this to work, ./configure has
   * to be told to use LAPACK.
   */
  void compute_eigenvalues_symmetric (const number        lower_bound,
                                      const number        upper_bound,
                                      const number        abs_accuracy,
                                      Vector<number>     &eigenvalues,
                                      FullMatrix<number> &eigenvectors);

  /**
   * Compute generalized eigenvalues and eigenvectors of a real generalized
   * symmetric eigenproblem of the form itype = 1: $Ax=\lambda B x$ itype = 2:
   * $ABx=\lambda x$ itype = 3: $BAx=\lambda x$, where A is this matrix.  A
   * and B are assumed to be symmetric, and B has to be positive definite.
   * Only eigenvalues in the interval (lower_bound, upper_bound] are computed
   * with the absolute tolerance abs_accuracy.  An approximate eigenvalue is
   * accepted as converged when it is determined to lie in an interval [a,b]
   * of width less than or equal to abs_accuracy + eps * max( |a|,|b| ), where
   * eps is the machine precision.  If abs_accuracy is less than or equal to
   * zero, then eps*|t| will be used in its place, where |t| is the 1-norm of
   * the tridiagonal matrix obtained by reducing A to tridiagonal form.
   * Eigenvalues will be computed most accurately when abs_accuracy is set to
   * twice the underflow threshold, not zero.  After this routine has been
   * called, all eigenvalues in (lower_bound, upper_bound] will be stored in
   * eigenvalues and the corresponding eigenvectors will be stored in
   * eigenvectors, whose dimension is set accordingly.
   *
   * @note Calls the LAPACK function Xsygvx. For this to work, ./configure has
   * to be told to use LAPACK.
   */
  void compute_generalized_eigenvalues_symmetric (LAPACKFullMatrix<number>     &B,
                                                  const number                  lower_bound,
                                                  const number                  upper_bound,
                                                  const number                  abs_accuracy,
                                                  Vector<number>               &eigenvalues,
                                                  std::vector<Vector<number> > &eigenvectors,
                                                  const int                     itype = 1);

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
   * @note Calls the LAPACK function Xsygv. For this to work, ./configure has
   * to be told to use LAPACK.
   */
  void compute_generalized_eigenvalues_symmetric (LAPACKFullMatrix<number>     &B,
                                                  std::vector<Vector<number> > &eigenvectors,
                                                  const int                     itype = 1);

  /**
   * Compute the singular value decomposition of the matrix using LAPACK
   * function Xgesdd.
   *
   * Requires that the #state is LAPACKSupport::matrix, fills the data members
   * #wr, #svd_u, and #svd_vt, and leaves the object in the #state
   * LAPACKSupport::svd.
   */
  void compute_svd ();

  /**
   * Compute the inverse of the matrix by singular value decomposition.
   *
   * Requires that #state is either LAPACKSupport::matrix or
   * LAPACKSupport::svd. In the first case, this function calls
   * compute_svd(). After this function, the object will have the #state
   * LAPACKSupport::inverse_svd.
   *
   * For a singular value decomposition, the inverse is simply computed by
   * replacing all singular values by their reciprocal values. If the matrix
   * does not have maximal rank, singular values 0 are not touched, thus
   * computing the minimal norm right inverse of the matrix.
   *
   * The parameter <tt>threshold</tt> determines, when a singular value should
   * be considered zero. It is the ratio of the smallest to the largest
   * nonzero singular value <i>s</i><sub>max</sub>. Thus, the inverses of all
   * singular values less than <i>s</i><sub>max</sub>/<tt>threshold</tt> will
   * be set to zero.
   */
  void compute_inverse_svd (const double threshold = 0.);

  /**
   * Retrieve eigenvalue after compute_eigenvalues() was called.
   */
  std::complex<number>
  eigenvalue (const size_type i) const;

  /**
   * Retrieve singular values after compute_svd() or compute_inverse_svd() was
   * called.
   */
  number
  singular_value (const size_type i) const;

  /**
   * Print the matrix and allow formatting of entries.
   *
   * The parameters allow for a flexible setting of the output format:
   *
   * @arg <tt>precision</tt> denotes the number of trailing digits.
   *
   * @arg <tt>scientific</tt> is used to determine the number format, where
   * <tt>scientific</tt> = <tt>false</tt> means fixed point notation.
   *
   * @arg <tt>width</tt> denotes the with of each column. A zero entry for
   * <tt>width</tt> makes the function compute a width, but it may be changed
   * to a positive value, if output is crude.
   *
   * @arg <tt>zero_string</tt> specifies a string printed for zero entries.
   *
   * @arg <tt>denominator</tt> Multiply the whole matrix by this common
   * denominator to get nicer numbers.
   *
   * @arg <tt>threshold</tt>: all entries with absolute value smaller than
   * this are considered zero.
  */
  void print_formatted (std::ostream       &out,
                        const unsigned int  presicion=3,
                        const bool          scientific  = true,
                        const unsigned int  width       = 0,
                        const char         *zero_string = " ",
                        const double        denominator = 1.,
                        const double        threshold   = 0.) const;

private:

  /**
   * Since LAPACK operations notoriously change the meaning of the
   * matrix entries, we record the current state after the last
   * operation here.
   */
  LAPACKSupport::State state;

  /**
   * Additional properties of the matrix which may help to select more
   * efficient LAPACK functions.
   */
  LAPACKSupport::Properties properties;

  /**
   * The working array used for some LAPACK functions.
   */
  mutable std::vector<number> work;

  /**
   * The vector storing the permutations applied for pivoting in the
   * LU-factorization.
   *
   * Also used as the scratch array IWORK for LAPACK functions needing it.
   */
  std::vector<int> ipiv;

  /**
   * Workspace for calculating the inverse matrix from an LU factorization.
   */
  std::vector<number> inv_work;

  /**
   * Real parts of eigenvalues or the singular values. Filled by
   * compute_eigenvalues() or compute_svd().
   */
  std::vector<number> wr;

  /**
   * Imaginary parts of eigenvalues. Filled by compute_eigenvalues.
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
   * The matrix <i>U</i> in the singular value decomposition
   * <i>USV<sup>T</sup></i>.
   */
  std_cxx11::shared_ptr<LAPACKFullMatrix<number> > svd_u;

  /**
   * The matrix <i>V<sup>T</sup></i> in the singular value decomposition
   * <i>USV<sup>T</sup></i>.
   */
  std_cxx11::shared_ptr<LAPACKFullMatrix<number> > svd_vt;
};



/**
 * A preconditioner based on the LU-factorization of LAPACKFullMatrix.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat, 2006
 */
template <typename number>
class PreconditionLU
  :
  public Subscriptor
{
public:
  void initialize(const LAPACKFullMatrix<number> &);
  void initialize(const LAPACKFullMatrix<number> &,
                  VectorMemory<Vector<number> > &);
  void vmult(Vector<number> &, const Vector<number> &) const;
  void Tvmult(Vector<number> &, const Vector<number> &) const;
  void vmult(BlockVector<number> &,
             const BlockVector<number> &) const;
  void Tvmult(BlockVector<number> &,
              const BlockVector<number> &) const;
private:
  SmartPointer<const LAPACKFullMatrix<number>,PreconditionLU<number> > matrix;
  SmartPointer<VectorMemory<Vector<number> >,PreconditionLU<number> > mem;
};

/*---------------------- Inline functions -----------------------------------*/

template <typename number>
inline
unsigned int
LAPACKFullMatrix<number>::m () const
{
  return this->n_rows ();
}

template <typename number>
inline
unsigned int
LAPACKFullMatrix<number>::n () const
{
  return this->n_cols ();
}

template <typename number>
template <class MATRIX>
inline void
LAPACKFullMatrix<number>::copy_from (const MATRIX &M)
{
  this->reinit (M.m(), M.n());

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // copy them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MATRIX::const_iterator end_row = M.end(row);
      for (typename MATRIX::const_iterator entry = M.begin(row);
           entry != end_row; ++entry)
        this->el(row, entry->column()) = entry->value();
    }

  state = LAPACKSupport::matrix;
}



template <typename number>
template <class MATRIX>
inline void
LAPACKFullMatrix<number>::fill (
  const MATRIX &M,
  const size_type dst_offset_i,
  const size_type dst_offset_j,
  const size_type src_offset_i,
  const size_type src_offset_j,
  const number factor,
  const bool transpose)
{
  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class
  for (size_type row = src_offset_i; row < M.m(); ++row)
    {
      const typename MATRIX::const_iterator end_row = M.end(row);
      for (typename MATRIX::const_iterator entry = M.begin(row);
           entry != end_row; ++entry)
        {
          const size_type i = transpose ? entry->column() : row;
          const size_type j = transpose ? row : entry->column();

          const size_type dst_i=dst_offset_i+i-src_offset_i;
          const size_type dst_j=dst_offset_j+j-src_offset_j;
          if (dst_i<this->n_rows() && dst_j<this->n_cols())
            (*this)(dst_i, dst_j) = factor * entry->value();
        }
    }

  state = LAPACKSupport::matrix;
}


template <typename number>
template <class VECTOR>
inline void
LAPACKFullMatrix<number>::vmult(VECTOR &, const VECTOR &, bool) const
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
template <class VECTOR>
inline void
LAPACKFullMatrix<number>::vmult_add(VECTOR &, const VECTOR &) const
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
template <class VECTOR>
inline void
LAPACKFullMatrix<number>::Tvmult(VECTOR &, const VECTOR &, bool) const
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
template <class VECTOR>
inline void
LAPACKFullMatrix<number>::Tvmult_add(VECTOR &, const VECTOR &) const
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
inline std::complex<number>
LAPACKFullMatrix<number>::eigenvalue (const size_type i) const
{
  Assert (state & LAPACKSupport::eigenvalues, ExcInvalidState());
  Assert (wr.size() == this->n_rows(), ExcInternalError());
  Assert (wi.size() == this->n_rows(), ExcInternalError());
  Assert (i<this->n_rows(), ExcIndexRange(i,0,this->n_rows()));

  return std::complex<number>(wr[i], wi[i]);
}


template <typename number>
inline number
LAPACKFullMatrix<number>::singular_value (const size_type i) const
{
  Assert (state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd, LAPACKSupport::ExcState(state));
  AssertIndexRange(i,wr.size());

  return wr[i];
}



DEAL_II_NAMESPACE_CLOSE

#endif
