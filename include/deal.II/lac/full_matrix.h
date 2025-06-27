// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_full_matrix_h
#define dealii_full_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/vector_operation.h>

#include <cstring>
#include <iomanip>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class LAPACKFullMatrix;
#endif

/**
 * @addtogroup Matrix1
 * @{
 */


/**
 * This class implements a rectangular, dense ("full") matrix, that is,
 * a matrix where we store every element whether it is zero or
 * nonzero.  This is in contrast to the SparseMatrix and related
 * classes, which store a "sparse" representation in which memory is
 * only allocated for those elements for which it is known that they
 * are (or, more precisely, could be) nonzero.
 *
 * The data type of the entries of the matrix is provided by the
 * template argument <tt>number</tt>.
 *
 * This class provides a number of operations that involve both the
 * current matrix and either another matrix or one or more
 * vectors. Each of these other objects may store their elements in a
 * different data type than the current object does (again, as
 * indicated by the corresponding template arguments of the other
 * objects). In these cases, internal calculations are usually done
 * with the accuracy of the vector argument, or with the more accurate
 * of the data types.
 *
 * @note Instantiations for this template are provided for <tt>@<float@>,
 * @<double@>, @<std::complex@<float@>@>,
 * @<std::complex@<double@>@></tt>. Others can be generated in application
 * programs, see
 * @ref Instantiations
 * for details.
 */
template <typename number>
class FullMatrix : public Table<2, number>
{
public:
  /**
   * This class only supports basic numeric types (i.e., we support double and
   * float but not automatically differentiated numbers).
   *
   * @note we test real_type here to get the underlying scalar type when using
   * std::complex.
   */
  static_assert(
    std::is_arithmetic_v<typename numbers::NumberTraits<number>::real_type>,
    "The FullMatrix class only supports basic numeric types. In particular, it "
    "does not support automatically differentiated numbers.");


  /**
   * A type of used to index into this container.
   */
  using size_type = std::size_t;

  /**
   * Type of matrix entries. This alias is analogous to <tt>value_type</tt>
   * in the standard library containers.
   */
  using value_type = number;

  /**
   * Use the base class mutable iterator type.
   */
  using iterator = typename Table<2, number>::iterator;

  /**
   * Use the base class constant iterator type.
   */
  using const_iterator = typename Table<2, number>::const_iterator;

  /**
   * Use the base class iterator functions.
   */
  using Table<2, number>::begin;

  /**
   * Use the base class iterator functions
   */
  using Table<2, number>::end;

  /**
   * Declare a type that has holds real-valued numbers with the same precision
   * as the template argument to this class. If the template argument of this
   * class is a real data type, then real_type equals the template argument.
   * If the template argument is a std::complex type then real_type equals the
   * type underlying the complex numbers.
   *
   * This alias is used to represent the return type of norms.
   */
  using real_type = typename numbers::NumberTraits<number>::real_type;

  /**
   * @name Constructors and initialization.  See also the base class Table.
   */
  /** @{ */

  /**
   * Constructor. Initialize the matrix as a square matrix with dimension
   * <tt>n</tt>.
   *
   * In order to avoid the implicit conversion of integers and other types to
   * a matrix, this constructor is declared <tt>explicit</tt>.
   *
   * By default, no memory is allocated.
   */
  explicit FullMatrix(const size_type n = 0);

  /**
   * Constructor. Initialize the matrix as a rectangular matrix.
   */
  FullMatrix(const size_type rows, const size_type cols);

  /**
   * Constructor initializing from an array of numbers. The array is arranged
   * line by line. No range checking is performed.
   */
  FullMatrix(const size_type rows, const size_type cols, const number *entries);

  /**
   * Construct a full matrix that equals the identity matrix of the size of
   * the argument. Using this constructor, one can easily create an identity
   * matrix of size <code>n</code> by saying
   * @code
   * FullMatrix<double> M(IdentityMatrix(n));
   * @endcode
   */
  FullMatrix(const IdentityMatrix &id);
  /**
   * @}
   */

  /**
   * @name Copying into and out of other matrices
   * @{
   */

  /**
   * Variable assignment operator.
   */
  template <typename number2>
  FullMatrix<number> &
  operator=(const FullMatrix<number2> &);

  /**
   * This operator assigns a scalar to a matrix. To avoid confusion with the
   * semantics of this function, zero is the only value allowed for
   * <tt>d</tt>, allowing you to clear a matrix in an intuitive way.
   *
   * @dealiiOperationIsMultithreaded
   */
  FullMatrix<number> &
  operator=(const number d);

  /**
   * Copy operator to create a full matrix that equals the identity matrix of
   * the size of the argument. This way, one can easily create an identity
   * matrix of size <code>n</code> by saying
   * @code
   *   M = IdentityMatrix(n);
   * @endcode
   */
  FullMatrix<number> &
  operator=(const IdentityMatrix &id);

  /**
   * Assignment operator for a LapackFullMatrix. The calling matrix must be of
   * the same size as the LAPACK matrix.
   */
  template <typename number2>
  FullMatrix<number> &
  operator=(const LAPACKFullMatrix<number2> &);


  /**
   * Assignment from different matrix classes. This assignment operator uses
   * iterators of the typename MatrixType. Therefore, sparse matrices are
   * possible sources.
   */
  template <typename MatrixType>
  void
  copy_from(const MatrixType &);

  /**
   * Transposing assignment from different matrix classes. This assignment
   * operator uses iterators of the typename MatrixType. Therefore, sparse
   * matrices are possible sources.
   */
  template <typename MatrixType>
  void
  copy_transposed(const MatrixType &);

  /**
   * Fill matrix with elements extracted from a tensor, taking rows included
   * between <tt>r_i</tt> and <tt>r_j</tt> and columns between <tt>c_i</tt>
   * and <tt>c_j</tt>. The resulting matrix is then inserted in the
   * destination matrix at position <tt>(dst_r, dst_c)</tt> Checks on the
   * indices are made.
   */
  template <int dim>
  void
  copy_from(const Tensor<2, dim> &T,
            const unsigned int    src_r_i = 0,
            const unsigned int    src_r_j = dim - 1,
            const unsigned int    src_c_i = 0,
            const unsigned int    src_c_j = dim - 1,
            const size_type       dst_r   = 0,
            const size_type       dst_c   = 0);

  /**
   * Insert a submatrix (also rectangular) into a tensor, putting its upper
   * left element at the specified position <tt>(dst_r, dst_c)</tt> and the
   * other elements consequently. Default values are chosen so that no
   * parameter needs to be specified if the size of the tensor and that of the
   * matrix coincide.
   */
  template <int dim>
  void
  copy_to(Tensor<2, dim>    &T,
          const size_type    src_r_i = 0,
          const size_type    src_r_j = dim - 1,
          const size_type    src_c_i = 0,
          const size_type    src_c_j = dim - 1,
          const unsigned int dst_r   = 0,
          const unsigned int dst_c   = 0) const;

  /**
   * Copy a subset of the rows and columns of another matrix into the current
   * object.
   *
   * @param matrix The matrix from which a subset is to be taken from.
   * @param row_index_set The set of rows of @p matrix from which to extract.
   * @param column_index_set The set of columns of @p matrix from which to
   * extract. @pre The number of elements in @p row_index_set and @p
   * column_index_set shall be equal to the number of rows and columns in the
   * current object. In other words, the current object is not resized for
   * this operation.
   */
  template <typename MatrixType, typename index_type>
  void
  extract_submatrix_from(const MatrixType              &matrix,
                         const std::vector<index_type> &row_index_set,
                         const std::vector<index_type> &column_index_set);

  /**
   * Copy the elements of the current matrix object into a specified set of
   * rows and columns of another matrix. Thus, this is a scatter operation.
   *
   * @param row_index_set The rows of @p matrix into which to write.
   * @param column_index_set The columns of @p matrix into which to write.
   * @param matrix The matrix within which certain elements are to be
   * replaced. @pre The number of elements in @p row_index_set and @p
   * column_index_set shall be equal to the number of rows and columns in the
   * current object. In other words, the current object is not resized for
   * this operation.
   */
  template <typename MatrixType, typename index_type>
  void
  scatter_matrix_to(const std::vector<index_type> &row_index_set,
                    const std::vector<index_type> &column_index_set,
                    MatrixType                    &matrix) const;

  /**
   * Fill rectangular block.
   *
   * A rectangular block of the matrix <tt>src</tt> is copied into
   * <tt>this</tt>. The upper left corner of the block being copied is
   * <tt>(src_offset_i,src_offset_j)</tt>.  The upper left corner of the
   * copied block is <tt>(dst_offset_i,dst_offset_j)</tt>.  The size of the
   * rectangular block being copied is the maximum size possible, determined
   * either by the size of <tt>this</tt> or <tt>src</tt>.
   */
  template <typename number2>
  void
  fill(const FullMatrix<number2> &src,
       const size_type            dst_offset_i = 0,
       const size_type            dst_offset_j = 0,
       const size_type            src_offset_i = 0,
       const size_type            src_offset_j = 0);


  /**
   * Make function of base class available.
   */
  template <typename number2>
  void
  fill(const number2 *);

  /**
   * Fill with permutation of another matrix.
   *
   * The matrix <tt>src</tt> is copied into the target. The two permutation
   * <tt>p_r</tt> and <tt>p_c</tt> operate in a way, such that <tt>result(i,j)
   * = src(p_r[i], p_c[j])</tt>.
   *
   * The vectors may also be a selection from a larger set of integers, if the
   * matrix <tt>src</tt> is bigger. It is also possible to duplicate rows or
   * columns by this method.
   */
  template <typename number2>
  void
  fill_permutation(const FullMatrix<number2>    &src,
                   const std::vector<size_type> &p_rows,
                   const std::vector<size_type> &p_cols);

  /**
   * Set a particular entry of the matrix to a value. Thus, calling
   * <code>A.set(1,2,3.141);</code> is entirely equivalent to the operation
   * <code>A(1,2) = 3.141;</code>. This function exists for compatibility with
   * the various sparse matrix objects.
   *
   * @param i The row index of the element to be set.
   * @param j The columns index of the element to be set.
   * @param value The value to be written into the element.
   */
  void
  set(const size_type i, const size_type j, const number value);
  /**
   * @}
   */
  /**
   * @name Non-modifying operators
   * @{
   */

  /**
   * Comparison operator. Be careful with this thing, it may eat up huge
   * amounts of computing time! It is most commonly used for internal
   * consistency checks of programs.
   */
  bool
  operator==(const FullMatrix<number> &) const;

  /**
   * Number of rows of this matrix.  Note that the matrix is of dimension <i>m
   * x n</i>.
   */
  size_type
  m() const;

  /**
   * Number of columns of this matrix.  Note that the matrix is of dimension
   * <i>m x n</i>.
   */
  size_type
  n() const;

  /**
   * Return whether the matrix contains only elements with value zero. This
   * function is mainly for internal consistency checks and should seldom be
   * used when not in debug mode since it uses quite some time.
   */
  bool
  all_zero() const;

  /**
   * Return the square of the norm of the vector <tt>v</tt> induced by this
   * matrix, i.e. <i>(v,Mv)</i>. This is useful, e.g. in the finite element
   * context, where the <i>L<sup>2</sup></i> norm of a function equals the
   * matrix norm with respect to the @ref GlossMassMatrix "mass matrix" of the vector representing
   * the nodal values of the finite element function.
   *
   * Obviously, the matrix needs to be quadratic for this operation, and for
   * the result to actually be a norm it also needs to be either real
   * symmetric or complex hermitian.
   *
   * The underlying template types of both this matrix and the given vector
   * should either both be real or complex-valued, but not mixed, for this
   * function to make sense.
   */
  template <typename number2>
  number2
  matrix_norm_square(const Vector<number2> &v) const;

  /**
   * Build the matrix scalar product <tt>u<sup>T</sup> M v</tt>. This function
   * is mostly useful when building the cellwise scalar product of two
   * functions in the finite element context.
   *
   * The underlying template types of both this matrix and the given vector
   * should either both be real or complex-valued, but not mixed, for this
   * function to make sense.
   */
  template <typename number2>
  number2
  matrix_scalar_product(const Vector<number2> &u,
                        const Vector<number2> &v) const;

  /**
   * Return the <i>l<sub>1</sub></i>-norm of the matrix, where $||M||_1 =
   * \max_j \sum_i |M_{ij}|$ (maximum of the sums over columns).
   */
  real_type
  l1_norm() const;

  /**
   * Return the $l_\infty$-norm of the matrix, where $||M||_\infty = \max_i
   * \sum_j |M_{ij}|$ (maximum of the sums over rows).
   */
  real_type
  linfty_norm() const;

  /**
   * Compute the Frobenius norm of the matrix.  Return value is the root of
   * the square sum of all matrix entries.
   *
   * @note For the timid among us: this norm is not the norm compatible with
   * the <i>l<sub>2</sub></i>-norm of the vector space.
   */
  real_type
  frobenius_norm() const;

  /**
   * Compute the relative norm of the skew-symmetric part. The return value is
   * the Frobenius norm of the skew-symmetric part of the matrix divided by
   * that of the matrix.
   *
   * Main purpose of this function is to check, if a matrix is symmetric
   * within a certain accuracy, or not.
   */
  real_type
  relative_symmetry_norm2() const;

  /**
   * Compute the determinant of a matrix.  This is only implemented for one,
   * two, and three dimensions, since for higher dimensions the numerical work
   * explodes.  Obviously, the matrix needs to be quadratic for this function.
   */
  number
  determinant() const;

  /**
   * Return the trace of the matrix, i.e. the sum of the diagonal values
   * (which happens to also equal the sum of the eigenvalues of a matrix).
   * Obviously, the matrix needs to be quadratic for this function.
   */
  number
  trace() const;

  /**
   * Output of the matrix in user-defined format given by the specified
   * precision and width. This function saves width and precision of the
   * stream before setting these given values for output, and restores the
   * previous values after output.
   */
  template <typename StreamType>
  void
  print(StreamType        &s,
        const unsigned int width     = 5,
        const unsigned int precision = 2) const;

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
   * @param separator specifies a string printed to separate row entries.
   */
  void
  print_formatted(std::ostream      &out,
                  const unsigned int precision   = 3,
                  const bool         scientific  = true,
                  const unsigned int width       = 0,
                  const char        *zero_string = " ",
                  const double       denominator = 1.,
                  const double       threshold   = 0.,
                  const char        *separator   = " ") const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  /** @} */
  /**
   * @name Iterator functions
   * @{
   */

  /**
   * Mutable iterator starting at the first entry of row <tt>r</tt>.
   */
  iterator
  begin(const size_type r);

  /**
   * One past the end mutable iterator of row <tt>r</tt>.
   */
  iterator
  end(const size_type r);

  /**
   * Constant iterator starting at the first entry of row <tt>r</tt>.
   */
  const_iterator
  begin(const size_type r) const;

  /**
   * One past the end constant iterator of row <tt>r</tt>.
   */
  const_iterator
  end(const size_type r) const;

  /** @} */
  /**
   * @name Modifying operators
   * @{
   */

  /**
   * Scale the entire matrix by a fixed factor.
   */
  FullMatrix &
  operator*=(const number factor);

  /**
   * Scale the entire matrix by the inverse of the given factor.
   */
  FullMatrix &
  operator/=(const number factor);

  /**
   * Simple addition of a scaled matrix, i.e. <tt>*this += a*A</tt>.
   *
   * The matrix <tt>A</tt> may be a full matrix over an arbitrary underlying
   * scalar type, as long as its data type is convertible to the data type of
   * this matrix.
   */
  template <typename number2>
  void
  add(const number a, const FullMatrix<number2> &A);

  /**
   * Multiple addition of scaled matrices, i.e. <tt>*this += a*A + b*B</tt>.
   *
   * The matrices <tt>A</tt> and <tt>B</tt> may be a full matrix over an
   * arbitrary underlying scalar type, as long as its data type is convertible
   * to the data type of this matrix.
   */
  template <typename number2>
  void
  add(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B);

  /**
   * Multiple addition of scaled matrices, i.e. <tt>*this += a*A + b*B +
   * c*C</tt>.
   *
   * The matrices <tt>A</tt>, <tt>B</tt> and <tt>C</tt> may be a full matrix
   * over an arbitrary underlying scalar type, as long as its data type is
   * convertible to the data type of this matrix.
   */
  template <typename number2>
  void
  add(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B,
      const number               c,
      const FullMatrix<number2> &C);

  /**
   * Add rectangular block.
   *
   * A rectangular block of the matrix <tt>src</tt> is added to <tt>this</tt>.
   * The upper left corner of the block being copied is
   * <tt>(src_offset_i,src_offset_j)</tt>.  The upper left corner of the
   * copied block is <tt>(dst_offset_i,dst_offset_j)</tt>.  The size of the
   * rectangular block being copied is the maximum size possible, determined
   * either by the size of <tt>this</tt> or <tt>src</tt> and the given
   * offsets.
   */
  template <typename number2>
  void
  add(const FullMatrix<number2> &src,
      const number               factor,
      const size_type            dst_offset_i = 0,
      const size_type            dst_offset_j = 0,
      const size_type            src_offset_i = 0,
      const size_type            src_offset_j = 0);

  /**
   * Weighted addition of the transpose of <tt>B</tt> to <tt>this</tt>.
   *
   * <i>A += s B<sup>T</sup></i>
   */
  template <typename number2>
  void
  Tadd(const number s, const FullMatrix<number2> &B);

  /**
   * Add transpose of a rectangular block.
   *
   * A rectangular block of the matrix <tt>src</tt> is transposed and
   * addedadded to <tt>this</tt>. The upper left corner of the block being
   * copied is <tt>(src_offset_i,src_offset_j)</tt> in the coordinates of the
   * <b>non</b>-transposed matrix.  The upper left corner of the copied block
   * is <tt>(dst_offset_i,dst_offset_j)</tt>.  The size of the rectangular
   * block being copied is the maximum size possible, determined either by the
   * size of <tt>this</tt> or <tt>src</tt>.
   */
  template <typename number2>
  void
  Tadd(const FullMatrix<number2> &src,
       const number               factor,
       const size_type            dst_offset_i = 0,
       const size_type            dst_offset_j = 0,
       const size_type            src_offset_i = 0,
       const size_type            src_offset_j = 0);

  /**
   * Add a single element at the given position.
   */
  void
  add(const size_type row, const size_type column, const number value);

  /**
   * Add an array of values given by <tt>values</tt> in the given global
   * matrix row at columns specified by col_indices in the full matrix. This
   * function is present for compatibility with the various sparse matrices in
   * deal.II. In particular, the two boolean fields @p elide_zero_values and
   * @p col_indices_are_sorted do not impact the performance of this routine,
   * as opposed to the sparse matrix case and are indeed ignored in the
   * implementation.
   */
  template <typename number2, typename index_type>
  void
  add(const size_type   row,
      const size_type   n_cols,
      const index_type *col_indices,
      const number2    *values,
      const bool        elide_zero_values      = true,
      const bool        col_indices_are_sorted = false);

  /**
   * <i>A(i,1...n) += s*A(j,1...n)</i>.  Simple addition of rows of this
   */
  void
  add_row(const size_type i, const number s, const size_type j);

  /**
   * <i>A(i,1...n) += s*A(j,1...n) + t*A(k,1...n)</i>.  Multiple addition of
   * rows of this.
   */
  void
  add_row(const size_type i,
          const number    s,
          const size_type j,
          const number    t,
          const size_type k);

  /**
   * <i>A(1...n,i) += s*A(1...n,j)</i>.  Simple addition of columns of this.
   */
  void
  add_col(const size_type i, const number s, const size_type j);

  /**
   * <i>A(1...n,i) += s*A(1...n,j) + t*A(1...n,k)</i>.  Multiple addition of
   * columns of this.
   */
  void
  add_col(const size_type i,
          const number    s,
          const size_type j,
          const number    t,
          const size_type k);

  /**
   * Swap <i>A(i,1...n) <-> A(j,1...n)</i>.  Swap rows i and j of this
   */
  void
  swap_row(const size_type i, const size_type j);

  /**
   * Swap <i>A(1...n,i) <-> A(1...n,j)</i>.  Swap columns i and j of this
   */
  void
  swap_col(const size_type i, const size_type j);

  /**
   * Add constant to diagonal elements of this, i.e. add a multiple of the
   * identity matrix.
   */
  void
  diagadd(const number s);

  /**
   * Assignment <tt>*this = a*A</tt>.
   */
  template <typename number2>
  void
  equ(const number a, const FullMatrix<number2> &A);

  /**
   * Assignment <tt>*this = a*A + b*B</tt>.
   */
  template <typename number2>
  void
  equ(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B);

  /**
   * Assignment <tt>*this = a*A + b*B + c*C</tt>.
   */
  template <typename number2>
  void
  equ(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B,
      const number               c,
      const FullMatrix<number2> &C);

  /**
   * Symmetrize the matrix by forming the mean value between the existing
   * matrix and its transpose, <i>A = 1/2(A+A<sup>T</sup>)</i>.
   *
   * Obviously the matrix must be quadratic for this operation.
   */
  void
  symmetrize();

  /**
   * A=Inverse(A). A must be a square matrix. Inversion of this matrix by
   * Gauss-Jordan algorithm with partial pivoting. This process is
   * well-behaved for positive definite matrices, but be aware of round-off
   * errors in the indefinite case.
   *
   * In case deal.II was configured with LAPACK, the functions Xgetrf and
   * Xgetri build an LU factorization and invert the matrix upon that
   * factorization, providing best performance up to matrices with a few
   * hundreds rows and columns.
   *
   * The numerical effort to invert an $n \times n$ matrix is of the order
   * $n^3$.
   */
  void
  gauss_jordan();

  /**
   * Assign the inverse of the given matrix to <tt>*this</tt>. This function
   * is hardcoded for quadratic matrices of dimension one to four. However,
   * since the amount of code needed grows quickly, the method gauss_jordan()
   * is invoked implicitly if the dimension is larger.
   */
  template <typename number2>
  void
  invert(const FullMatrix<number2> &M);

  /**
   * Assign the Cholesky decomposition $A=:L L^T$ of the given matrix $A$ to
   * <tt>*this</tt>, where $L$ is lower triangular matrix. The given matrix must
   * be symmetric positive definite.
   *
   * ExcMatrixNotPositiveDefinite will be thrown in the case that the matrix
   * is not positive definite.
   */
  template <typename number2>
  void
  cholesky(const FullMatrix<number2> &A);

  /**
   * <tt>*this(i,j)</tt> = $V(i) W(j)$ where $V,W$ are vectors of the same
   * length.
   */
  template <typename number2>
  void
  outer_product(const Vector<number2> &V, const Vector<number2> &W);

  /**
   * Assign the left_inverse of the given matrix to <tt>*this</tt>. The
   * calculation being performed is <i>(A<sup>T</sup>*A)<sup>-1</sup>
   * *A<sup>T</sup></i>.
   */
  template <typename number2>
  void
  left_invert(const FullMatrix<number2> &M);

  /**
   * Assign the right_inverse of the given matrix to <tt>*this</tt>. The
   * calculation being performed is <i>A<sup>T</sup>*(A*A<sup>T</sup>)
   * <sup>-1</sup></i>.
   */
  template <typename number2>
  void
  right_invert(const FullMatrix<number2> &M);

  /** @} */
  /**
   * @name Multiplications
   * @{
   */

  /**
   * Matrix-matrix-multiplication.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding) <i>C += A*B</i>
   *
   * if (!adding) <i>C = A*B</i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm if the product of the three
   * matrix dimensions is larger than 300 and BLAS was detected during
   * configuration. Using BLAS usually results in considerable performance
   * gains.
   */
  template <typename number2>
  void
  mmult(FullMatrix<number2>       &C,
        const FullMatrix<number2> &B,
        const bool                 adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt>.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding) <i>C += A<sup>T</sup>*B</i>
   *
   * if (!adding) <i>C = A<sup>T</sup>*B</i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm if the product of the three
   * matrix dimensions is larger than 300 and BLAS was detected during
   * configuration. Using BLAS usually results in considerable performance
   * gains.
   */
  template <typename number2>
  void
  Tmmult(FullMatrix<number2>       &C,
         const FullMatrix<number2> &B,
         const bool                 adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>B</tt>.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding) <i>C += A*B<sup>T</sup></i>
   *
   * if (!adding) <i>C = A*B<sup>T</sup></i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm if the product of the three
   * matrix dimensions is larger than 300 and BLAS was detected during
   * configuration. Using BLAS usually results in considerable performance
   * gains.
   */
  template <typename number2>
  void
  mTmult(FullMatrix<number2>       &C,
         const FullMatrix<number2> &B,
         const bool                 adding = false) const;

  /**
   * Matrix-matrix-multiplication using transpose of <tt>this</tt> and
   * <tt>B</tt>.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>C</tt> or added to <tt>C</tt>.
   *
   * if (adding) <i>C += A<sup>T</sup>*B<sup>T</sup></i>
   *
   * if (!adding) <i>C = A<sup>T</sup>*B<sup>T</sup></i>
   *
   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
   * <tt>C</tt> already has the right size.
   *
   * This function uses the BLAS function Xgemm if the product of the three
   * matrix dimensions is larger than 300 and BLAS was detected during
   * configuration. Using BLAS usually results in considerable performance
   * gains.
   */
  template <typename number2>
  void
  TmTmult(FullMatrix<number2>       &C,
          const FullMatrix<number2> &B,
          const bool                 adding = false) const;

  /**
   * Add to the current matrix the triple product <b>B A D</b>. Optionally,
   * use the transposes of the matrices <b>B</b> and <b>D</b>. The scaling
   * factor scales the whole product, which is helpful when adding a multiple
   * of the triple product to the matrix.
   *
   * This product was written with the Schur complement <b>B<sup>T</sup>
   * A<sup>-1</sup> D</b> in mind.  Note that in this case the argument for
   * <tt>A</tt> must be the inverse of the matrix <b>A</b>.
   */
  void
  triple_product(const FullMatrix<number> &A,
                 const FullMatrix<number> &B,
                 const FullMatrix<number> &D,
                 const bool                transpose_B = false,
                 const bool                transpose_D = false,
                 const number              scaling     = number(1.));


  /**
   * @brief Compute the Kronecker product of two matrices.
   *
   * This function computes the Kronecker product of two matrices A and B, and
   * stores the result in the current matrix. The Kronecker product of an
   * m x n matrix A and a p x q matrix B is an (m*p) x (n*q) matrix defined as:
   *
   * ```
   * A ⊗ B = | a11*B   a12*B   ...   a1n*B |
   *         | a21*B   a22*B   ...   a2n*B |
   *         | ...     ...     ...   ...   |
   *         | am1*B   am2*B   ...   amn*B |
   * ```
   *
   * where aij are the elements of the matrix A.
   *
   * @param A The first matrix (m x n).
   * @param B The second matrix (p x q).
   * @param adding If `true`, the result is added to the current matrix. If
   *               `false` (default), the current matrix is overwritten.
   */
  void
  kronecker_product(const FullMatrix<number> &A,
                    const FullMatrix<number> &B,
                    const bool                adding = false);

  /**
   * Matrix-vector-multiplication.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>w</tt> or added to <tt>w</tt>.
   *
   * if (adding) <i>w += A*v</i>
   *
   * if (!adding) <i>w = A*v</i>
   *
   * Source and destination must not be the same vector.
   */
  template <typename number2>
  void
  vmult(Vector<number2>       &w,
        const Vector<number2> &v,
        const bool             adding = false) const;

  /**
   * Adding Matrix-vector-multiplication.  <i>w += A*v</i>
   *
   * Source and destination must not be the same vector.
   */
  template <typename number2>
  void
  vmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * Transpose matrix-vector-multiplication.
   *
   * The optional parameter <tt>adding</tt> determines, whether the result is
   * stored in <tt>w</tt> or added to <tt>w</tt>.
   *
   * if (adding) <i>w += A<sup>T</sup>*v</i>
   *
   * if (!adding) <i>w = A<sup>T</sup>*v</i>
   *
   *
   * Source and destination must not be the same vector.
   */
  template <typename number2>
  void
  Tvmult(Vector<number2>       &w,
         const Vector<number2> &v,
         const bool             adding = false) const;

  /**
   * Adding transpose matrix-vector-multiplication.  <i>w +=
   * A<sup>T</sup>*v</i>
   *
   * Source and destination must not be the same vector.
   */
  template <typename number2>
  void
  Tvmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * Apply the Jacobi preconditioner, which multiplies every element of the
   * <tt>src</tt> vector by the inverse of the respective diagonal element and
   * multiplies the result with the damping factor <tt>omega</tt>.
   */
  template <typename somenumber>
  void
  precondition_Jacobi(Vector<somenumber>       &dst,
                      const Vector<somenumber> &src,
                      const number              omega = 1.) const;

  /**
   * <i>dst=b-A*x</i>. Residual calculation, returns the
   * <i>l<sub>2</sub></i>-norm |<i>dst</i>|.
   *
   * Source <i>x</i> and destination <i>dst</i> must not be the same vector.
   */
  template <typename number2, typename number3>
  number
  residual(Vector<number2>       &dst,
           const Vector<number2> &x,
           const Vector<number3> &b) const;

  /**
   * Forward elimination of lower triangle.  Inverts the lower triangle of a
   * rectangular matrix for a given right hand side.
   *
   * If the matrix has more columns than rows, this function only operates on
   * the left quadratic submatrix. If there are more rows, the upper quadratic
   * part of the matrix is considered.
   *
   * @note It is safe to use the same object for @p dst and @p src.
   */
  template <typename number2>
  void
  forward(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * Backward elimination of upper triangle.
   *
   * See forward()
   *
   * @note It is safe to use the same object for @p dst and @p src.
   */
  template <typename number2>
  void
  backward(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * Dummy function for compatibility with distributed, parallel matrices.
   */
  void compress(VectorOperation::values);

  /** @} */

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException0(ExcEmptyMatrix);

  /**
   * Exception
   */
  DeclException1(
    ExcNotRegular,
    number,
    << "The maximal pivot is " << arg1
    << ", which is below the threshold. The matrix may be singular.");
  /**
   * Exception
   */
  DeclException3(ExcInvalidDestination,
                 size_type,
                 size_type,
                 size_type,
                 << "Target region not in matrix: size in this direction="
                 << arg1 << ", size of new matrix=" << arg2
                 << ", offset=" << arg3);
  /**
   * Exception
   */
  DeclExceptionMsg(ExcSourceEqualsDestination,
                   "You are attempting an operation on two vectors that "
                   "are the same object, but the operation requires that the "
                   "two objects are in fact different.");
  /**
   * Exception
   */
  DeclException0(ExcMatrixNotPositiveDefinite);
  /** @} */
};

/** @} */

#ifndef DOXYGEN
/*-------------------------Inline functions -------------------------------*/



template <typename number>
inline typename FullMatrix<number>::size_type
FullMatrix<number>::m() const
{
  return this->n_rows();
}



template <typename number>
inline typename FullMatrix<number>::size_type
FullMatrix<number>::n() const
{
  return this->n_cols();
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator=(const number d)
{
  Assert(d == number(0), ExcScalarAssignmentOnlyForZeroValue());

  if (this->n_elements() != 0)
    this->reset_values();

  return *this;
}



template <typename number>
template <typename number2>
inline void
FullMatrix<number>::fill(const number2 *src)
{
  Table<2, number>::fill(src);
}



template <typename number>
template <typename MatrixType>
void
FullMatrix<number>::copy_from(const MatrixType &M)
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
}



template <typename number>
template <int dim>
void
FullMatrix<number>::copy_from(const Tensor<2, dim> &T,
                              const unsigned int    src_r_i,
                              const unsigned int    src_r_j,
                              const unsigned int    src_c_i,
                              const unsigned int    src_c_j,
                              const size_type       dst_r,
                              const size_type       dst_c)
{
  Assert(!this->empty(), ExcEmptyMatrix());
  AssertIndexRange(src_r_j - src_r_i, this->m() - dst_r);
  AssertIndexRange(src_c_j - src_c_i, this->n() - dst_c);
  AssertIndexRange(src_r_j, dim);
  AssertIndexRange(src_c_j, dim);
  AssertIndexRange(src_r_i, src_r_j + 1);
  AssertIndexRange(src_c_i, src_c_j + 1);

  for (size_type i = 0; i < src_r_j - src_r_i + 1; ++i)
    for (size_type j = 0; j < src_c_j - src_c_i + 1; ++j)
      {
        const unsigned int src_r_index = static_cast<unsigned int>(i + src_r_i);
        const unsigned int src_c_index = static_cast<unsigned int>(j + src_c_i);
        (*this)(i + dst_r, j + dst_c)  = number(T[src_r_index][src_c_index]);
      }
}



template <typename number>
template <int dim>
void
FullMatrix<number>::copy_to(Tensor<2, dim>    &T,
                            const size_type    src_r_i,
                            const size_type    src_r_j,
                            const size_type    src_c_i,
                            const size_type    src_c_j,
                            const unsigned int dst_r,
                            const unsigned int dst_c) const
{
  Assert(!this->empty(), ExcEmptyMatrix());
  AssertIndexRange(src_r_j - src_r_i, dim - dst_r);
  AssertIndexRange(src_c_j - src_c_i, dim - dst_c);
  AssertIndexRange(src_r_j, this->m());
  AssertIndexRange(src_r_j, this->n());
  AssertIndexRange(src_r_i, src_r_j + 1);
  AssertIndexRange(src_c_j, src_c_j + 1);

  for (size_type i = 0; i < src_r_j - src_r_i + 1; ++i)
    for (size_type j = 0; j < src_c_j - src_c_i + 1; ++j)
      {
        const unsigned int dst_r_index = static_cast<unsigned int>(i + dst_r);
        const unsigned int dst_c_index = static_cast<unsigned int>(j + dst_c);
        T[dst_r_index][dst_c_index] = double((*this)(i + src_r_i, j + src_c_i));
      }
}



template <typename number>
template <typename MatrixType>
void
FullMatrix<number>::copy_transposed(const MatrixType &M)
{
  this->reinit(M.n(), M.m());

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // copy them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        this->el(entry->column(), row) = entry->value();
    }
}



template <typename number>
template <typename MatrixType, typename index_type>
inline void
FullMatrix<number>::extract_submatrix_from(
  const MatrixType              &matrix,
  const std::vector<index_type> &row_index_set,
  const std::vector<index_type> &column_index_set)
{
  AssertDimension(row_index_set.size(), this->n_rows());
  AssertDimension(column_index_set.size(), this->n_cols());

  const size_type n_rows_submatrix = row_index_set.size();
  const size_type n_cols_submatrix = column_index_set.size();

  for (size_type sub_row = 0; sub_row < n_rows_submatrix; ++sub_row)
    for (size_type sub_col = 0; sub_col < n_cols_submatrix; ++sub_col)
      (*this)(sub_row, sub_col) =
        matrix.el(row_index_set[sub_row], column_index_set[sub_col]);
}



template <typename number>
template <typename MatrixType, typename index_type>
inline void
FullMatrix<number>::scatter_matrix_to(
  const std::vector<index_type> &row_index_set,
  const std::vector<index_type> &column_index_set,
  MatrixType                    &matrix) const
{
  AssertDimension(row_index_set.size(), this->n_rows());
  AssertDimension(column_index_set.size(), this->n_cols());

  const size_type n_rows_submatrix = row_index_set.size();
  const size_type n_cols_submatrix = column_index_set.size();

  for (size_type sub_row = 0; sub_row < n_rows_submatrix; ++sub_row)
    for (size_type sub_col = 0; sub_col < n_cols_submatrix; ++sub_col)
      matrix.set(row_index_set[sub_row],
                 column_index_set[sub_col],
                 (*this)(sub_row, sub_col));
}


template <typename number>
inline void
FullMatrix<number>::set(const size_type i,
                        const size_type j,
                        const number    value)
{
  (*this)(i, j) = value;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::vmult_add(Vector<number2>       &w,
                              const Vector<number2> &v) const
{
  vmult(w, v, true);
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::Tvmult_add(Vector<number2>       &w,
                               const Vector<number2> &v) const
{
  Tvmult(w, v, true);
}


//---------------------------------------------------------------------------
template <typename number>
inline typename FullMatrix<number>::iterator
FullMatrix<number>::begin(const size_type r)
{
  AssertIndexRange(r, m());
  return iterator(this, r, 0);
}



template <typename number>
inline typename FullMatrix<number>::iterator
FullMatrix<number>::end(const size_type r)
{
  AssertIndexRange(r, m());
  return iterator(this, r + 1, 0);
}



template <typename number>
inline typename FullMatrix<number>::const_iterator
FullMatrix<number>::begin(const size_type r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r, 0);
}



template <typename number>
inline typename FullMatrix<number>::const_iterator
FullMatrix<number>::end(const size_type r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r + 1, 0);
}



template <typename number>
inline void
FullMatrix<number>::add(const size_type r, const size_type c, const number v)
{
  AssertIndexRange(r, this->m());
  AssertIndexRange(c, this->n());

  this->operator()(r, c) += v;
}



template <typename number>
template <typename number2, typename index_type>
inline void
FullMatrix<number>::add(const size_type   row,
                        const size_type   n_cols,
                        const index_type *col_indices,
                        const number2    *values,
                        const bool,
                        const bool)
{
  AssertIndexRange(row, this->m());
  for (size_type col = 0; col < n_cols; ++col)
    {
      AssertIndexRange(col_indices[col], this->n());
      this->operator()(row, col_indices[col]) += values[col];
    }
}


template <typename number>
template <typename StreamType>
inline void
FullMatrix<number>::print(StreamType        &s,
                          const unsigned int w,
                          const unsigned int p) const
{
  Assert(!this->empty(), ExcEmptyMatrix());

  // save the state of out stream
  const std::streamsize old_precision = s.precision(p);
  const std::streamsize old_width     = s.width(w);

  for (size_type i = 0; i < this->m(); ++i)
    {
      for (size_type j = 0; j < this->n(); ++j)
        {
          s.width(w);
          s.precision(p);
          s << this->el(i, j);
        }
      s << std::endl;
    }

  // reset output format
  s.precision(old_precision);
  s.width(old_width);
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
