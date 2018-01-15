// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_scalapack_h
#define dealii_scalapack_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SCALAPACK

#include <deal.II/base/exceptions.h>
#include <deal.II/base/process_grid.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/thread_management.h>
#include <mpi.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * A wrapper class around ScaLAPACK parallel dense linear algebra.
 *
 * ScaLAPACK assumes that matrices are distributed according to the
 * block-cyclic decomposition scheme. An $M$ by $N$ matrix is first decomposed
 * into $MB$ by $NB$ blocks which are then uniformly distributed across
 * the 2D process grid $p*q \le Np$, where $p,q$ are grid dimensions and
 * $Np$ is the total number of processes.
 *
 * For example, a global real symmetric matrix of size $9\times 9$ is stored in
 * upper storage mode with block sizes 4 × 4:
 * @code
 *                0                       1                2
 *     ┌                                                       ┐
 *     | -6.0  0.0  0.0  0.0  |   0.0 -2.0 -2.0  0.0  |   -2.0 |
 *     |   .  -6.0 -2.0  0.0  |  -2.0 -4.0  0.0 -4.0  |   -2.0 |
 * 0   |   .    .  -6.0 -2.0  |  -2.0  0.0  2.0  0.0  |    6.0 |
 *     |   .    .    .  -6.0  |   2.0  0.0  2.0  0.0  |    2.0 |
 *     | ---------------------|-----------------------|------- |
 *     |   .    .    .    .   |  -8.0 -4.0  0.0 -2.0  |    0.0 |
 *     |   .    .    .    .   |    .  -6.0  0.0 -4.0  |   -6.0 |
 * 1   |   .    .    .    .   |    .    .  -4.0  0.0  |    0.0 |
 *     |   .    .    .    .   |    .    .    .  -4.0  |   -4.0 |
 *     | ---------------------|-----------------------|------- |
 * 2   |   .    .    .    .   |    .    .    .    .   |  -16.0 |
 *     └                                                       ┘
 * @endcode
 * may be distributed using the 2x2 process grid:
 * @code
 *      |   0 2   |   1
 * -----| ------- |-----
 * 0    |   P00   |  P01
 * 2    |         |
 * -----| ------- |-----
 * 1    |   P10   |  P11
 * @endcode
 * with the following local arrays:
 * @code
 * p,q  |             0              |           1
 * -----|----------------------------|----------------------
 *      | -6.0  0.0  0.0  0.0  -2.0  |   0.0 -2.0 -2.0  0.0
 *      |   .  -6.0 -2.0  0.0  -2.0  |  -2.0 -4.0  0.0 -4.0
 *  0   |   .    .  -6.0 -2.0   6.0  |  -2.0  0.0  2.0  0.0
 *      |   .    .    .  -6.0   2.0  |   2.0  0.0  2.0  0.0
 *      |   .    .    .    .  -16.0  |    .    .    .    .
 * -----|----------------------------|----------------------
 *      |   .    .    .    .    0.0  |  -8.0 -4.0  0.0 -2.0
 *      |   .    .    .    .   -6.0  |    .  -6.0  0.0 -4.0
 *  1   |   .    .    .    .    0.0  |    .    .  -4.0  0.0
 *      |   .    .    .    .   -4.0  |    .    .    .  -4.0
 * @endcode
 * Note how processes $(0,0)$ and $(1,0)$ of the process grid store an
 * extra column to represent the last column of the original matrix that
 * did not fit the decomposition into $4\times 4$ sub-blocks.
 *
 * The choice of the block size is a compromise between a sufficiently large
 * size for efficient local/serial BLAS, but one that is also small enough to achieve
 * good parallel load balance.
 *
 * Below we show a strong scaling example of ScaLAPACKMatrix::invert()
 * on up to 5 nodes each composed of two Intel Xeon 2660v2 IvyBridge sockets
 * 2.20GHz, 10 cores/socket. Calculations are performed on square processor
 * grids 1x1, 2x2, 3x3, 4x4, 5x5, 6x6, 7x7, 8x8, 9x9, 10x10.
 *
 * @image html scalapack_invert.png
 *
 * @ingroup Matrix1
 * @author Denis Davydov, Benjamin Brands, 2017
 */
template <typename NumberType>
class ScaLAPACKMatrix : protected TransposeTable<NumberType>
{
public:

  /**
   * Declare the type for container size.
   */
  typedef unsigned int size_type;

  /**
   * Constructor for a rectangular matrix with @p n_rows and @p n_cols
   * and distributed using the grid @p process_grid.
   */
  ScaLAPACKMatrix(const size_type n_rows,
                  const size_type n_columns,
                  const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
                  const size_type row_block_size = 32,
                  const size_type column_block_size = 32,
                  const LAPACKSupport::Property property = LAPACKSupport::Property::general);

  /**
   * Constructor for a square matrix of size @p size, and distributed
   * using the process grid in @p process_grid.
   */
  ScaLAPACKMatrix(const size_type size,
                  const std::shared_ptr<const Utilities::MPI::ProcessGrid> process_grid,
                  const size_type block_size = 32,
                  const LAPACKSupport::Property property = LAPACKSupport::Property::symmetric);

  /**
   * Destructor
   */
  ~ScaLAPACKMatrix() = default;

  /**
   * Assign @p property to this matrix.
   */
  void set_property(const LAPACKSupport::Property property);

  /**
   * Assignment operator from a regular FullMatrix.
   *
   * @note This function should only be used for relatively small matrix
   * dimensions. It is primarily intended for debugging purposes.
   */
  ScaLAPACKMatrix<NumberType> &
  operator = (const FullMatrix<NumberType> &);

  /**
   * Copy the contents of the distributed matrix into @p matrix.
   *
   * @note This function should only be used for relatively small matrix
   * dimensions. It is primarily intended for debugging purposes.
   */
  void copy_to (FullMatrix<NumberType> &matrix) const;


  /**
   * Copy the contents of the distributed matrix into a differently distributed matrix @p dest.
   * The function also works for matrices with different process grids
   * or block-cyclic distributions.
   */
  void copy_to (ScaLAPACKMatrix<NumberType> &dest) const;


  /**
   * Compute the Cholesky factorization of the matrix using ScaLAPACK
   * function <code>pXpotrf</code>. The result of the factorization is stored in this object.
   */
  void compute_cholesky_factorization ();

  /**
   * Invert the matrix by first computing a Cholesky factorization and then
   * building the actual inverse using <code>pXpotri</code>. The inverse is stored
   * in this object.
   */
  void invert();



  /**
   * Function to compute selected eigenvalues and, optionally, the eigenvectors.
   * If the function is called with the default arguments all eigenvalues are computed but no eigenvectors.
   * The eigenvalues/eigenvectors are selected by either prescribing a range of indices @p index_limits
   * or a range of values @p value_limits for the eigenvalues. The funtion will throw an exception
   * if both ranges are prescribed (meaning that both ranges differ from the default value)
   * as this ambiguity is prohibited.
   * If successful, the computed eigenvalues are arranged in ascending order.
   * The eigenvectors are stored in the columns of the matrix, thereby
   * overwriting the original content of the matrix.
   */
  std::vector<NumberType> eigenpairs_symmetric(const bool compute_eigenvectors=false,
                                               const std::pair<int,int> &index_limits = std::make_pair(-1,-1),
                                               const std::pair<NumberType,NumberType> &value_limits = std::make_pair(-1,-1));



  /**
  * Funcion to compute the singular value decomposition (SVD) of an
  * M-by-N matrix A, optionally computing the left and/or right
  * singular vectors. The SVD is written as A = U * SIGMA * transpose(V)
  * where SIGMA is an M-by-N diagonal matrix, @p U is an M-by-M orthogonal matrix,
  * and @p V is an N-by-N orthogonal matrix. The diagonal elements of SIGMA
  * are the singular values of A and the columns of U and V are the
  * corresponding left and right singular vectors, respectively. The
  * singular values are returned in decreasing order and only the first min(M,N)
  * columns of U and rows of VT = transpose(V) are computed.
  * Upon return the content of the matrix is unusable.
  * The matrix A must have identical block cyclic distribution for the rows and column
  * If left singular vectors are required matrices A and U
  * have to be constructed with the same process grid and block cyclic distribution.
  * If right singular vectors are required matrices A and VT
  * have to be constructed with the same process grid  and block cyclic distribution.
   */
  std::vector<NumberType> compute_SVD(ScaLAPACKMatrix<NumberType> &U,
                                      ScaLAPACKMatrix<NumberType> &VT,
                                      const bool left_singluar_vectors=false,
                                      const bool right_singluar_vectors=false);

  /**
   * Estimate the the condition number of a SPD matrix in the $l_1$-norm.
   * The matrix has to be in the Cholesky state (see compute_cholesky_factorization()).
   * The reciprocal of the
   * condition number is returned in order to avoid the possibility of
   * overflow when the condition number is very large.
   *
   * @p a_norm must contain the $l_1$-norm of the matrix prior to calling
   * Cholesky factorization (see l1_norm()).
   *
   * @note An alternative is to compute the inverse of the matrix
   * explicitly and manually construct $k_1 = ||A||_1 ||A^{-1}||_1$.
   */
  NumberType reciprocal_condition_number(const NumberType a_norm) const;

  /**
   * Compute the $l_1$-norm of the matrix.
   */
  NumberType l1_norm() const;

  /**
   * Compute the $l_{\infty}$ norm of the matrix.
   */
  NumberType linfty_norm() const;

  /**
   * Compute the Frobenius norm of the matrix.
   */
  NumberType frobenius_norm() const;

  /**
   * Number of rows of the $M \times N$ matrix.
   */
  size_type m() const;

  /**
   * Number of columns of the $M \times N$ matrix.
   */
  size_type n() const;

  /**
   * Number of local rows on this MPI processes.
   */
  unsigned int local_m() const;

  /**
   * Number of local columns on this MPI process.
   */
  unsigned int local_n() const;

  /**
   * Return the global row number for the given local row @p loc_row .
   */
  unsigned int global_row(const unsigned int loc_row) const;

  /**
   * Return the global column number for the given local column @p loc_column.
   */
  unsigned int global_column(const unsigned int loc_column) const;

  /**
   * Read access to local element.
   */
  NumberType local_el(const unsigned int loc_row, const unsigned int loc_column) const;

  /**
   * Write access to local element.
   */
  NumberType &local_el(const unsigned int loc_row, const unsigned int loc_column);

private:

  /**
   * Calculate the norm of a distributed dense matrix using ScaLAPACK's
   * internal function.
   */
  NumberType norm(const char type) const;

  /**
   * Since ScaLAPACK operations notoriously change the meaning of the matrix
   * entries, we record the current state after the last operation here.
   */
  LAPACKSupport::State state;

  /**
   * Additional property of the matrix which may help to select more
   * efficient ScaLAPACK functions.
   */
  LAPACKSupport::Property property;

  /**
   * A shared pointer to a Utilities::MPI::ProcessGrid object which contains a BLACS context
   * and a MPI communicator, as well as other necessary data structures.
   */
  std::shared_ptr<const Utilities::MPI::ProcessGrid> grid;

  /**
   * Number of rows in the matrix.
   */
  int n_rows;

  /**
   * Number of columns in the matrix.
   */
  int n_columns;

  /**
   * Row block size.
   */
  int row_block_size;

  /**
   * Column block size.
   */
  int column_block_size;

  /**
   * Number of rows in the matrix owned by the current process.
   */
  int n_local_rows;

  /**
   * Number of columns in the matrix owned by the current process.
   */
  int n_local_columns;

  /**
   * ScaLAPACK description vector.
   */
  int descriptor[9];

  /**
   * Workspace array.
   */
  mutable std::vector<NumberType> work;

  /**
   * Integer workspace array.
   */
  mutable std::vector<int> iwork;

  /**
   * A character to define where elements are stored in case
   * ScaLAPACK operations support this.
   */
  const char uplo;

  /**
   * The process row of the process grid over which the first row
   * of the global matrix is distributed.
   */
  const int first_process_row;

  /**
   * The process column of the process grid over which the first column
   * of the global matrix is distributed.
   */
  const int first_process_column;

  /**
   * Global row index that determines where to start a submatrix.
   * Currently this equals unity, as we don't use submatrices.
   */
  const int submatrix_row;

  /**
   * Global column index that determines where to start a submatrix.
   * Currently this equals unity, as we don't use submatrices.
   */
  const int submatrix_column;

  /**
   * Thread mutex.
   */
  mutable Threads::Mutex mutex;
};

// ----------------------- inline functions ----------------------------

#ifndef DOXYGEN

template <typename NumberType>
inline
NumberType
ScaLAPACKMatrix<NumberType>::local_el(const unsigned int loc_row, const unsigned int loc_column) const
{
  return (*this)(loc_row,loc_column);
}



template <typename NumberType>
inline
NumberType &
ScaLAPACKMatrix<NumberType>::local_el(const unsigned int loc_row, const unsigned int loc_column)
{
  return (*this)(loc_row,loc_column);
}


template <typename NumberType>
inline
unsigned int
ScaLAPACKMatrix<NumberType>::m() const
{
  return n_rows;
}



template <typename NumberType>
inline
unsigned int
ScaLAPACKMatrix<NumberType>::n() const
{
  return n_columns;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::local_m() const
{
  return n_local_rows;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::local_n() const
{
  return n_local_columns;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK

#endif
