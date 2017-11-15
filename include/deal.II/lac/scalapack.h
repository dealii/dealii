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

// FIXME: #ifdef DEAL_II_WITH_SCALAPACK (<--- Lapack+MPI+Scalapack)

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <mpi.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A wrapper class around ScaLAPACK parallel dense linear algebra.
 *
 * ScaLAPACK assumes that matrices to be distributed according to the
 * block-cyclic decomposition scheme. A $M$ by $N$ matrix is first decomposed
 * into $MB$ by $NB$ blocks which are then uniformly distributed across the 2D process grid
 * $p*q <= Np$.
 *
 * For example, global real symmetric matrix of order 9 is stored in upper storage mode with block sizes 4 × 4:
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
 *
 * Currently, only symmetric real matrices are supported.
 *
 * Here is a strong scaling example of ScaLAPACKMatrix::invert()
 * on up to 5 nodes each composed of two Intel Xeon 2660v2 IvyBridge sockets
 * 2.20GHz, 10 cores/socket. Calculate are performed on square processor
 * grids 1x1, 2x2, 3x3, 4x4, 5x5, 6x6, 7x7, 8x8, 9x9, 10x10.
 *
 * @image html scalapack_invert.png
 *
 * @ingroup Matrix1
 * @author Denis Davydov, 2017
 */
template <typename NumberType>
class ScaLAPACKMatrix : protected TransposeTable<NumberType>
{
public:

  /**
    * Declare type for container size.
    */
  typedef unsigned int size_type;

  /**
   * Constructor for a rectangular matrix with @p rows and @p columns, distributed
   * in a given @p mpi_communicator .
   *
   * The choice of the block size @p block_size is a compromize between a large
   * enough sizes for efficient local BLAS and small enough get
   * good parallel load balance.
   */
  ScaLAPACKMatrix(const size_type rows, const size_type columns,
                  MPI_Comm mpi_communicator,
                  const unsigned int block_size_row = 32, const unsigned int block_size_column = 32,
                  const LAPACKSupport::Property property = LAPACKSupport::Property::general);

  /**
   * Constructor for a square matrix of size @p size, distributed
   * in a given @p mpi_communicator .
   *
   * The choice of the block size @p block_size is a compromize between a large
   * enough sizes for efficient local BLAS and small enough get
   * good parallel load balance.
   */
  ScaLAPACKMatrix(const size_type size,
                  MPI_Comm mpi_communicator,
                  const unsigned int block_size = 32);

  /**
   * Destructor
   */
  virtual ~ScaLAPACKMatrix();

  /**
   * Assignment operator from a regular FullMatrix.
   */
  ScaLAPACKMatrix<NumberType> &
  operator = (const FullMatrix<NumberType> &);

  /**
   * Copy the content of the distributed matrix into @p matrix.
   */
  void copy_to (FullMatrix<NumberType> &matrix) const;

  /**
   * Compute the Cholesky factorization of the matrix using ScaLAPACK function pXpotrf.
   */
  void compute_cholesky_factorization ();

  /**
   * Invert the matrix by first computing Cholesky factorization and then
   * building the actual inverse using pXpotri.
   */
  void invert();

  /**
   * Compute all eigenvalues of real symmetric matrix using pdsyev
   */
  void eigenvalues_symmetric (std::vector<NumberType> &eigenvalues);

  /**
   * Compute all eigenpairs of real symmetric matrix using pdsyev
   */
  void eigenpairs_symmetric (std::vector<NumberType> &eigenvalues);

  /**
   * Return the number of rows in processes grid.
   */
  int get_process_grid_rows() const;

  /**
   * Return the number of columns in processes grid.
   */
  int get_process_grid_columns() const;

  /**
   * Estimate the the condition number of a SPD matrix in $l_1$-norm.
   * The matrix has to be in Cholesky state. The reciprocal of the
   * condition number is returned in order to avoid the possibility of
   * overflow when the condition number is very large.
   *
   * @p a_norm shall contain $l_1$-norm of the matrix prior to calling
   * Cholesky factorization.
   *
   * @note an alternative is to compute the inverse of the matrix
   * explicitly and manually constructor $k_1 = ||A||_1 ||A^{-1}||_1$.
   */
  NumberType reciprocal_condition_number(const NumberType a_norm) const;

  /**
   * Compute $l_1$-norm of the matrix.
   * @return
   */
  NumberType l1_norm() const;

  /**
   * Compute $l_{\infty}$ norm of the matrix.
   */
  NumberType linfty_norm() const;

  /**
   * Compute the Frobenius norm of the matrix.
   */
  NumberType frobenius_norm() const;

  /**
   * Number of rows of the matrix $mxn$.
   */
  int m() const;

  /**
   * Number of columns of the matrix $m x n$.
   */
  int n() const;

  /**
   * Number of local rows on this MPI processes.
   */
  int local_m() const;

  /**
   * Number of local columns on this MPI process.
   */
  int local_n() const;

  /**
   * Return global row number for the given local row @p loc_row .
   */
  int global_row(const int loc_row) const;

  /**
   * Return global column number for the given local column @p loc_column.
   */
  int global_column(const int loc_column) const;

  /**
   * Read access to local element.
   */
  NumberType local_el(const int loc_row, const int loc_column) const;

  /**
   * Write access to local element.
   */
  NumberType &local_el(const int loc_row, const int loc_column);

private:

  /**
   * Calculate norm of a distributed dense matrix using ScaLAPACK's
   * internal function.
   */
  NumberType norm(const char type) const;

  /**
   * Send @p value from process with rank zero to processes which
   * are not in the process grid.
   */
  void send_to_inactive(NumberType &value) const;

  /**
   * Since ScaLAPACK operations notoriously change the meaning of the matrix
   * entries, we record the current state after the last operation here.
   */
  LAPACKSupport::State state;

  /**
   * Additional properties of the matrix which may help to select more
   * efficient ScaLAPACK functions.
   */
  LAPACKSupport::Property properties;

  /**
   * MPI communicator with all processes
   */
  MPI_Comm mpi_communicator;

  /**
   * MPI communicator with inactive processes and the root
   */
  MPI_Comm mpi_communicator_inactive_with_root;

  /**
   * Number of rows in the matrix
   */
  int n_rows;

  /**
   * Number of columns in the matrix
   */
  int n_columns;

  /**
   * Row block size
   */
  int row_block_size;

  /**
   * Column block size
   */
  int column_block_size;

  /**
   * BLACS context
   */
  int blacs_context;

  /**
   * ID of this MPI process
   */
  int this_mpi_process;

  /**
   * Total number of MPI processes
   */
  int n_mpi_processes;

  /**
   * Number of rows in processes grid
   */
  int n_process_rows;

  /**
   * Number of columns in processes grid
   */
  int n_process_columns;

  /**
   * Row of this process in the grid
   */
  int this_process_row;

  /**
   * Column of this process in the grid
   */
  int this_process_column;

  /**
   * Number of rows in the matrix owned by the current process
   */
  int n_local_rows;

  /**
   * Number of columns in the matrix owned by the current process
   */
  int n_local_columns;

  /**
   * Scalapack description vector.
   */
  int descriptor[9];

  /**
   * Workspace array
   */
  mutable std::vector<NumberType> work;

  /**
   * Integer workspace array
   */
  mutable std::vector<int> iwork;

  /**
   * A flag which is true for processes within the 2D process grid
   */
  bool active;

  /**
   * A character to define where elements are stored in case
   * ScaLAPACK operation support this.
   *
   * FIXME: switch to enum UpperOrLower: LOWER UPPER
   */
  const char uplo;

  /**
   * The process row of the process grid over which the first row
   * of the global matrix is distributed.
   */
  const int first_process_row;

  /**
   * The process column of the process grid over which the first column
   * of the global matrix is distributed
   */
  const int first_process_column;

  /**
   * Global indices where to start a submatrix. Equals to unity as we
   * don't use submatrices.
   */
  const int submatrix_row;

  /**
   * Global indices where to start a submatrix. Equals to unity as we
   * don't use submatrices.
   */
  const int submatrix_column;

};

// ----------------------- inline functions ----------------------------

template <typename NumberType>
inline
NumberType
ScaLAPACKMatrix<NumberType>::local_el(const int loc_row, const int loc_column) const
{
  return (*this)(loc_row,loc_column);
}



template <typename NumberType>
inline
NumberType &
ScaLAPACKMatrix<NumberType>::local_el(const int loc_row, const int loc_column)
{
  return (*this)(loc_row,loc_column);
}

DEAL_II_NAMESPACE_CLOSE

#endif
