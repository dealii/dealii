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


#ifndef dealii_superlu_mt_h
#define dealii_superlu_mt_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "slu_mt_ddefs.h"
//#include "slu_mt_sdefs.h"

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Wrapper for SuperLU-MT direct linear solver.
 * \todo Complete the docu here.
 * *
 * @ingroup Solvers Preconditioners
 *
 * @author Paras Kumar, 2020
 */
class SuperLU_MT : public Subscriptor
{
public:
  /**
   * Declare type for container sizes.
   */
  using size_type = types::global_dof_index;


  // \todo remove this issue below
  /**
   * in order to avoid the error with
   * <c> invalid conversion from ‘unsigned int*’ to ‘int_t* {aka int*}’ </c>
   * we create a temporary type. Ideally this should be same as \c size_type
   */
  using my_int = int_t;

  /**
   * Declare types for matrix and vectorm for extension to block sparse matrices
   * (to be replaced by instantiation later)
   */
  using number_type = double;
  using matrix_type = SparseMatrix<number_type>;
  using vector_type = Vector<number_type>;

  /**
   * To contain additional data such as number of process for multi-threading.
   */
  struct AdditionalData //: public superlumt_options_t
  {
  public:
    AdditionalData()
      : column_ordering(0)
    {
      options_.nprocs    = MultithreadInfo::n_threads();
      options_.fact      = DOFACT;
      options_.PrintStat = NO;
    }

    explicit AdditionalData(const unsigned int &num_threads,
                            const unsigned int &factorize,
                            const size_type &   col_order,
                            const bool          printSolverStats)
      : column_ordering(col_order)
    {
      options_.nprocs = num_threads;
      options_.fact   = static_cast<fact_t>(factorize);
      if (printSolverStats)
        options_.PrintStat = YES;
      else
        options_.PrintStat = NO;
    }
    /**
     * Object to store additional options for the solver
     */
    superlumt_options_t options_;
    /**
     * Parameter to define the type of column ordering for LU factorization.
     * See \c permc_spec argument of \c get_perm_c routine for details.
     */
    size_type column_ordering = 0;
  };

  /**
   * Construcotr.
   */
  SuperLU_MT(const AdditionalData &data = AdditionalData());

  /**
   * Desctructor
   */
  ~SuperLU_MT() override;

  /**
   * Function to intialize the datastructures required by the SuperLU-MT
   * solver.
   */
  void
  initialize(const matrix_type &A, const vector_type &b);

  void
  solve(vector_type &x);

private:
  /**
   * Function to clear the memory allocated.
   */
  void
  clear();

  /**
   * Function to generate the \ac SuperMatrix object as required by the
   * SuperLU-MT solver.
   */
  void
  create_supermatrices();
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
   * \c SuperMatrix objects to store the data as required by SuperLU.
   * SuperLU is actually meant to solve several systems of linear
   * equations of the form <c> AX = B <\c>, where \c X and \c B are dense
   * matrices comprising of \c nrhs different right hand sides adn solutions
   * respectively. Thus, we need to store the right hand side (which will be
   * overwritten by the solution) as a \c Supermatrix with one column.
   */
  SuperMatrix A;
  SuperMatrix b;

  /**
   * Object to store additional input arguments for the solver.
   */
  AdditionalData additional_data;

  /**
   * Arrays to store the sparse matrix in the compressed row (SLU_NR)
   * format of SuperLU.
   * \c colind: array of column indices of the nonzero values
   * \c rowptr: array containing beginning of columns in \c nzval_NR and
   * \c colind. Note: Zero-based indexing is used; \c rowptr has
   * n_rows+1 entries, the last one pointing beyond the last row,
   * so that rowptr[n_rows] = nnz.
   * \c nzval_NR: array of nonzero values, packed by row
   */
  std::vector<my_int>      colind;
  std::vector<my_int>      rowptr;
  std::vector<number_type> nzval;

  /**
   * Since, the right hand side gets overwritten by solution we need to
   * make a copy of it.
   */
  std::vector<number_type> rhs_vals;
};


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_superlu_mt_h
