//----------------------------  petsc_sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_sparse_matrix.h  ---------------------------
#ifndef __deal2__petsc_sparse_matrix_h
#define __deal2__petsc_sparse_matrix_h


#include <base/config.h>
#include <base/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <lac/petsc_matrix_base.h>

#include <vector>



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
 * @author Wolfgang Bangerth, 2004
 */
  class SparseMatrix : public MatrixBase
  {
    public:
                                       /**
                                        * Create a sparse matrix of dimensions
                                        * <tt>m</tt> times <tt>n</tt>, with an
                                        * initial guess of
                                        * <tt>n_nonzero_per_row</tt> nonzero
                                        * elements per row. PETSc is able to
                                        * cope with the situation that more
                                        * than this number of elements is
                                        * later allocated for a row, but this
                                        * involves copying data, and is thus
                                        * expensive.
                                        *
                                        * The @p{is_symmetric} flag determines
                                        * whether we should tell PETSc that
                                        * the matrix is going to be symmetric
                                        * (as indicated by the call
                                        * @p{MatSetOption(mat,MAT_SYMMETRIC)}. Note
                                        * that the PETSc documentation states
                                        * that one cannot form an ILU
                                        * decomposition of a matrix for which
                                        * this flag has been set. The default
                                        * value of this flag is @p{false}.
                                        */
      SparseMatrix (const unsigned int m,
                    const unsigned int n,
                    const unsigned int n_nonzero_per_row,
                    const bool         is_symmetric = false);

                                       /**
                                        * Initialize a rectangular matrix with
                                        * <tt>m</tt> rows and <tt>n</tt>
                                        * columns.  The maximal number of
                                        * nonzero entries for each row
                                        * separately is given by the
                                        * <tt>row_lengths</tt> array.
                                        *
                                        * Just as for the other constructors:
                                        * PETSc is able to cope with the
                                        * situation that more than this number
                                        * of elements is later allocated for a
                                        * row, but this involves copying data,
                                        * and is thus expensive.
                                        *
                                        * The @p{is_symmetric} flag determines
                                        * whether we should tell PETSc that
                                        * the matrix is going to be symmetric
                                        * (as indicated by the call
                                        * @p{MatSetOption(mat,MAT_SYMMETRIC)}. Note
                                        * that the PETSc documentation states
                                        * that one cannot form an ILU
                                        * decomposition of a matrix for which
                                        * this flag has been set. The default
                                        * value of this flag is @p{false}.
                                        */
      SparseMatrix (const unsigned int               m,
                    const unsigned int               n,
                    const std::vector<unsigned int> &row_lengths,
                    const bool                       is_symmetric = false);
      
  };
}

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_sparse_matrix.h     ---------------------------*/

#endif
/*----------------------------   petsc_sparse_matrix.h     ---------------------------*/
