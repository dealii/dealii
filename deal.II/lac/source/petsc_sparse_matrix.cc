//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------


#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_vector.h>

#ifdef DEAL_II_USE_PETSC


namespace PETScWrappers
{
  SparseMatrix::SparseMatrix ()
  {
    const int m=0, n=0, n_nonzero_per_row=0;
    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                        0, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  SparseMatrix::SparseMatrix (const unsigned int m,
                              const unsigned int n,
                              const unsigned int n_nonzero_per_row,
                              const bool         is_symmetric)
  {
    do_reinit (m, n, n_nonzero_per_row, is_symmetric);
  }



  SparseMatrix::SparseMatrix (const unsigned int m,
                              const unsigned int n,
                              const std::vector<unsigned int> &row_lengths,
                              const bool         is_symmetric)
  {
    do_reinit (m, n, row_lengths, is_symmetric);
  }



  SparseMatrix &
  SparseMatrix::operator = (const double d)
  {
    MatrixBase::operator = (d);
    return *this;
  }



  void
  SparseMatrix::reinit (const unsigned int m,
                        const unsigned int n,
                        const unsigned int n_nonzero_per_row,
                        const bool         is_symmetric)
  {
                                     // get rid of old matrix and generate a
                                     // new one
    const int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));    
    
    do_reinit (m, n, n_nonzero_per_row, is_symmetric);
  }



  void
  SparseMatrix::reinit (const unsigned int m,
                        const unsigned int n,
                        const std::vector<unsigned int> &row_lengths,
                        const bool         is_symmetric)
  {
                                     // get rid of old matrix and generate a
                                     // new one
    const int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));    

    do_reinit (m, n, row_lengths, is_symmetric);
  }  



  void
  SparseMatrix::do_reinit (const unsigned int m,
                           const unsigned int n,
                           const unsigned int n_nonzero_per_row,
                           const bool         is_symmetric)
  {
                                     // use the call sequence indicating only
                                     // a maximal number of elements per row
                                     // for all rows globally
    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                        0, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // set symmetric flag, if so requested
    if (is_symmetric == true)
      {
        const int ierr
          = MatSetOption (matrix, MAT_SYMMETRIC);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }
  }



  void
  SparseMatrix::do_reinit (const unsigned int m,
                           const unsigned int n,
                           const std::vector<unsigned int> &row_lengths,
                           const bool         is_symmetric)
  {
    Assert (row_lengths.size() == m,
            ExcDimensionMismatch (row_lengths.size(), m));
    
                                     // use the call sequence indicating a
                                     // maximal number of elements for each
                                     // row individually. annoyingly, we
                                     // always use unsigned ints for cases
                                     // like this, while PETSc wants to see
                                     // signed integers. so we have to
                                     // convert, unless we want to play dirty
                                     // tricks with conversions of pointers
    const std::vector<signed int> int_row_lengths (row_lengths.begin(),
                                                   row_lengths.end());
    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0,
                        &int_row_lengths[0], &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // set symmetric flag, if so requested
    if (is_symmetric == true)
      {
        const int ierr
          = MatSetOption (matrix, MAT_SYMMETRIC);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }    
  }
}


#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
