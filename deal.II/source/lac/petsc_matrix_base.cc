//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/lac/petsc_matrix_base.h>

#ifdef DEAL_II_USE_PETSC

#  include <deal.II/lac/petsc_full_matrix.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>
#  include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#  include <deal.II/lac/petsc_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MatrixIterators
  {
    void
    MatrixBase::const_iterator::Accessor::
    visit_present_row ()
    {
                                       // if we are asked to visit the
                                       // past-the-end line, then simply
                                       // release all our caches and go on
                                       // with life
      if (this->a_row == matrix->m())
        {
          colnum_cache.reset ();
          value_cache.reset ();

          return;
        }

                                       // otherwise first flush PETSc caches
      matrix->compress ();

                                       // get a representation of the present
                                       // row
#if DEAL_II_PETSC_VERSION_LT(2,2,1)
      int          ncols;
      int         *colnums;
      PetscScalar *values;
#else
      PetscInt           ncols;
      const PetscInt    *colnums;
      const PetscScalar *values;
#endif

      int ierr;
      ierr = MatGetRow(*matrix, this->a_row, &ncols, &colnums, &values);
      AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));

                                       // copy it into our caches if the line
                                       // isn't empty. if it is, then we've
                                       // done something wrong, since we
                                       // shouldn't have initialized an
                                       // iterator for an empty line (what
                                       // would it point to?)
      Assert (ncols != 0, ExcInternalError());
      colnum_cache.reset (new std::vector<unsigned int> (colnums,
                                                         colnums+ncols));
      value_cache.reset (new std::vector<PetscScalar> (values, values+ncols));

                                       // and finally restore the matrix
      ierr = MatRestoreRow(*matrix, this->a_row, &ncols, &colnums, &values);
      AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));
    }
  }



  MatrixBase::MatrixBase ()
                  :
                  last_action (LastAction::none)
  {}



  MatrixBase::~MatrixBase ()
  {
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr = MatDestroy (matrix);
#else
    const int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::clear ()
  {
                                     // destroy the matrix...
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    int ierr = MatDestroy (matrix);
#else
    int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));
                                     // ...and replace it by an empty
                                     // sequential matrix
    const int m=0, n=0, n_nonzero_per_row=0;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                           0, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  MatrixBase &
  MatrixBase::operator = (const value_type d)
  {
    Assert (d==value_type(), ExcScalarAssignmentOnlyForZeroValue());

                                     // flush previously cached elements. this
                                     // seems to be necessary since petsc
                                     // 2.2.1, at least for parallel vectors
                                     // (see test bits/petsc_64)
    compress ();

    const int ierr = MatZeroEntries (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  void
  MatrixBase::clear_row (const unsigned int row,
                         const PetscScalar  new_diag_value)
  {
    compress ();

                                     // now set all the entries of this row to
                                     // zero
#if DEAL_II_PETSC_VERSION_LT(2,2,1)
    const int petsc_row      = row;
#else
    const PetscInt petsc_row = row;
#endif

    IS index_set;
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    ISCreateGeneral (get_mpi_communicator(), 1, &petsc_row, &index_set);
#else
    ISCreateGeneral (get_mpi_communicator(), 1, &petsc_row, PETSC_COPY_VALUES, &index_set);
#endif


#if DEAL_II_PETSC_VERSION_LT(2,3,0)
    const int ierr
      = MatZeroRows(matrix, index_set, &new_diag_value);
#else
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value);
#else
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value, PETSC_NULL, PETSC_NULL);
#endif
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    ISDestroy (index_set);
#else
    ISDestroy (&index_set);
#endif

    compress ();
  }



  void
  MatrixBase::clear_rows (const std::vector<unsigned int> &rows,
                          const PetscScalar                new_diag_value)
  {
    compress ();

                                     // now set all the entries of these rows
                                     // to zero
#if DEAL_II_PETSC_VERSION_LT(2,2,1)
    const std::vector<int>      petsc_rows (rows.begin(), rows.end());
#else
    const std::vector<PetscInt> petsc_rows (rows.begin(), rows.end());
#endif

                                     // call the functions. note that we have
                                     // to call them even if #rows is empty,
                                     // since this is a collective operation
    IS index_set;

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    ISCreateGeneral (get_mpi_communicator(), rows.size(),
                     &petsc_rows[0], &index_set);
#else
    ISCreateGeneral (get_mpi_communicator(), rows.size(),
                     &petsc_rows[0], PETSC_COPY_VALUES, &index_set);
#endif

#if DEAL_II_PETSC_VERSION_LT(2,3,0)
    const int ierr
      = MatZeroRows(matrix, index_set, &new_diag_value);
#else
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value);
#else
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value, PETSC_NULL, PETSC_NULL);
#endif
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    ISDestroy (index_set);
#else
    ISDestroy (&index_set);
#endif

    compress ();
  }



  PetscScalar
  MatrixBase::el (const unsigned int i,
                  const unsigned int j) const
  {
#ifdef PETSC_USE_64BIT_INDICES
    PetscInt
#else
    int
#endif
      petsc_i = i, petsc_j = j;
    PetscScalar value;

    const int ierr
      = MatGetValues (matrix, 1, &petsc_i, 1, &petsc_j,
                      &value);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return value;
  }



  PetscScalar
  MatrixBase::diag_element (const unsigned int i) const
  {
    Assert (m() == n(), ExcNotQuadratic());

                                     // this doesn't seem to work any
                                     // different than any other element
    return el(i,i);
  }



  void
  MatrixBase::compress ()
  {
                                     // flush buffers
    int ierr;
    ierr = MatAssemblyBegin (matrix,MAT_FINAL_ASSEMBLY);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = MatAssemblyEnd (matrix,MAT_FINAL_ASSEMBLY);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    last_action = LastAction::none;
  }



  unsigned int
  MatrixBase::m () const
  {
#ifdef PETSC_USE_64BIT_INDICES
    PetscInt
#else
    int
#endif
      n_rows, n_cols;
    int ierr = MatGetSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }



  unsigned int
  MatrixBase::n () const
  {
#ifdef PETSC_USE_64BIT_INDICES
    PetscInt
#else
    int
#endif
      n_rows, n_cols;
    int ierr = MatGetSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_cols;
  }



  unsigned int
  MatrixBase::local_size () const
  {
#ifdef PETSC_USE_64BIT_INDICES
    PetscInt
#else
    int
#endif
      n_rows, n_cols;
    int ierr = MatGetLocalSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }



  std::pair<unsigned int, unsigned int>
  MatrixBase::local_range () const
  {
#ifdef PETSC_USE_64BIT_INDICES
    PetscInt
#else
    int
#endif
      begin, end;
    const int ierr = MatGetOwnershipRange (static_cast<const Mat &>(matrix),
					   &begin, &end);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return std::make_pair (begin, end);
  }



  unsigned int
  MatrixBase::n_nonzero_elements () const
  {
    MatInfo mat_info;
    const int ierr
      = MatGetInfo (matrix, MAT_GLOBAL_SUM, &mat_info);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return static_cast<unsigned int>(mat_info.nz_used);
  }



  unsigned int
  MatrixBase::
  row_length (const unsigned int row) const
  {
//TODO: this function will probably only work if compress() was called on the
//matrix previously. however, we can't do this here, since it would impose
//global communication and one would have to make sure that this function is
//called the same number of times from all processors, something that is
//unreasonable. there should simply be a way in PETSc to query the number of
//entries in a row bypassing the call to compress(), but I can't find one
    Assert (row < m(), ExcInternalError());

                                     // get a representation of the present
                                     // row
#if DEAL_II_PETSC_VERSION_LT(2,2,1)
    int ncols;
    int         *colnums;
    PetscScalar *values;
#else
    PetscInt ncols;
    const PetscInt    *colnums;
    const PetscScalar *values;
#endif

//TODO: this is probably horribly inefficient; we should lobby for a way to
//query this information from PETSc
    int ierr;
    ierr = MatGetRow(*this, row, &ncols, &colnums, &values);
    AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));

                                     // then restore the matrix and return the
                                     // number of columns in this row as
                                     // queried previously
    ierr = MatRestoreRow(*this, row, &ncols, &colnums, &values);
    AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));

    return ncols;
  }


  PetscReal
  MatrixBase::l1_norm () const
  {
    PetscReal result;

    const int ierr
      = MatNorm (matrix, NORM_1, &result);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return result;
  }



  PetscReal
  MatrixBase::linfty_norm () const
  {
    PetscReal result;

    const int ierr
      = MatNorm (matrix, NORM_INFINITY, &result);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return result;
  }



  PetscReal
  MatrixBase::frobenius_norm () const
  {
    PetscReal result;

    const int ierr
      = MatNorm (matrix, NORM_FROBENIUS, &result);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return result;
  }



#if DEAL_II_PETSC_VERSION_GTE(3,1,0)
  PetscReal
  MatrixBase::trace () const
  {
    PetscReal result;

    const int ierr
      = MatGetTrace (matrix, &result);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return result;
  }
#endif



  MatrixBase &
  MatrixBase::operator *= (const PetscScalar a)
  {
#if DEAL_II_PETSC_VERSION_LT(2,3,0)
    const int ierr = MatScale (&a, matrix);
#else
    const int ierr = MatScale (matrix, a);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  MatrixBase &
  MatrixBase::operator /= (const PetscScalar a)
  {
    const PetscScalar factor = 1./a;

#if DEAL_II_PETSC_VERSION_LT(2,3,0)
    const int ierr = MatScale (&factor, matrix);
#else
    const int ierr = MatScale (matrix, factor);
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }


  void
  MatrixBase::vmult (VectorBase       &dst,
                     const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    const int ierr = MatMult (matrix, src, dst);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::Tvmult (VectorBase       &dst,
                      const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    const int ierr = MatMultTranspose (matrix, src, dst);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::vmult_add (VectorBase       &dst,
                         const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    const int ierr = MatMultAdd (matrix, src, dst, dst);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::Tvmult_add (VectorBase       &dst,
                          const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    const int ierr = MatMultTransposeAdd (matrix, src, dst, dst);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  PetscScalar
  MatrixBase::matrix_norm_square (const VectorBase &v) const
  {
    Vector tmp(v.size());
    vmult (tmp, v);
    return tmp*v;
  }



  PetscScalar
  MatrixBase::matrix_scalar_product (const VectorBase &u,
                                     const VectorBase &v) const
  {
    Vector tmp(v.size());
    vmult (tmp, v);
    return u*tmp;
  }



  PetscScalar
  MatrixBase::residual (VectorBase       &dst,
                        const VectorBase &x,
                        const VectorBase &b) const
  {
                                     // avoid the use of a temporary, and
                                     // rather do one negation pass more than
                                     // necessary
    vmult (dst, x);
    dst -= b;
    dst *= -1;

    return dst.l2_norm();
  }



  MatrixBase::operator Mat () const
  {
    return matrix;
  }

  void
  MatrixBase::transpose ()
  {
    int ierr;

#if DEAL_II_PETSC_VERSION_LT(3,0,0)
    ierr = MatTranspose(matrix, PETSC_NULL);
#else
    ierr = MatTranspose(matrix, MAT_REUSE_MATRIX, &matrix);
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
  PetscTruth
#else
  PetscBool
#endif
  MatrixBase::is_symmetric (const double tolerance)
  {
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    PetscTruth 
#else
      PetscBool
#endif
      truth;
                                       // First flush PETSc caches
    compress ();
    MatIsSymmetric (matrix, tolerance, &truth);
    return truth;
  }

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
  PetscTruth
#else
  PetscBool
#endif
  MatrixBase::is_hermitian (const double tolerance)
  {
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    PetscTruth 
#else
      PetscBool
#endif
      truth;

				     // First flush PETSc caches
    compress ();

#if DEAL_II_PETSC_VERSION_LT(3,0,0)
				     // avoid warning about unused variables
    (void) tolerance;

    MatIsHermitian (matrix, &truth);
#else
    MatIsHermitian (matrix, tolerance, &truth);
#endif
    return truth;
  }

  void
  MatrixBase::write_ascii ()
  {
                                       // First flush PETSc caches
    compress ();

                                       // Set options
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
    PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
			  PETSC_VIEWER_ASCII_DEFAULT);
#else
    PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
			  PETSC_VIEWER_DEFAULT);
#endif
                                       // Write to screen
    MatView (matrix,PETSC_VIEWER_STDOUT_WORLD);
  }



  std::size_t
  MatrixBase::memory_consumption() const
  {
    MatInfo info;
    MatGetInfo(matrix, MAT_LOCAL, &info);

    return sizeof(*this) + static_cast<unsigned int>(info.memory);
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
