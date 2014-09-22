// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#include <deal.II/lac/petsc_matrix_base.h>

#ifdef DEAL_II_WITH_PETSC

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
      PetscInt           ncols;
      const PetscInt    *colnums;
      const PetscScalar *values;

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
      colnum_cache.reset (new std::vector<size_type> (colnums, colnums+ncols));
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
  MatrixBase::clear_row (const size_type   row,
                         const PetscScalar new_diag_value)
  {
    compress ();

    // now set all the entries of this row to
    // zero
    const PetscInt petsc_row = row;

    IS index_set;
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    ISCreateGeneral (get_mpi_communicator(), 1, &petsc_row, &index_set);
#else
    ISCreateGeneral (get_mpi_communicator(), 1, &petsc_row, PETSC_COPY_VALUES, &index_set);
#endif

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value);
#else
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value, PETSC_NULL, PETSC_NULL);
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
  MatrixBase::clear_rows (const std::vector<size_type> &rows,
                          const PetscScalar             new_diag_value)
  {
    compress ();

    // now set all the entries of these rows
    // to zero
    const std::vector<PetscInt> petsc_rows (rows.begin(), rows.end());

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

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value);
#else
    const int ierr
      = MatZeroRowsIS(matrix, index_set, new_diag_value, PETSC_NULL, PETSC_NULL);
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
  MatrixBase::el (const size_type i,
                  const size_type j) const
  {
    PetscInt petsc_i = i, petsc_j = j;

    PetscScalar value;

    const int ierr
      = MatGetValues (matrix, 1, &petsc_i, 1, &petsc_j,
                      &value);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return value;
  }



  PetscScalar
  MatrixBase::diag_element (const size_type i) const
  {
    Assert (m() == n(), ExcNotQuadratic());

    // this doesn't seem to work any
    // different than any other element
    return el(i,i);
  }



  void
  MatrixBase::compress (::dealii::VectorOperation::values)
  {
    // flush buffers
    int ierr;
    ierr = MatAssemblyBegin (matrix,MAT_FINAL_ASSEMBLY);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = MatAssemblyEnd (matrix,MAT_FINAL_ASSEMBLY);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    last_action = LastAction::none;
  }



  void
  MatrixBase::compress ()
  {
    compress(::dealii::VectorOperation::unknown);
  }



  MatrixBase::size_type
  MatrixBase::m () const
  {
    PetscInt n_rows, n_cols;

    int ierr = MatGetSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }



  MatrixBase::size_type
  MatrixBase::n () const
  {
    PetscInt n_rows, n_cols;

    int ierr = MatGetSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_cols;
  }



  MatrixBase::size_type
  MatrixBase::local_size () const
  {
    PetscInt n_rows, n_cols;

    int ierr = MatGetLocalSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }



  std::pair<MatrixBase::size_type, MatrixBase::size_type>
  MatrixBase::local_range () const
  {
    PetscInt begin, end;

    const int ierr = MatGetOwnershipRange (static_cast<const Mat &>(matrix),
                                           &begin, &end);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return std::make_pair (begin, end);
  }



  MatrixBase::size_type
  MatrixBase::n_nonzero_elements () const
  {
    MatInfo mat_info;
    const int ierr
      = MatGetInfo (matrix, MAT_GLOBAL_SUM, &mat_info);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return static_cast<size_type>(mat_info.nz_used);
  }



  MatrixBase::size_type
  MatrixBase::
  row_length (const size_type row) const
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
    PetscInt ncols;
    const PetscInt    *colnums;
    const PetscScalar *values;

//TODO: this is probably horribly inefficient; we should lobby for a way to
//query this information from PETSc
    int ierr;
    ierr = MatGetRow(*this, row, &ncols, &colnums, &values);
    AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));

    // then restore the matrix and return the number of columns in this row as
    // queried previously. Starting with PETSc 3.4, MatRestoreRow actually
    // resets the last three arguments to zero/NULL, to avoid abuse of pointers
    // now dangling. as a consequence, we need to save the size of the array
    // and return the saved value.
    const PetscInt ncols_saved = ncols;
    ierr = MatRestoreRow(*this, row, &ncols, &colnums, &values);
    AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));

    return ncols_saved;
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


#if DEAL_II_PETSC_VERSION_GTE(3,1,0)
  PetscScalar
  MatrixBase::trace () const
  {
    PetscScalar result;

    const int ierr
      = MatGetTrace (matrix, &result);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return result;
  }
#endif



  MatrixBase &
  MatrixBase::operator *= (const PetscScalar a)
  {
    const int ierr = MatScale (matrix, a);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  MatrixBase &
  MatrixBase::operator /= (const PetscScalar a)
  {
    const PetscScalar factor = 1./a;
    const int ierr = MatScale (matrix, factor);

    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }


  MatrixBase &
  MatrixBase::add (const MatrixBase &other,
                   const PetscScalar factor)
  {
    const int ierr = MatAXPY (matrix, factor,
                              other, DIFFERENT_NONZERO_PATTERN);

    Assert (ierr == 0, ExcPETScError(ierr));

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
    int ierr = MatTranspose(matrix, MAT_REUSE_MATRIX, &matrix);
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
    MatIsHermitian (matrix, tolerance, &truth);

    return truth;
  }

  void
  MatrixBase::write_ascii (const PetscViewerFormat format)
  {
    // First flush PETSc caches
    compress ();

    // Set options
    PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                          format);

    // Write to screen
    MatView (matrix, PETSC_VIEWER_STDOUT_WORLD);
  }

  void
  MatrixBase::print (std::ostream &out,
                     const bool    alternative_output) const
  {
    std::pair<MatrixBase::size_type, MatrixBase::size_type>
    loc_range = local_range();

    PetscInt ncols;
    const PetscInt    *colnums;
    const PetscScalar *values;

    MatrixBase::size_type row;
    for (row = loc_range.first; row < loc_range.second; ++row)
      {
        int ierr;
        ierr = MatGetRow(*this, row, &ncols, &colnums, &values);
        AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));

        for (PetscInt col = 0; col < ncols; ++col)
          {
            out << "(" << row << "," << colnums[col] << ") " << values[col] << std::endl;
          }

        ierr = MatRestoreRow(*this, row, &ncols, &colnums, &values);
        AssertThrow (ierr == 0, MatrixBase::ExcPETScError(ierr));
      }

    AssertThrow (out, ExcIO());
  }



  std::size_t
  MatrixBase::memory_consumption() const
  {
    MatInfo info;
    MatGetInfo(matrix, MAT_LOCAL, &info);

    return sizeof(*this) + static_cast<size_type>(info.memory);
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
