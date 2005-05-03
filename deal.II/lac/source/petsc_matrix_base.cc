//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_matrix_base.h>
#include <lac/petsc_full_matrix.h>
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_vector.h>

#ifdef DEAL_II_USE_PETSC


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
      int ncols;

#if (PETSC_VERSION_MAJOR <= 2) && \
    ((PETSC_VERSION_MINOR < 2) ||  \
     ((PETSC_VERSION_MINOR == 2) && (PETSC_VERSION_SUBMINOR == 0)))
      int         *colnums;
      PetscScalar *values;
#else
      const int         *colnums;
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
    const int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));    
  }



  void
  MatrixBase::clear ()
  {
                                     // destroy the matrix...
    int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));    
                                     // ...and replace it by an empty
                                     // sequential matrix
    const int m=0, n=0, n_nonzero_per_row=0;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                           0, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  


  MatrixBase &
  MatrixBase::operator = (const double d)
  {
    Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());    

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
  MatrixBase::set (const unsigned int i,
                   const unsigned int j,
                   const PetscScalar value)
  {
    if (last_action != LastAction::insert)
      {
        int ierr;
        ierr = MatAssemblyBegin(matrix,MAT_FLUSH_ASSEMBLY);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        ierr = MatAssemblyEnd(matrix,MAT_FLUSH_ASSEMBLY);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }
    
    const signed int petsc_i = i;
    const signed int petsc_j = j;
          
    const int ierr
      = MatSetValues (matrix, 1, &petsc_i, 1, &petsc_j,
                      &value, INSERT_VALUES);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    last_action = LastAction::insert;
  }



  void
  MatrixBase::add (const unsigned int i,
                   const unsigned int j,
                   const PetscScalar value)
  {
    if (last_action != LastAction::add)
      {
        int ierr;
        ierr = MatAssemblyBegin(matrix,MAT_FLUSH_ASSEMBLY);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        ierr = MatAssemblyEnd(matrix,MAT_FLUSH_ASSEMBLY);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        last_action = LastAction::add;
      }
    
                                     // we have to do above actions in any
                                     // case to be consistent with the MPI
                                     // communication model (see the
                                     // comments in the documentation of
                                     // PETScWrappers::MPI::Vector), but we
                                     // can save some work if the addend is
                                     // zero
    if (value == 0)
      return;

    const signed int petsc_i = i;
    const signed int petsc_j = j;
          
    const int ierr
      = MatSetValues (matrix, 1, &petsc_i, 1, &petsc_j,
                      &value, ADD_VALUES);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


  PetscScalar
  MatrixBase::el (const unsigned int i,
                  const unsigned int j) const
  {
    const signed int petsc_i = i;
    const signed int petsc_j = j;
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

                                     // try to compress the representation
    ierr = MatCompress (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  unsigned int
  MatrixBase::m () const
  {
    int n_rows, n_cols;
    int ierr = MatGetSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }
  


  unsigned int
  MatrixBase::n () const
  {
    int n_rows, n_cols;
    int ierr = MatGetSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_cols;
  }



  unsigned int
  MatrixBase::local_size () const
  {
    int n_rows, n_cols;
    int ierr = MatGetLocalSize (matrix, &n_rows, &n_cols);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }


  
  std::pair<unsigned int, unsigned int>
  MatrixBase::local_range () const
  {
    int begin, end;
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
    int ncols;

#if (PETSC_VERSION_MAJOR <= 2) && \
    ((PETSC_VERSION_MINOR < 2) ||  \
     ((PETSC_VERSION_MINOR == 2) && (PETSC_VERSION_SUBMINOR == 0)))
    int         *colnums;
    PetscScalar *values;
#else
    const int         *colnums;
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



  MatrixBase &
  MatrixBase::operator *= (const PetscScalar a)
  {
    const int ierr
      = MatScale (&a, matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  MatrixBase &
  MatrixBase::operator /= (const PetscScalar a)
  {
    const PetscScalar factor = 1./a;
    
    const int ierr
      = MatScale (&factor, matrix);
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
  

  
  MatrixBase::operator const Mat () const
  {
    return matrix;
  }  
}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
