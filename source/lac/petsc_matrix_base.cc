// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2022 by the deal.II authors
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

#include <deal.II/lac/petsc_matrix_base.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_full_matrix.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>
#  include <deal.II/lac/petsc_vector_base.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MatrixIterators
  {
#  ifndef DOXYGEN
    void
    MatrixBase::const_iterator::Accessor::visit_present_row()
    {
      // if we are asked to visit the past-the-end line (or a line that is not
      // stored on the current processor), then simply release all our caches
      // and go on with life
      if (matrix->in_local_range(this->a_row) == false)
        {
          colnum_cache.reset();
          value_cache.reset();

          return;
        }

      // get a representation of the present row
      PetscInt           ncols;
      const PetscInt *   colnums;
      const PetscScalar *values;

      PetscErrorCode ierr =
        MatGetRow(*matrix, this->a_row, &ncols, &colnums, &values);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      // copy it into our caches if the line
      // isn't empty. if it is, then we've
      // done something wrong, since we
      // shouldn't have initialized an
      // iterator for an empty line (what
      // would it point to?)
      Assert(ncols != 0, ExcInternalError());
      colnum_cache =
        std::make_shared<std::vector<size_type>>(colnums, colnums + ncols);
      value_cache =
        std::make_shared<std::vector<PetscScalar>>(values, values + ncols);

      // and finally restore the matrix
      ierr = MatRestoreRow(*matrix, this->a_row, &ncols, &colnums, &values);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }
#  endif
  } // namespace MatrixIterators



  MatrixBase::MatrixBase()
    : matrix(nullptr)
    , last_action(VectorOperation::unknown)
  {}


  MatrixBase::MatrixBase(const Mat &A)
    : matrix(A)
    , last_action(VectorOperation::unknown)
  {
    const PetscErrorCode ierr =
      PetscObjectReference(reinterpret_cast<PetscObject>(matrix));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  MatrixBase::assign_petsc_matrix(Mat A)
  {
    AssertThrow(last_action == ::dealii::VectorOperation::unknown,
                ExcMessage("Cannot assign a new Mat"));
    PetscErrorCode ierr =
      PetscObjectReference(reinterpret_cast<PetscObject>(A));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    destroy_matrix(matrix);
    matrix = A;
  }

  MatrixBase::~MatrixBase()
  {
    destroy_matrix(matrix);
  }

  void
  MatrixBase::clear()
  {
    // destroy the matrix...
    {
      const PetscErrorCode ierr = destroy_matrix(matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }

    // ...and replace it by an empty
    // sequential matrix
    const int            m = 0, n = 0, n_nonzero_per_row = 0;
    const PetscErrorCode ierr = MatCreateSeqAIJ(
      PETSC_COMM_SELF, m, n, n_nonzero_per_row, nullptr, &matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  MatrixBase &
  MatrixBase::operator=(const value_type d)
  {
    (void)d;
    Assert(d == value_type(), ExcScalarAssignmentOnlyForZeroValue());

    assert_is_compressed();

    const PetscErrorCode ierr = MatZeroEntries(matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  void
  MatrixBase::clear_row(const size_type row, const PetscScalar new_diag_value)
  {
    std::vector<size_type> rows(1, row);
    clear_rows(rows, new_diag_value);
  }



  void
  MatrixBase::clear_rows(const std::vector<size_type> &rows,
                         const PetscScalar             new_diag_value)
  {
    assert_is_compressed();

    // now set all the entries of these rows
    // to zero
    const std::vector<PetscInt> petsc_rows(rows.begin(), rows.end());

    // call the functions. note that we have
    // to call them even if #rows is empty,
    // since this is a collective operation
    IS index_set;

    ISCreateGeneral(get_mpi_communicator(),
                    rows.size(),
                    petsc_rows.data(),
                    PETSC_COPY_VALUES,
                    &index_set);

    const PetscErrorCode ierr =
      MatZeroRowsIS(matrix, index_set, new_diag_value, nullptr, nullptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ISDestroy(&index_set);
  }



  PetscScalar
  MatrixBase::el(const size_type i, const size_type j) const
  {
    PetscInt petsc_i = i, petsc_j = j;

    PetscScalar value;

    const PetscErrorCode ierr =
      MatGetValues(matrix, 1, &petsc_i, 1, &petsc_j, &value);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return value;
  }



  PetscScalar
  MatrixBase::diag_element(const size_type i) const
  {
    Assert(m() == n(), ExcNotQuadratic());

    // this doesn't seem to work any
    // different than any other element
    return el(i, i);
  }



  void
  MatrixBase::compress(const VectorOperation::values operation)
  {
    {
#  ifdef DEBUG
#    ifdef DEAL_II_WITH_MPI
      // Check that all processors agree that last_action is the same (or none!)

      int my_int_last_action = last_action;
      int all_int_last_action;

      const int ierr = MPI_Allreduce(&my_int_last_action,
                                     &all_int_last_action,
                                     1,
                                     MPI_INT,
                                     MPI_BOR,
                                     get_mpi_communicator());
      AssertThrowMPI(ierr);

      AssertThrow(all_int_last_action !=
                    (VectorOperation::add | VectorOperation::insert),
                  ExcMessage("Error: not all processors agree on the last "
                             "VectorOperation before this compress() call."));
#    endif
#  endif
    }

    AssertThrow(
      last_action == VectorOperation::unknown || last_action == operation,
      ExcMessage(
        "Missing compress() or calling with wrong VectorOperation argument."));

    // flush buffers
    PetscErrorCode ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    last_action = VectorOperation::unknown;
  }



  MatrixBase::size_type
  MatrixBase::m() const
  {
    PetscInt n_rows, n_cols;

    const PetscErrorCode ierr = MatGetSize(matrix, &n_rows, &n_cols);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }



  MatrixBase::size_type
  MatrixBase::n() const
  {
    PetscInt n_rows, n_cols;

    const PetscErrorCode ierr = MatGetSize(matrix, &n_rows, &n_cols);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return n_cols;
  }



  MatrixBase::size_type
  MatrixBase::local_size() const
  {
    PetscInt n_rows, n_cols;

    const PetscErrorCode ierr = MatGetLocalSize(matrix, &n_rows, &n_cols);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return n_rows;
  }



  std::pair<MatrixBase::size_type, MatrixBase::size_type>
  MatrixBase::local_range() const
  {
    PetscInt begin, end;

    const PetscErrorCode ierr =
      MatGetOwnershipRange(static_cast<const Mat &>(matrix), &begin, &end);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return std::make_pair(begin, end);
  }



  std::uint64_t
  MatrixBase::n_nonzero_elements() const
  {
    MatInfo              mat_info;
    const PetscErrorCode ierr = MatGetInfo(matrix, MAT_GLOBAL_SUM, &mat_info);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // MatInfo logs quantities as PetscLogDouble. So we need to cast it to match
    // our interface.
    return static_cast<std::uint64_t>(mat_info.nz_used);
  }



  MatrixBase::size_type
  MatrixBase::row_length(const size_type row) const
  {
    // TODO: this function will probably only work if compress() was called on
    // the matrix previously. however, we can't do this here, since it would
    // impose global communication and one would have to make sure that this
    // function is called the same number of times from all processors,
    // something that is unreasonable. there should simply be a way in PETSc to
    // query the number of entries in a row bypassing the call to compress(),
    // but I can't find one
    Assert(row < m(), ExcInternalError());

    // get a representation of the present
    // row
    PetscInt           ncols;
    const PetscInt *   colnums;
    const PetscScalar *values;

    // TODO: this is probably horribly inefficient; we should lobby for a way to
    // query this information from PETSc
    PetscErrorCode ierr = MatGetRow(*this, row, &ncols, &colnums, &values);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then restore the matrix and return the number of columns in this row as
    // queried previously. Starting with PETSc 3.4, MatRestoreRow actually
    // resets the last three arguments to nullptr, to avoid abuse of pointers
    // now dangling. as a consequence, we need to save the size of the array
    // and return the saved value.
    const PetscInt ncols_saved = ncols;
    ierr = MatRestoreRow(*this, row, &ncols, &colnums, &values);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return ncols_saved;
  }


  PetscReal
  MatrixBase::l1_norm() const
  {
    PetscReal result;

    const PetscErrorCode ierr = MatNorm(matrix, NORM_1, &result);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return result;
  }



  PetscReal
  MatrixBase::linfty_norm() const
  {
    PetscReal result;

    const PetscErrorCode ierr = MatNorm(matrix, NORM_INFINITY, &result);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return result;
  }



  PetscReal
  MatrixBase::frobenius_norm() const
  {
    PetscReal result;

    const PetscErrorCode ierr = MatNorm(matrix, NORM_FROBENIUS, &result);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return result;
  }


  PetscScalar
  MatrixBase::matrix_norm_square(const VectorBase &v) const
  {
    VectorBase tmp(v);
    vmult(tmp, v);
    return tmp * v;
  }


  PetscScalar
  MatrixBase::matrix_scalar_product(const VectorBase &u,
                                    const VectorBase &v) const
  {
    VectorBase tmp(u);
    vmult(tmp, v);
    return u * tmp;
  }


  PetscScalar
  MatrixBase::trace() const
  {
    PetscScalar result;

    const PetscErrorCode ierr = MatGetTrace(matrix, &result);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return result;
  }



  MatrixBase &
  MatrixBase::operator*=(const PetscScalar a)
  {
    const PetscErrorCode ierr = MatScale(matrix, a);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  MatrixBase &
  MatrixBase::operator/=(const PetscScalar a)
  {
    const PetscScalar    factor = 1. / a;
    const PetscErrorCode ierr   = MatScale(matrix, factor);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  MatrixBase &
  MatrixBase::add(const PetscScalar factor, const MatrixBase &other)
  {
    const PetscErrorCode ierr =
      MatAXPY(matrix, factor, other, DIFFERENT_NONZERO_PATTERN);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }


  void
  MatrixBase::vmult(VectorBase &dst, const VectorBase &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());

    const PetscErrorCode ierr = MatMult(matrix, src, dst);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::Tvmult(VectorBase &dst, const VectorBase &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());

    const PetscErrorCode ierr = MatMultTranspose(matrix, src, dst);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::vmult_add(VectorBase &dst, const VectorBase &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());

    const PetscErrorCode ierr = MatMultAdd(matrix, src, dst, dst);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  MatrixBase::Tvmult_add(VectorBase &dst, const VectorBase &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());

    const PetscErrorCode ierr = MatMultTransposeAdd(matrix, src, dst, dst);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  namespace internals
  {
    void
    perform_mmult(const MatrixBase &inputleft,
                  const MatrixBase &inputright,
                  MatrixBase &      result,
                  const VectorBase &V,
                  const bool        transpose_left)
    {
      const bool use_vector = (V.size() == inputright.m() ? true : false);
      if (transpose_left == false)
        {
          Assert(inputleft.n() == inputright.m(),
                 ExcDimensionMismatch(inputleft.n(), inputright.m()));
        }
      else
        {
          Assert(inputleft.m() == inputright.m(),
                 ExcDimensionMismatch(inputleft.m(), inputright.m()));
        }

      result.clear();

      PetscErrorCode ierr;

      if (use_vector == false)
        {
          if (transpose_left)
            {
              ierr = MatTransposeMatMult(inputleft,
                                         inputright,
                                         MAT_INITIAL_MATRIX,
                                         PETSC_DEFAULT,
                                         &result.petsc_matrix());
              AssertThrow(ierr == 0, ExcPETScError(ierr));
            }
          else
            {
              ierr = MatMatMult(inputleft,
                                inputright,
                                MAT_INITIAL_MATRIX,
                                PETSC_DEFAULT,
                                &result.petsc_matrix());
              AssertThrow(ierr == 0, ExcPETScError(ierr));
            }
        }
      else
        {
          Mat tmp;
          ierr = MatDuplicate(inputleft, MAT_COPY_VALUES, &tmp);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          if (transpose_left)
            {
#  if DEAL_II_PETSC_VERSION_LT(3, 8, 0)
              ierr = MatTranspose(tmp, MAT_REUSE_MATRIX, &tmp);
#  else
              ierr = MatTranspose(tmp, MAT_INPLACE_MATRIX, &tmp);
#  endif
              AssertThrow(ierr == 0, ExcPETScError(ierr));
            }
          ierr = MatDiagonalScale(tmp, nullptr, V);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          ierr = MatMatMult(tmp,
                            inputright,
                            MAT_INITIAL_MATRIX,
                            PETSC_DEFAULT,
                            &result.petsc_matrix());
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          ierr = PETScWrappers::destroy_matrix(tmp);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
        }
    }
  } // namespace internals

  void
  MatrixBase::mmult(MatrixBase &      C,
                    const MatrixBase &B,
                    const VectorBase &V) const
  {
    internals::perform_mmult(*this, B, C, V, false);
  }

  void
  MatrixBase::Tmmult(MatrixBase &      C,
                     const MatrixBase &B,
                     const VectorBase &V) const
  {
    internals::perform_mmult(*this, B, C, V, true);
  }

  PetscScalar
  MatrixBase::residual(VectorBase &      dst,
                       const VectorBase &x,
                       const VectorBase &b) const
  {
    // avoid the use of a temporary, and
    // rather do one negation pass more than
    // necessary
    vmult(dst, x);
    dst -= b;
    dst *= -1;

    return dst.l2_norm();
  }



  MatrixBase::operator Mat() const
  {
    return matrix;
  }

  Mat &
  MatrixBase::petsc_matrix()
  {
    return matrix;
  }

  void
  MatrixBase::transpose()
  {
#  if DEAL_II_PETSC_VERSION_LT(3, 8, 0)
    const PetscErrorCode ierr = MatTranspose(matrix, MAT_REUSE_MATRIX, &matrix);
#  else
    const PetscErrorCode ierr =
      MatTranspose(matrix, MAT_INPLACE_MATRIX, &matrix);
#  endif
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  PetscBool
  MatrixBase::is_symmetric(const double tolerance)
  {
    PetscBool truth;
    assert_is_compressed();
    const PetscErrorCode ierr = MatIsSymmetric(matrix, tolerance, &truth);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    return truth;
  }

  PetscBool
  MatrixBase::is_hermitian(const double tolerance)
  {
    PetscBool truth;

    assert_is_compressed();
    const PetscErrorCode ierr = MatIsHermitian(matrix, tolerance, &truth);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return truth;
  }

  void
  MatrixBase::write_ascii(const PetscViewerFormat format)
  {
    assert_is_compressed();
    MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(matrix));

    // Set options
    PetscErrorCode ierr =
      PetscViewerSetFormat(PETSC_VIEWER_STDOUT_(comm), format);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // Write to screen
    ierr = MatView(matrix, PETSC_VIEWER_STDOUT_(comm));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  MatrixBase::print(std::ostream &out, const bool /*alternative_output*/) const
  {
    PetscBool has;

    PetscErrorCode ierr = MatHasOperation(matrix, MATOP_GET_ROW, &has);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    Mat vmatrix = matrix;
    if (!has)
      {
        ierr = MatConvert(matrix, MATAIJ, MAT_INITIAL_MATRIX, &vmatrix);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    std::pair<MatrixBase::size_type, MatrixBase::size_type> loc_range =
      local_range();

    PetscInt           ncols;
    const PetscInt *   colnums;
    const PetscScalar *values;

    MatrixBase::size_type row;
    for (row = loc_range.first; row < loc_range.second; ++row)
      {
        ierr = MatGetRow(vmatrix, row, &ncols, &colnums, &values);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        for (PetscInt col = 0; col < ncols; ++col)
          {
            out << "(" << row << "," << colnums[col] << ") " << values[col]
                << std::endl;
          }

        ierr = MatRestoreRow(vmatrix, row, &ncols, &colnums, &values);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    if (vmatrix != matrix)
      {
        ierr = PETScWrappers::destroy_matrix(vmatrix);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    AssertThrow(out.fail() == false, ExcIO());
  }



  std::size_t
  MatrixBase::memory_consumption() const
  {
    MatInfo              info;
    const PetscErrorCode ierr = MatGetInfo(matrix, MAT_LOCAL, &info);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return sizeof(*this) + static_cast<size_type>(info.memory);
  }

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
