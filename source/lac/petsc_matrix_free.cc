// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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


#include <deal.II/lac/petsc_matrix_free.h>

#ifdef DEAL_II_WITH_PETSC

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  MatrixFree::MatrixFree ()
    : communicator (PETSC_COMM_SELF)
  {
    const int m=0;
    do_reinit (m, m, m, m);
  }



  MatrixFree::MatrixFree (const MPI_Comm     &communicator,
                          const unsigned int  m,
                          const unsigned int  n,
                          const unsigned int  local_rows,
                          const unsigned int  local_columns)
    : communicator (communicator)
  {
    do_reinit (m, n, local_rows, local_columns);
  }



  MatrixFree::MatrixFree (const MPI_Comm     &communicator,
                          const unsigned int  m,
                          const unsigned int  n,
                          const std::vector<unsigned int> &local_rows_per_process,
                          const std::vector<unsigned int> &local_columns_per_process,
                          const unsigned int  this_process)
    : communicator (communicator)
  {
    Assert (local_rows_per_process.size() == local_columns_per_process.size(),
            ExcDimensionMismatch (local_rows_per_process.size(),
                                  local_columns_per_process.size()));
    Assert (this_process < local_rows_per_process.size(),
            ExcInternalError());

    do_reinit (m, n,
               local_rows_per_process[this_process],
               local_columns_per_process[this_process]);
  }



  MatrixFree::MatrixFree (const unsigned int  m,
                          const unsigned int  n,
                          const unsigned int  local_rows,
                          const unsigned int  local_columns)
    : communicator (MPI_COMM_WORLD)
  {
    do_reinit (m, n, local_rows, local_columns);
  }



  MatrixFree::MatrixFree (const unsigned int  m,
                          const unsigned int  n,
                          const std::vector<unsigned int> &local_rows_per_process,
                          const std::vector<unsigned int> &local_columns_per_process,
                          const unsigned int  this_process)
    : communicator (MPI_COMM_WORLD)
  {
    Assert (local_rows_per_process.size() == local_columns_per_process.size(),
            ExcDimensionMismatch (local_rows_per_process.size(),
                                  local_columns_per_process.size()));
    Assert (this_process < local_rows_per_process.size(),
            ExcInternalError());

    do_reinit (m, n,
               local_rows_per_process[this_process],
               local_columns_per_process[this_process]);
  }



  void MatrixFree::reinit (const MPI_Comm     &communicator,
                           const unsigned int  m,
                           const unsigned int  n,
                           const unsigned int  local_rows,
                           const unsigned int  local_columns)
  {
    this->communicator = communicator;

    // destroy the matrix and
    // generate a new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    int ierr = MatDestroy (matrix);
#else
    int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (m, n, local_rows, local_columns);
  }



  void MatrixFree::reinit (const MPI_Comm     &communicator,
                           const unsigned int  m,
                           const unsigned int  n,
                           const std::vector<unsigned int> &local_rows_per_process,
                           const std::vector<unsigned int> &local_columns_per_process,
                           const unsigned int  this_process)
  {
    Assert (local_rows_per_process.size() == local_columns_per_process.size(),
            ExcDimensionMismatch (local_rows_per_process.size(),
                                  local_columns_per_process.size()));
    Assert (this_process < local_rows_per_process.size(),
            ExcInternalError());

    this->communicator = communicator;

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    int ierr = MatDestroy (matrix);
#else
    int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (m, n,
               local_rows_per_process[this_process],
               local_columns_per_process[this_process]);
  }



  void MatrixFree::reinit (const unsigned int  m,
                           const unsigned int  n,
                           const unsigned int  local_rows,
                           const unsigned int  local_columns)
  {
    reinit (MPI_COMM_WORLD, m, n, local_rows, local_columns);
  }



  void MatrixFree::reinit (const unsigned int  m,
                           const unsigned int  n,
                           const std::vector<unsigned int> &local_rows_per_process,
                           const std::vector<unsigned int> &local_columns_per_process,
                           const unsigned int  this_process)
  {
    reinit (MPI_COMM_WORLD, m, n, local_rows_per_process, local_columns_per_process, this_process);
  }



  void MatrixFree::clear ()
  {
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    int ierr = MatDestroy (matrix);
#else
    int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    const int m=0;
    do_reinit (m, m, m, m);
  }



  void MatrixFree::vmult (Vec  &dst, const Vec  &src) const
  {

//TODO: Translate the given PETSc Vec* vector into a deal.II
// vector so we can call the vmult function with the usual
// interface; then convert back. This could be much more
// efficient, if the PETScWrappers::*::Vector classes
// had a way to simply generate such a vector object from
// a given PETSc Vec* object without allocating new memory
// and without taking ownership of the Vec*

    VectorBase  *x = 0;
    VectorBase  *y = 0;
    // because we do not know,
    // if dst and src are sequential
    // or distributed vectors,
    // we ask for the vector-type
    // and reinit x and y with
    // dealii::PETScWrappers::*::Vector:
    const char  *vec_type;
    int ierr = VecGetType (src, &vec_type);

    PetscInt  local_size;
    ierr = VecGetLocalSize (src, &local_size);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    if (strcmp(vec_type,"mpi") == 0)
      {
        PetscInt  size;
        ierr = VecGetSize (src, &size);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        x = new PETScWrappers::MPI::Vector (this->get_mpi_communicator (), size, local_size);
        y = new PETScWrappers::MPI::Vector (this->get_mpi_communicator (), size, local_size);
      }
    else if (strcmp(vec_type,"seq") == 0)
      {
        x = new PETScWrappers::Vector (local_size);
        y = new PETScWrappers::Vector (local_size);
      }
    else
      AssertThrow (false, ExcMessage("PETScWrappers::MPI::MatrixFree::do_matrix_vector_action: "
                                     "This only works for Petsc Vec Type = VECMPI | VECSEQ"));

    // copy src to x
    x->equ(1., PETScWrappers::VectorBase(src));
    // and call vmult(x,y) which must
    // be reimplemented in derived classes
    vmult (*y, *x);

    // copy the result back to dst
    ierr = VecCopy (static_cast<const Vec &>(*y), dst);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    delete (x);
    delete (y);
  }



  int MatrixFree::matrix_free_mult (Mat  A, Vec  src, Vec  dst)
  {
    // create a pointer to this MatrixFree
    // object and link the given matrix A
    // to the matrix-vector multiplication
    // of this MatrixFree object,
    void  *this_object;
    int ierr = MatShellGetContext (A, &this_object);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    // call vmult of this object:
    reinterpret_cast<MatrixFree *>(this_object)->vmult (dst, src);

    return (0);
  }



  void MatrixFree::do_reinit (const unsigned int  m,
                              const unsigned int  n,
                              const unsigned int  local_rows,
                              const unsigned int  local_columns)
  {
    Assert (local_rows <= m, ExcDimensionMismatch (local_rows, m));
    Assert (local_columns <= n, ExcDimensionMismatch (local_columns, n));

    int ierr;
    // create a PETSc MatShell matrix-type
    // object of dimension m x n and local size
    // local_rows x local_columns
    ierr = MatCreateShell(communicator, local_rows, local_columns, m, n, (void *)this, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    // register the MatrixFree::matrix_free_mult function
    // as the matrix multiplication used by this matrix
    ierr = MatShellSetOperation (matrix, MATOP_MULT,
                                 (void( *)(void))&dealii::PETScWrappers::MatrixFree::matrix_free_mult);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = MatSetFromOptions (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
