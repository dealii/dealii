// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <cerrno>
#include <iostream>
#include <list>
#include <typeinfo>
#include <vector>


DEAL_II_NAMESPACE_OPEN


// include UMFPACK file.
#ifdef DEAL_II_WITH_UMFPACK
#  include <umfpack.h>
#endif



SparseDirectUMFPACK::~SparseDirectUMFPACK ()
{
  clear ();
}


void
SparseDirectUMFPACK::
initialize (const SparsityPattern &)
{}


#ifdef DEAL_II_WITH_UMFPACK

SparseDirectUMFPACK::SparseDirectUMFPACK ()
  :
  symbolic_decomposition (0),
  numeric_decomposition (0),
  control (UMFPACK_CONTROL)
{
  umfpack_dl_defaults (&control[0]);
}



void
SparseDirectUMFPACK::clear ()
{
  // delete objects that haven't been deleted yet
  if (symbolic_decomposition != 0)
    {
      umfpack_dl_free_symbolic (&symbolic_decomposition);
      symbolic_decomposition = 0;
    }

  if (numeric_decomposition != 0)
    {
      umfpack_dl_free_numeric (&numeric_decomposition);
      numeric_decomposition = 0;
    }

  {
    std::vector<long int> tmp;
    tmp.swap (Ap);
  }

  {
    std::vector<long int> tmp;
    tmp.swap (Ai);
  }

  {
    std::vector<double> tmp;
    tmp.swap (Ax);
  }

  umfpack_dl_defaults (&control[0]);
}



template <typename number>
void
SparseDirectUMFPACK::
sort_arrays (const SparseMatrix<number> &matrix)
{
  // do the copying around of entries so that the diagonal entry is in the
  // right place. note that this is easy to detect: since all entries apart
  // from the diagonal entry are sorted, we know that the diagonal entry is
  // in the wrong place if and only if its column index is larger than the
  // column index of the second entry in a row
  //
  // ignore rows with only one or no entry
  for (size_type row=0; row<matrix.m(); ++row)
    {
      // we may have to move some elements that are left of the diagonal
      // but presently after the diagonal entry to the left, whereas the
      // diagonal entry has to move to the right. we could first figure out
      // where to move everything to, but for simplicity we just make a
      // series of swaps instead (this is kind of a single run of
      // bubble-sort, which gives us the desired result since the array is
      // already "almost" sorted)
      //
      // in the first loop, the condition in the while-header also checks
      // that the row has at least two entries and that the diagonal entry
      // is really in the wrong place
      long int cursor = Ap[row];
      while ((cursor < Ap[row+1]-1) &&
             (Ai[cursor] > Ai[cursor+1]))
        {
          std::swap (Ai[cursor], Ai[cursor+1]);
          std::swap (Ax[cursor], Ax[cursor+1]);
          ++cursor;
        }
    }
}



template <typename number>
void
SparseDirectUMFPACK::
sort_arrays (const SparseMatrixEZ<number> &matrix)
{
  //same thing for SparseMatrixEZ
  for (size_type row=0; row<matrix.m(); ++row)
    {
      long int cursor = Ap[row];
      while ((cursor < Ap[row+1]-1) &&
             (Ai[cursor] > Ai[cursor+1]))
        {
          std::swap (Ai[cursor], Ai[cursor+1]);
          std::swap (Ax[cursor], Ax[cursor+1]);
          ++cursor;
        }
    }
}



template <typename number>
void
SparseDirectUMFPACK::
sort_arrays (const BlockSparseMatrix<number> &matrix)
{
  // the case for block matrices is a bit more difficult, since all we know
  // is that *within each block*, the diagonal of that block may come
  // first. however, that means that there may be as many entries per row
  // in the wrong place as there are block columns. we can do the same
  // thing as above, but we have to do it multiple times
  for (size_type row=0; row<matrix.m(); ++row)
    {
      long int cursor = Ap[row];
      for (size_type block=0; block<matrix.n_block_cols(); ++block)
        {
          // find the next out-of-order element
          while ((cursor < Ap[row+1]-1) &&
                 (Ai[cursor] < Ai[cursor+1]))
            ++cursor;

          // if there is none, then just go on
          if (cursor == Ap[row+1]-1)
            break;

          // otherwise swap this entry with successive ones as long as
          // necessary
          long int element = cursor;
          while ((element < Ap[row+1]-1) &&
                 (Ai[element] > Ai[element+1]))
            {
              std::swap (Ai[element], Ai[element+1]);
              std::swap (Ax[element], Ax[element+1]);
              ++element;
            }
        }
    }
}



template <class Matrix>
void
SparseDirectUMFPACK::
factorize (const Matrix &matrix)
{
  Assert (matrix.m() == matrix.n(), ExcNotQuadratic())

  clear ();

  const size_type N = matrix.m();

  // copy over the data from the matrix to the data structures UMFPACK
  // wants. note two things: first, UMFPACK wants compressed column storage
  // whereas we always do compressed row storage; we work around this by,
  // rather than shuffling things around, copy over the data we have, but
  // then call the umfpack_dl_solve function with the UMFPACK_At argument,
  // meaning that we want to solve for the transpose system
  //
  // second: the data we have in the sparse matrices is "almost" right
  // already; UMFPACK wants the entries in each row (i.e. really: column)
  // to be sorted in ascending order. we almost have that, except that we
  // usually store the diagonal first in each row to allow for some
  // optimizations. thus, we have to resort things a little bit, but only
  // within each row
  //
  // final note: if the matrix has entries in the sparsity pattern that are
  // actually occupied by entries that have a zero numerical value, then we
  // keep them anyway. people are supposed to provide accurate sparsity
  // patterns.
  Ap.resize (N+1);
  Ai.resize (matrix.n_nonzero_elements());
  Ax.resize (matrix.n_nonzero_elements());

  // first fill row lengths array
  Ap[0] = 0;
  for (size_type row=1; row<=N; ++row)
    Ap[row] = Ap[row-1] + matrix.get_row_length(row-1);
  Assert (static_cast<size_type>(Ap.back()) == Ai.size(),
          ExcInternalError());

  // then copy over matrix elements. note that for sparse matrices,
  // iterators are sorted so that they traverse each row from start to end
  // before moving on to the next row. however, this isn't true for block
  // matrices, so we have to do a bit of book keeping
  {
    // have an array that for each row points to the first entry not yet
    // written to
    std::vector<long int> row_pointers = Ap;

    // loop over the elements of the matrix row by row, as suggested in the
    // documentation of the sparse matrix iterator class
    for (size_type row = 0; row < matrix.m(); ++row)
      {
        for (typename Matrix::const_iterator p=matrix.begin(row);
             p!=matrix.end(row); ++p)
          {
            // write entry into the first free one for this row
            Ai[row_pointers[row]] = p->column();
            Ax[row_pointers[row]] = p->value();

            // then move pointer ahead
            ++row_pointers[row];
          }
      }

    // at the end, we should have written all rows completely
    for (size_type i=0; i<Ap.size()-1; ++i)
      Assert (row_pointers[i] == Ap[i+1], ExcInternalError());
  }

  // make sure that the elements in each row are sorted. we have to be more
  // careful for block sparse matrices, so ship this task out to a
  // different function
  sort_arrays (matrix);

  int status;
  status = umfpack_dl_symbolic (N, N,
                                &Ap[0], &Ai[0], &Ax[0],
                                &symbolic_decomposition,
                                &control[0], 0);
  AssertThrow (status == UMFPACK_OK,
               ExcUMFPACKError("umfpack_dl_symbolic", status));

  status = umfpack_dl_numeric (&Ap[0], &Ai[0], &Ax[0],
                               symbolic_decomposition,
                               &numeric_decomposition,
                               &control[0], 0);
  AssertThrow (status == UMFPACK_OK,
               ExcUMFPACKError("umfpack_dl_numeric", status));

  umfpack_dl_free_symbolic (&symbolic_decomposition) ;
}



void
SparseDirectUMFPACK::solve (Vector<double> &rhs_and_solution,
                            bool            transpose /*=false*/) const
{
  // make sure that some kind of factorize() call has happened before
  Assert (Ap.size() != 0, ExcNotInitialized());
  Assert (Ai.size() != 0, ExcNotInitialized());
  Assert (Ai.size() == Ax.size(), ExcNotInitialized());

  Vector<double> rhs (rhs_and_solution.size());
  rhs = rhs_and_solution;

  // solve the system. note that since UMFPACK wants compressed column
  // storage instead of the compressed row storage format we use in
  // deal.II's SparsityPattern classes, we solve for UMFPACK's A^T instead

  // Conversely, if we solve for the transpose, we have to use UMFPACK_A
  // instead.
  const int status
    = umfpack_dl_solve (transpose ? UMFPACK_A : UMFPACK_At,
                        &Ap[0], &Ai[0], &Ax[0],
                        rhs_and_solution.begin(), rhs.begin(),
                        numeric_decomposition,
                        &control[0], 0);
  AssertThrow (status == UMFPACK_OK, ExcUMFPACKError("umfpack_dl_solve", status));
}


void
SparseDirectUMFPACK::solve (BlockVector<double> &rhs_and_solution,
                            bool                 transpose /*=false*/) const
{
  // the UMFPACK functions want a contiguous array of elements, so
  // there is no way around copying data around. thus, just copy the
  // data into a regular vector and back
  Vector<double> tmp (rhs_and_solution.size());
  tmp = rhs_and_solution;
  solve (tmp, transpose);
  rhs_and_solution = tmp;
}



template <class Matrix>
void
SparseDirectUMFPACK::solve (const Matrix   &matrix,
                            Vector<double> &rhs_and_solution,
                            bool            transpose /*=false*/)
{
  factorize (matrix);
  solve (rhs_and_solution, transpose);
}


template <class Matrix>
void
SparseDirectUMFPACK::solve (const Matrix        &matrix,
                            BlockVector<double> &rhs_and_solution,
                            bool                 transpose /*=false*/)
{
  factorize (matrix);
  solve (rhs_and_solution, transpose);
}


#else


SparseDirectUMFPACK::SparseDirectUMFPACK ()
  :
  symbolic_decomposition (0),
  numeric_decomposition (0),
  control (0)
{}


void
SparseDirectUMFPACK::clear ()
{}


template <class Matrix>
void SparseDirectUMFPACK::factorize (const Matrix &)
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but configured deal.II without passing the necessary switch to 'cmake'. Please consult the installation instructions in doc/readme.html."));
}


void
SparseDirectUMFPACK::solve (Vector<double> &, bool) const
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but configured deal.II without passing the necessary switch to 'cmake'. Please consult the installation instructions in doc/readme.html."));
}



void
SparseDirectUMFPACK::solve (BlockVector<double> &, bool) const
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but configured deal.II without passing the necessary switch to 'cmake'. Please consult the installation instructions in doc/readme.html."));
}


template <class Matrix>
void
SparseDirectUMFPACK::solve (const Matrix &,
                            Vector<double> &,
                            bool)
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but configured deal.II without passing the necessary switch to 'cmake'. Please consult the installation instructions in doc/readme.html."));
}



template <class Matrix>
void
SparseDirectUMFPACK::solve (const Matrix &,
                            BlockVector<double> &,
                            bool)
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but configured deal.II without passing the necessary switch to 'cmake'. Please consult the installation instructions in doc/readme.html."));
}

#endif


template <class Matrix>
void
SparseDirectUMFPACK::initialize (const Matrix        &M,
                                 const AdditionalData)
{
  this->factorize(M);
}


void
SparseDirectUMFPACK::vmult (
  Vector<double>       &dst,
  const Vector<double> &src) const
{
  dst = src;
  this->solve(dst);
}



void
SparseDirectUMFPACK::vmult (
  BlockVector<double>       &dst,
  const BlockVector<double> &src) const
{
  dst = src;
  this->solve(dst);
}


void
SparseDirectUMFPACK::Tvmult (
  Vector<double> &dst,
  const Vector<double> &src) const
{
  dst = src;
  this->solve(dst, /*transpose=*/ true);
}



void
SparseDirectUMFPACK::Tvmult (
  BlockVector<double>       &dst,
  const BlockVector<double> &src) const
{
  dst = src;
  this->solve(dst, /*transpose=*/ true);
}


void
SparseDirectUMFPACK::vmult_add (
  Vector<double> &,
  const Vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}


void
SparseDirectUMFPACK::Tvmult_add (
  Vector<double> &,
  const Vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}




#ifdef DEAL_II_WITH_MUMPS

SparseDirectMUMPS::SparseDirectMUMPS ()
  :
  initialize_called (false)
{}

SparseDirectMUMPS::~SparseDirectMUMPS ()
{
  id.job = -2;
  dmumps_c (&id);

  // Do some cleaning
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      delete[] a;
      delete[] irn;
      delete[] jcn;
    }
}

template <class Matrix>
void SparseDirectMUMPS::initialize_matrix (const Matrix &matrix)
{
  Assert(matrix.n() == matrix.m(), ExcMessage("Matrix needs to be square."));


  // Check we haven't been here before:
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());

  // Initialize MUMPS instance:
  id.job = -1;
  id.par =  1;
  id.sym =  0;

  // Use MPI_COMM_WORLD as communicator
  id.comm_fortran = -987654;
  dmumps_c (&id);

  // Hand over matrix and right-hand side
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      // Objects denoting a MUMPS data structure:
      //
      // Set number of unknowns
      n   = matrix.n ();

      // number of nonzero elements in matrix
      nz  = matrix.n_actually_nonzero_elements ();

      // representation of the matrix
      a   = new double[nz];

      // matrix indices pointing to the row and column dimensions
      // respectively of the matrix representation above (a): ie. a[k] is
      // the matrix element (irn[k], jcn[k])
      irn = new int[nz];
      jcn = new int[nz];

      size_type index = 0;

      // loop over the elements of the matrix row by row, as suggested in
      // the documentation of the sparse matrix iterator class
      for (size_type row = 0; row < matrix.m(); ++row)
        {
          for (typename Matrix::const_iterator ptr = matrix.begin (row);
               ptr != matrix.end (row); ++ptr)
            if (std::abs (ptr->value ()) > 0.0)
              {
                a[index]   = ptr->value ();
                irn[index] = row + 1;
                jcn[index] = ptr->column () + 1;
                ++index;
              }
        }

      id.n   = n;
      id.nz  = nz;
      id.irn = irn;
      id.jcn = jcn;
      id.a   = a;
    }

  // No outputs
  id.icntl[0] = -1;
  id.icntl[1] = -1;
  id.icntl[2] = -1;
  id.icntl[3] =  0;

  // Exit by setting this flag:
  initialize_called = true;
}

template <class Matrix>
void SparseDirectMUMPS::initialize (const Matrix &matrix,
                                    const Vector<double>       &vector)
{
  // Hand over matrix and right-hand side
  initialize_matrix (matrix);

  copy_rhs_to_mumps(vector);
}

void SparseDirectMUMPS::copy_rhs_to_mumps(const Vector<double> &new_rhs)
{
  Assert(n == new_rhs.size(), ExcMessage("Matrix size and rhs length must be equal."));

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      rhs.resize(n);
      for (size_type i = 0; i < n; ++i)
        rhs[i] = new_rhs (i);

      id.rhs = &rhs[0];
    }
}

void SparseDirectMUMPS::copy_solution (Vector<double> &vector)
{
  Assert(n == vector.size(), ExcMessage("Matrix size and solution vector length must be equal."));
  Assert(n == rhs.size(), ExcMessage("Class not initialized with a rhs vector."));

  // Copy solution into the given vector
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      for (size_type i=0; i<n; ++i)
        vector(i) = rhs[i];

      rhs.resize(0); // remove rhs again
    }
}

template <class Matrix>
void SparseDirectMUMPS::initialize (const Matrix &matrix)
{
  // Initialize MUMPS instance:
  initialize_matrix (matrix);
  // Start factorization
  id.job = 4;
  dmumps_c (&id);
}

void SparseDirectMUMPS::solve (Vector<double> &vector)
{
  //TODO: this could be implemented similar to SparseDirectUMFPACK where
  // the given vector will be used as the RHS. Sadly, there is no easy
  // way to do this without breaking the interface.

  // Check that the solver has been initialized by the routine above:
  Assert (initialize_called == true, ExcNotInitialized());

  // and that the matrix has at least one nonzero element:
  Assert (nz != 0, ExcNotInitialized());

  // Start solver
  id.job = 6; // 6 = analysis, factorization, and solve
  dmumps_c (&id);
  copy_solution (vector);
}

void SparseDirectMUMPS::vmult (Vector<double>       &dst,
                               const Vector<double> &src)
{
  // Check that the solver has been initialized by the routine above:
  Assert (initialize_called == true, ExcNotInitialized());

  // and that the matrix has at least one nonzero element:
  Assert (nz != 0, ExcNotInitialized());

  Assert(n == dst.size(), ExcMessage("Destination vector has the wrong size."));
  Assert(n == src.size(), ExcMessage("Source vector has the wrong size."));

  // Hand over right-hand side
  copy_rhs_to_mumps(src);

  // Start solver
  id.job = 3;
  dmumps_c (&id);
  copy_solution (dst);
}

#endif // DEAL_II_WITH_MUMPS


// explicit instantiations for SparseMatrixUMFPACK
#define InstantiateUMFPACK(MATRIX)                        \
  template                                                \
  void SparseDirectUMFPACK::factorize (const MATRIX &);   \
  template                                                \
  void SparseDirectUMFPACK::solve (const MATRIX   &,      \
                                   Vector<double> &,      \
                                   bool);                 \
  template                                                \
  void SparseDirectUMFPACK::solve (const MATRIX   &,      \
                                   BlockVector<double> &, \
                                   bool);                 \
  template                                                \
  void SparseDirectUMFPACK::initialize (const MATRIX &,   \
                                        const AdditionalData);

InstantiateUMFPACK(SparseMatrix<double>)
InstantiateUMFPACK(SparseMatrix<float>)
InstantiateUMFPACK(SparseMatrixEZ<double>)
InstantiateUMFPACK(SparseMatrixEZ<float>)
InstantiateUMFPACK(BlockSparseMatrix<double>)
InstantiateUMFPACK(BlockSparseMatrix<float>)

// explicit instantiations for SparseDirectMUMPS
#ifdef DEAL_II_WITH_MUMPS
#define InstantiateMUMPS(MATRIX)                          \
  template                                                \
  void SparseDirectMUMPS::initialize (const MATRIX &, const Vector<double> &);\
  template                                                \
  void SparseDirectMUMPS::initialize (const MATRIX &);

InstantiateMUMPS(SparseMatrix<double>)
InstantiateMUMPS(SparseMatrix<float>)
// InstantiateMUMPS(SparseMatrixEZ<double>)
// InstantiateMUMPS(SparseMatrixEZ<float>)
InstantiateMUMPS(BlockSparseMatrix<double>)
InstantiateMUMPS(BlockSparseMatrix<float>)
#endif

DEAL_II_NAMESPACE_CLOSE
