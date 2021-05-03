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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <cerrno>
#include <complex>
#include <iostream>
#include <list>
#include <typeinfo>
#include <vector>


DEAL_II_NAMESPACE_OPEN


// include UMFPACK file.
#ifdef DEAL_II_WITH_UMFPACK
#  include <umfpack.h>
#endif



SparseDirectUMFPACK::~SparseDirectUMFPACK()
{
  clear();
}


void
SparseDirectUMFPACK::initialize(const SparsityPattern &)
{}


#ifdef DEAL_II_WITH_UMFPACK

SparseDirectUMFPACK::SparseDirectUMFPACK()
  : n_rows(0)
  , n_cols(0)
  , symbolic_decomposition(nullptr)
  , numeric_decomposition(nullptr)
  , control(UMFPACK_CONTROL)
{
  umfpack_dl_defaults(control.data());
}



void
SparseDirectUMFPACK::clear()
{
  // delete objects that haven't been deleted yet
  if (symbolic_decomposition != nullptr)
    {
      umfpack_dl_free_symbolic(&symbolic_decomposition);
      symbolic_decomposition = nullptr;
    }

  if (numeric_decomposition != nullptr)
    {
      umfpack_dl_free_numeric(&numeric_decomposition);
      numeric_decomposition = nullptr;
    }

  {
    std::vector<long int> tmp;
    tmp.swap(Ap);
  }

  {
    std::vector<long int> tmp;
    tmp.swap(Ai);
  }

  {
    std::vector<double> tmp;
    tmp.swap(Ax);
  }

  {
    std::vector<double> tmp;
    tmp.swap(Az);
  }

  umfpack_dl_defaults(control.data());
}



template <typename number>
void
SparseDirectUMFPACK::sort_arrays(const SparseMatrix<number> &matrix)
{
  // do the copying around of entries so that the diagonal entry is in the
  // right place. note that this is easy to detect: since all entries apart
  // from the diagonal entry are sorted, we know that the diagonal entry is
  // in the wrong place if and only if its column index is larger than the
  // column index of the second entry in a row
  //
  // ignore rows with only one or no entry
  for (size_type row = 0; row < matrix.m(); ++row)
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
      while ((cursor < Ap[row + 1] - 1) && (Ai[cursor] > Ai[cursor + 1]))
        {
          std::swap(Ai[cursor], Ai[cursor + 1]);

          std::swap(Ax[cursor], Ax[cursor + 1]);
          if (numbers::NumberTraits<number>::is_complex == true)
            std::swap(Az[cursor], Az[cursor + 1]);

          ++cursor;
        }
    }
}



template <typename number>
void
SparseDirectUMFPACK::sort_arrays(const SparseMatrixEZ<number> &matrix)
{
  // same thing for SparseMatrixEZ
  for (size_type row = 0; row < matrix.m(); ++row)
    {
      long int cursor = Ap[row];
      while ((cursor < Ap[row + 1] - 1) && (Ai[cursor] > Ai[cursor + 1]))
        {
          std::swap(Ai[cursor], Ai[cursor + 1]);

          std::swap(Ax[cursor], Ax[cursor + 1]);
          if (numbers::NumberTraits<number>::is_complex == true)
            std::swap(Az[cursor], Az[cursor + 1]);

          ++cursor;
        }
    }
}



template <typename number>
void
SparseDirectUMFPACK::sort_arrays(const BlockSparseMatrix<number> &matrix)
{
  // the case for block matrices is a bit more difficult, since all we know
  // is that *within each block*, the diagonal of that block may come
  // first. however, that means that there may be as many entries per row
  // in the wrong place as there are block columns. we can do the same
  // thing as above, but we have to do it multiple times
  for (size_type row = 0; row < matrix.m(); ++row)
    {
      long int cursor = Ap[row];
      for (size_type block = 0; block < matrix.n_block_cols(); ++block)
        {
          // find the next out-of-order element
          while ((cursor < Ap[row + 1] - 1) && (Ai[cursor] < Ai[cursor + 1]))
            ++cursor;

          // if there is none, then just go on
          if (cursor == Ap[row + 1] - 1)
            break;

          // otherwise swap this entry with successive ones as long as
          // necessary
          long int element = cursor;
          while ((element < Ap[row + 1] - 1) && (Ai[element] > Ai[element + 1]))
            {
              std::swap(Ai[element], Ai[element + 1]);

              std::swap(Ax[element], Ax[element + 1]);
              if (numbers::NumberTraits<number>::is_complex == true)
                std::swap(Az[cursor], Az[cursor + 1]);

              ++element;
            }
        }
    }
}



template <class Matrix>
void
SparseDirectUMFPACK::factorize(const Matrix &matrix)
{
  Assert(matrix.m() == matrix.n(), ExcNotQuadratic());

  clear();

  using number = typename Matrix::value_type;

  n_rows = matrix.m();
  n_cols = matrix.n();

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
  Ap.resize(N + 1);
  Ai.resize(matrix.n_nonzero_elements());
  Ax.resize(matrix.n_nonzero_elements());
  if (numbers::NumberTraits<number>::is_complex == true)
    Az.resize(matrix.n_nonzero_elements());

  // first fill row lengths array
  Ap[0] = 0;
  for (size_type row = 1; row <= N; ++row)
    Ap[row] = Ap[row - 1] + matrix.get_row_length(row - 1);
  Assert(static_cast<size_type>(Ap.back()) == Ai.size(), ExcInternalError());

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
        for (typename Matrix::const_iterator p = matrix.begin(row);
             p != matrix.end(row);
             ++p)
          {
            // write entry into the first free one for this row
            Ai[row_pointers[row]] = p->column();
            Ax[row_pointers[row]] = std::real(p->value());
            if (numbers::NumberTraits<number>::is_complex == true)
              Az[row_pointers[row]] = std::imag(p->value());

            // then move pointer ahead
            ++row_pointers[row];
          }
      }

    // at the end, we should have written all rows completely
    for (size_type i = 0; i < Ap.size() - 1; ++i)
      Assert(row_pointers[i] == Ap[i + 1], ExcInternalError());
  }

  // make sure that the elements in each row are sorted. we have to be more
  // careful for block sparse matrices, so ship this task out to a
  // different function
  sort_arrays(matrix);

  int status;
  if (numbers::NumberTraits<number>::is_complex == false)
    status = umfpack_dl_symbolic(N,
                                 N,
                                 Ap.data(),
                                 Ai.data(),
                                 Ax.data(),
                                 &symbolic_decomposition,
                                 control.data(),
                                 nullptr);
  else
    status = umfpack_zl_symbolic(N,
                                 N,
                                 Ap.data(),
                                 Ai.data(),
                                 Ax.data(),
                                 Az.data(),
                                 &symbolic_decomposition,
                                 control.data(),
                                 nullptr);
  AssertThrow(status == UMFPACK_OK,
              ExcUMFPACKError("umfpack_dl_symbolic", status));

  if (numbers::NumberTraits<number>::is_complex == false)
    status = umfpack_dl_numeric(Ap.data(),
                                Ai.data(),
                                Ax.data(),
                                symbolic_decomposition,
                                &numeric_decomposition,
                                control.data(),
                                nullptr);
  else
    status = umfpack_zl_numeric(Ap.data(),
                                Ai.data(),
                                Ax.data(),
                                Az.data(),
                                symbolic_decomposition,
                                &numeric_decomposition,
                                control.data(),
                                nullptr);
  AssertThrow(status == UMFPACK_OK,
              ExcUMFPACKError("umfpack_dl_numeric", status));

  umfpack_dl_free_symbolic(&symbolic_decomposition);
}



void
SparseDirectUMFPACK::solve(Vector<double> &rhs_and_solution,
                           const bool      transpose /*=false*/) const
{
  // make sure that some kind of factorize() call has happened before
  Assert(Ap.size() != 0, ExcNotInitialized());
  Assert(Ai.size() != 0, ExcNotInitialized());
  Assert(Ai.size() == Ax.size(), ExcNotInitialized());

  Assert(Az.size() == 0,
         ExcMessage("You have previously factored a matrix using this class "
                    "that had complex-valued entries. This then requires "
                    "applying the factored matrix to a complex-valued "
                    "vector, but you are only providing a real-valued vector "
                    "here."));

  Vector<double> rhs(rhs_and_solution.size());
  rhs = rhs_and_solution;

  // solve the system. note that since UMFPACK wants compressed column
  // storage instead of the compressed row storage format we use in
  // deal.II's SparsityPattern classes, we solve for UMFPACK's A^T instead

  // Conversely, if we solve for the transpose, we have to use UMFPACK_A
  // instead.
  const int status = umfpack_dl_solve(transpose ? UMFPACK_A : UMFPACK_At,
                                      Ap.data(),
                                      Ai.data(),
                                      Ax.data(),
                                      rhs_and_solution.begin(),
                                      rhs.begin(),
                                      numeric_decomposition,
                                      control.data(),
                                      nullptr);
  AssertThrow(status == UMFPACK_OK,
              ExcUMFPACKError("umfpack_dl_solve", status));
}



void
SparseDirectUMFPACK::solve(Vector<std::complex<double>> &rhs_and_solution,
                           const bool transpose /*=false*/) const
{
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
  // make sure that some kind of factorize() call has happened before
  Assert(Ap.size() != 0, ExcNotInitialized());
  Assert(Ai.size() != 0, ExcNotInitialized());
  Assert(Ai.size() == Ax.size(), ExcNotInitialized());

  // First see whether the matrix that was factorized was complex-valued.
  // If so, just apply the complex factorization to the vector.
  if (Az.size() != 0)
    {
      Assert(Ax.size() == Az.size(), ExcInternalError());

      // It would be nice if we could just present a pointer to the
      // first element of the complex-valued solution vector and let
      // UMFPACK fill both the real and imaginary parts of the solution
      // vector at that address. UMFPACK calls this the 'packed' format,
      // and in those cases it only takes one pointer to the entire
      // vector, rather than a pointer to the real and one pointer to
      // the imaginary parts of the vector. The problem is that if we
      // want to pack, then we also need to pack the matrix, and the
      // functions above have already decided that we don't want to pack
      // the matrix but instead deal with split format for the matrix,
      // and then UMFPACK decides that it can't deal with a split matrix
      // and a packed vector. We have to choose one or the other, not
      // mix.
      //
      // So create four vectors, one each for the real and imaginary parts
      // of the right hand side and solution.
      Vector<double> rhs_re(rhs_and_solution.size());
      Vector<double> rhs_im(rhs_and_solution.size());
      for (unsigned int i = 0; i < rhs_and_solution.size(); ++i)
        {
          rhs_re(i) = std::real(rhs_and_solution(i));
          rhs_im(i) = std::imag(rhs_and_solution(i));
        }

      Vector<double> solution_re(rhs_and_solution.size());
      Vector<double> solution_im(rhs_and_solution.size());

      // Solve the system. note that since UMFPACK wants compressed column
      // storage instead of the compressed row storage format we use in
      // deal.II's SparsityPattern classes, we solve for UMFPACK's A^T instead
      //
      // Conversely, if we solve for the transpose, we have to use UMFPACK_A
      // instead.
      //
      // Note that for the complex case here, the transpose is selected using
      // UMFPACK_Aat, not UMFPACK_At.
      const int status = umfpack_zl_solve(transpose ? UMFPACK_A : UMFPACK_Aat,
                                          Ap.data(),
                                          Ai.data(),
                                          Ax.data(),
                                          Az.data(),
                                          solution_re.data(),
                                          solution_im.data(),
                                          rhs_re.data(),
                                          rhs_im.data(),
                                          numeric_decomposition,
                                          control.data(),
                                          nullptr);
      AssertThrow(status == UMFPACK_OK,
                  ExcUMFPACKError("umfpack_dl_solve", status));

      // Now put things back together into the output vector
      for (unsigned int i = 0; i < rhs_and_solution.size(); ++i)
        rhs_and_solution(i) = {solution_re(i), solution_im(i)};
    }
  else
    {
      // We have factorized a real-valued matrix, but the rhs and solution
      // vectors are complex-valued. UMFPACK does not natively support this
      // case, but we can just apply the factorization to real and imaginary
      // parts of the right hand side separately
      const Vector<std::complex<double>> rhs = rhs_and_solution;

      // Get the real part of the right hand side, solve with it, and copy the
      // results into the result vector by just copying the real output
      // into the complex-valued result vector (which implies setting the
      // imaginary parts to zero):
      Vector<double> rhs_real_or_imag(rhs_and_solution.size());
      for (unsigned int i = 0; i < rhs.size(); ++i)
        rhs_real_or_imag(i) = std::real(rhs(i));

      solve(rhs_real_or_imag, transpose);

      rhs_and_solution = rhs_real_or_imag;

      // Then repeat the whole thing with the imaginary part. The copying step
      // is more complicated now because we can only touch the imaginary
      // component of the output vector:
      for (unsigned int i = 0; i < rhs.size(); ++i)
        rhs_real_or_imag(i) = std::imag(rhs(i));

      solve(rhs_real_or_imag, transpose);

      for (unsigned int i = 0; i < rhs.size(); ++i)
        rhs_and_solution(i).imag(rhs_real_or_imag(i));
    }

#  else

  (void)rhs_and_solution;
  (void)transpose;
  Assert(false,
         ExcMessage(
           "This function can't be called if deal.II has been configured "
           "with DEAL_II_WITH_COMPLEX_VALUES=FALSE."));
#  endif
}


void
SparseDirectUMFPACK::solve(BlockVector<double> &rhs_and_solution,
                           const bool           transpose /*=false*/) const
{
  // the UMFPACK functions want a contiguous array of elements, so
  // there is no way around copying data around. thus, just copy the
  // data into a regular vector and back
  Vector<double> tmp(rhs_and_solution.size());
  tmp = rhs_and_solution;
  solve(tmp, transpose);
  rhs_and_solution = tmp;
}



void
SparseDirectUMFPACK::solve(BlockVector<std::complex<double>> &rhs_and_solution,
                           const bool transpose /*=false*/) const
{
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
  // the UMFPACK functions want a contiguous array of elements, so
  // there is no way around copying data around. thus, just copy the
  // data into a regular vector and back
  Vector<std::complex<double>> tmp(rhs_and_solution.size());
  tmp = rhs_and_solution;
  solve(tmp, transpose);
  rhs_and_solution = tmp;

#  else
  (void)rhs_and_solution;
  (void)transpose;
  Assert(false,
         ExcMessage(
           "This function can't be called if deal.II has been configured "
           "with DEAL_II_WITH_COMPLEX_VALUES=FALSE."));
#  endif
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &  matrix,
                           Vector<double> &rhs_and_solution,
                           const bool      transpose /*=false*/)
{
  factorize(matrix);
  solve(rhs_and_solution, transpose);
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &                matrix,
                           Vector<std::complex<double>> &rhs_and_solution,
                           const bool                    transpose /*=false*/)
{
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
  factorize(matrix);
  solve(rhs_and_solution, transpose);

#  else

  (void)matrix;
  (void)rhs_and_solution;
  (void)transpose;
  Assert(false,
         ExcMessage(
           "This function can't be called if deal.II has been configured "
           "with DEAL_II_WITH_COMPLEX_VALUES=FALSE."));
#  endif
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &       matrix,
                           BlockVector<double> &rhs_and_solution,
                           const bool           transpose /*=false*/)
{
  factorize(matrix);
  solve(rhs_and_solution, transpose);
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &                     matrix,
                           BlockVector<std::complex<double>> &rhs_and_solution,
                           const bool transpose /*=false*/)
{
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
  factorize(matrix);
  solve(rhs_and_solution, transpose);

#  else

  (void)matrix;
  (void)rhs_and_solution;
  (void)transpose;
  Assert(false,
         ExcMessage(
           "This function can't be called if deal.II has been configured "
           "with DEAL_II_WITH_COMPLEX_VALUES=FALSE."));
#  endif
}


#else


SparseDirectUMFPACK::SparseDirectUMFPACK()
  : n_rows(0)
  , n_cols(0)
  , symbolic_decomposition(nullptr)
  , numeric_decomposition(nullptr)
  , control(0)
{}


void
SparseDirectUMFPACK::clear()
{}


template <class Matrix>
void
SparseDirectUMFPACK::factorize(const Matrix &)
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}


void
SparseDirectUMFPACK::solve(Vector<double> &, const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



void
SparseDirectUMFPACK::solve(Vector<std::complex<double>> &, const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



void
SparseDirectUMFPACK::solve(BlockVector<double> &, const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



void
SparseDirectUMFPACK::solve(BlockVector<std::complex<double>> &,
                           const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &, Vector<double> &, const bool)
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &,
                           Vector<std::complex<double>> &,
                           const bool)
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &, BlockVector<double> &, const bool)
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix &,
                           BlockVector<std::complex<double>> &,
                           const bool)
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions in doc/readme.html."));
}

#endif


template <class Matrix>
void
SparseDirectUMFPACK::initialize(const Matrix &M, const AdditionalData)
{
  this->factorize(M);
}


void
SparseDirectUMFPACK::vmult(Vector<double> &dst, const Vector<double> &src) const
{
  dst = src;
  this->solve(dst);
}



void
SparseDirectUMFPACK::vmult(BlockVector<double> &      dst,
                           const BlockVector<double> &src) const
{
  dst = src;
  this->solve(dst);
}


void
SparseDirectUMFPACK::Tvmult(Vector<double> &      dst,
                            const Vector<double> &src) const
{
  dst = src;
  this->solve(dst, /*transpose=*/true);
}



void
SparseDirectUMFPACK::Tvmult(BlockVector<double> &      dst,
                            const BlockVector<double> &src) const
{
  dst = src;
  this->solve(dst, /*transpose=*/true);
}

SparseDirectUMFPACK::size_type
SparseDirectUMFPACK::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

SparseDirectUMFPACK::size_type
SparseDirectUMFPACK::n() const
{
  Assert(n_cols != 0, ExcNotInitialized());
  return n_cols;
}


// explicit instantiations for SparseMatrixUMFPACK
#define InstantiateUMFPACK(MatrixType)                                     \
  template void SparseDirectUMFPACK::factorize(const MatrixType &);        \
  template void SparseDirectUMFPACK::solve(const MatrixType &,             \
                                           Vector<double> &,               \
                                           const bool);                    \
  template void SparseDirectUMFPACK::solve(const MatrixType &,             \
                                           Vector<std::complex<double>> &, \
                                           const bool);                    \
  template void SparseDirectUMFPACK::solve(const MatrixType &,             \
                                           BlockVector<double> &,          \
                                           const bool);                    \
  template void SparseDirectUMFPACK::solve(                                \
    const MatrixType &, BlockVector<std::complex<double>> &, const bool);  \
  template void SparseDirectUMFPACK::initialize(const MatrixType &,        \
                                                const AdditionalData)

// Instantiate everything for real-valued matrices
InstantiateUMFPACK(SparseMatrix<double>);
InstantiateUMFPACK(SparseMatrix<float>);
InstantiateUMFPACK(SparseMatrixEZ<double>);
InstantiateUMFPACK(SparseMatrixEZ<float>);
InstantiateUMFPACK(BlockSparseMatrix<double>);
InstantiateUMFPACK(BlockSparseMatrix<float>);

// Now also for complex-valued matrices
#ifdef DEAL_II_WITH_COMPLEX_VALUES
InstantiateUMFPACK(SparseMatrix<std::complex<double>>);
InstantiateUMFPACK(SparseMatrix<std::complex<float>>);
InstantiateUMFPACK(BlockSparseMatrix<std::complex<double>>);
InstantiateUMFPACK(BlockSparseMatrix<std::complex<float>>);
#endif

DEAL_II_NAMESPACE_CLOSE
