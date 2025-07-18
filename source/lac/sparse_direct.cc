// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <complex>
#include <vector>

#ifdef DEAL_II_WITH_UMFPACK
#  include <umfpack.h>
#endif

#ifdef DEAL_II_WITH_MUMPS
#  include <dmumps_c.h>
#endif

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  class SparseMatrix;
  namespace MPI
  {
    class SparseMatrix;
    class Vector;
  } // namespace MPI
} // namespace TrilinosWrappers
namespace PETScWrappers
{
  namespace MPI
  {
    class SparseMatrix;
    class Vector;
  } // namespace MPI
} // namespace PETScWrappers

namespace
{
  /**
   * For a given matrix, compute a reasonable "grain size" when parallelizing
   * some copying and sorting operations. The grain size is the minimal number
   * of rows each thread should work on. We define it by assuming that
   * dealing with 1000 matrix entries is a reasonable lower bound for parallel
   * operations.
   */
  template <typename SparseMatrixType>
  unsigned int
  parallel_grainsize(const SparseMatrixType &matrix)
  {
    const unsigned int avg_entries_per_row =
      matrix.n_nonzero_elements() / matrix.m();
    return std::max(1000 / avg_entries_per_row, 1u);
  }
} // namespace


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
    std::vector<types::suitesparse_index> tmp;
    tmp.swap(Ap);
  }

  {
    std::vector<types::suitesparse_index> tmp;
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
  parallel::apply_to_subranges(
    0,
    matrix.m(),
    [this](const size_type row_begin, const size_type row_end) {
      for (size_type row = row_begin; row < row_end; ++row)
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
    },
    parallel_grainsize(matrix));
}



template <typename number>
void
SparseDirectUMFPACK::sort_arrays(const SparseMatrixEZ<number> &matrix)
{
  // same thing for SparseMatrixEZ
  parallel::apply_to_subranges(
    0,
    matrix.m(),
    [this](const size_type row_begin, const size_type row_end) {
      for (size_type row = row_begin; row < row_end; ++row)
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
    },
    parallel_grainsize(matrix));
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
  parallel::apply_to_subranges(
    0,
    matrix.m(),
    [this, &matrix](const size_type row_begin, const size_type row_end) {
      for (size_type row = row_begin; row < row_end; ++row)
        {
          long int cursor = Ap[row];
          for (size_type block = 0; block < matrix.n_block_cols(); ++block)
            {
              // find the next out-of-order element
              while ((cursor < Ap[row + 1] - 1) &&
                     (Ai[cursor] < Ai[cursor + 1]))
                ++cursor;

              // if there is none, then just go on
              if (cursor == Ap[row + 1] - 1)
                break;

              // otherwise swap this entry with successive ones as long as
              // necessary
              long int element = cursor;
              while ((element < Ap[row + 1] - 1) &&
                     (Ai[element] > Ai[element + 1]))
                {
                  std::swap(Ai[element], Ai[element + 1]);

                  std::swap(Ax[element], Ax[element + 1]);
                  if (numbers::NumberTraits<number>::is_complex == true)
                    std::swap(Az[element], Az[element + 1]);

                  ++element;
                }
            }
        }
    },
    parallel_grainsize(matrix));
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
  // before moving on to the next row.
  parallel::apply_to_subranges(
    0,
    matrix.m(),
    [this, &matrix](const size_type row_begin, const size_type row_end) {
      for (size_type row = row_begin; row < row_end; ++row)
        {
          long int index = Ap[row];
          for (typename Matrix::const_iterator p = matrix.begin(row);
               p != matrix.end(row);
               ++p)
            {
              // write entry into the first free one for this row
              Ai[index] = p->column();
              Ax[index] = std::real(p->value());
              if (numbers::NumberTraits<number>::is_complex == true)
                Az[index] = std::imag(p->value());

              // then move pointer ahead
              ++index;
            }
          Assert(index == Ap[row + 1], ExcInternalError());
        }
    },
    parallel_grainsize(matrix));

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

  // Clean up before we deal with the error code from the calls above:
  umfpack_dl_free_symbolic(&symbolic_decomposition);
  if (status == UMFPACK_WARNING_singular_matrix)
    {
      // UMFPACK sometimes warns that the matrix is singular, but that a
      // factorization was successful nonetheless. Report this by
      // throwing an exception that can be caught at a higher level:
      AssertThrow(false,
                  ExcMessage(
                    "UMFPACK reports that the matrix is singular, "
                    "but that the factorization was successful anyway. "
                    "You can try and see whether you can still "
                    "solve a linear system with such a factorization "
                    "by catching and ignoring this exception, "
                    "though in practice this will typically not "
                    "work."));
    }
  else
    AssertThrow(status == UMFPACK_OK,
                ExcUMFPACKError("umfpack_dl_numeric", status));
}



void
SparseDirectUMFPACK::solve(Vector<double> &rhs_and_solution,
                           const bool      transpose /*=false*/) const
{
  // make sure that some kind of factorize() call has happened before
  Assert(Ap.size() != 0, ExcNotInitialized());
  Assert(Ai.size() != 0, ExcNotInitialized());
  Assert(Ai.size() == Ax.size(), ExcNotInitialized());

  Assert(Az.empty(),
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
SparseDirectUMFPACK::solve(const Matrix   &matrix,
                           Vector<double> &rhs_and_solution,
                           const bool      transpose /*=false*/)
{
  factorize(matrix);
  solve(rhs_and_solution, transpose);
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix                 &matrix,
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
SparseDirectUMFPACK::solve(const Matrix        &matrix,
                           BlockVector<double> &rhs_and_solution,
                           const bool           transpose /*=false*/)
{
  factorize(matrix);
  solve(rhs_and_solution, transpose);
}



template <class Matrix>
void
SparseDirectUMFPACK::solve(const Matrix                      &matrix,
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
      "installation instructions at https://dealii.org/current/readme.html"));
}


void
SparseDirectUMFPACK::solve(Vector<double> &, const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions at https://dealii.org/current/readme.html"));
}



void
SparseDirectUMFPACK::solve(Vector<std::complex<double>> &, const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions at https://dealii.org/current/readme.html"));
}



void
SparseDirectUMFPACK::solve(BlockVector<double> &, const bool) const
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need UMFPACK, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions at https://dealii.org/current/readme.html"));
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
      "installation instructions at https://dealii.org/current/readme.html"));
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
      "installation instructions at https://dealii.org/current/readme.html"));
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
      "installation instructions at https://dealii.org/current/readme.html"));
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
      "installation instructions at https://dealii.org/current/readme.html"));
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
      "installation instructions at https://dealii.org/current/readme.html"));
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
SparseDirectUMFPACK::vmult(BlockVector<double>       &dst,
                           const BlockVector<double> &src) const
{
  dst = src;
  this->solve(dst);
}


void
SparseDirectUMFPACK::Tvmult(Vector<double>       &dst,
                            const Vector<double> &src) const
{
  dst = src;
  this->solve(dst, /*transpose=*/true);
}



void
SparseDirectUMFPACK::Tvmult(BlockVector<double>       &dst,
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



#ifdef DEAL_II_WITH_MUMPS

SparseDirectMUMPS::SparseDirectMUMPS(const AdditionalData &data,
                                     const MPI_Comm       &communicator)
  : additional_data(data)
  , mpi_communicator(communicator)
{
  // Initialize MUMPS instance:
  id.job = -1;
  id.par = 1;

  Assert(!(additional_data.symmetric == false &&
           additional_data.posdef == true),
         ExcMessage(
           "You can't have a positive definite matrix that is not symmetric."));

  if (additional_data.symmetric == true)
    {
      if (additional_data.posdef == true)
        id.sym = 1;
      else
        id.sym = 2;
    }
  else
    id.sym = 0;

  id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(mpi_communicator);
  dmumps_c(&id);

  if (additional_data.output_details == false)
    {
      // No outputs
      id.icntl[0] = -1;
      id.icntl[1] = -1;
      id.icntl[2] = -1;
      id.icntl[3] = 0;
    }

  if (additional_data.error_statistics == true)
    id.icntl[10] = 2;

  if (additional_data.blr_factorization == true)
    {
      id.icntl[34] = 2;
      if (additional_data.blr.blr_ucfs == true)
        id.icntl[35] = 1;
      Assert(additional_data.blr.lowrank_threshold > 0,
             ExcMessage("Lowrank threshold must be positive."));
      id.cntl[6] = additional_data.blr.lowrank_threshold;
    }
}



SparseDirectMUMPS::~SparseDirectMUMPS()
{
  // MUMPS destructor
  id.job = -2;
  dmumps_c(&id);
}



template <class Matrix>
void
SparseDirectMUMPS::initialize_matrix(const Matrix &matrix)
{
  Assert(matrix.n() == matrix.m(), ExcMessage("Matrix needs to be square."));

  n    = matrix.n();
  id.n = n;

  if constexpr (std::is_same_v<Matrix, SparseMatrix<double>>)
    {
      // Serial matrix: hand over matrix to MUMPS as centralized assembled
      // matrix
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          // number of nonzero elements in matrix
          nnz = matrix.n_actually_nonzero_elements();

          // representation of the matrix
          a = std::make_unique<double[]>(nnz);

          // matrix indices pointing to the row and column dimensions
          // respectively of the matrix representation above (a): ie. a[k] is
          // the matrix element (irn[k], jcn[k])
          irn = std::make_unique<MUMPS_INT[]>(nnz);
          jcn = std::make_unique<MUMPS_INT[]>(nnz);

          size_type n_non_zero_elements = 0;

          // loop over the elements of the matrix row by row, as suggested in
          // the documentation of the sparse matrix iterator class
          if (additional_data.symmetric == true)
            {
              for (size_type row = 0; row < matrix.m(); ++row)
                {
                  for (typename Matrix::const_iterator ptr = matrix.begin(row);
                       ptr != matrix.end(row);
                       ++ptr)
                    if (std::abs(ptr->value()) > 0.0 && ptr->column() >= row)
                      {
                        a[n_non_zero_elements]   = ptr->value();
                        irn[n_non_zero_elements] = row + 1;
                        jcn[n_non_zero_elements] = ptr->column() + 1;

                        ++n_non_zero_elements;
                      }
                }
            }
          else
            {
              for (size_type row = 0; row < matrix.m(); ++row)
                {
                  for (typename Matrix::const_iterator ptr = matrix.begin(row);
                       ptr != matrix.end(row);
                       ++ptr)
                    if (std::abs(ptr->value()) > 0.0)
                      {
                        a[n_non_zero_elements]   = ptr->value();
                        irn[n_non_zero_elements] = row + 1;
                        jcn[n_non_zero_elements] = ptr->column() + 1;
                        ++n_non_zero_elements;
                      }
                }
            }
          id.n   = n;
          id.nnz = n_non_zero_elements;
          id.irn = irn.get();
          id.jcn = jcn.get();
          id.a   = a.get();
        }
    }
  else if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                     std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    {
      int result;
      MPI_Comm_compare(mpi_communicator,
                       matrix.get_mpi_communicator(),
                       &result);
      AssertThrow(result == MPI_IDENT,
                  ExcMessage("The matrix communicator must match the MUMPS "
                             "communicator."));

      // Distributed matrix case
      id.icntl[17]               = 3; // distributed matrix assembly
      id.nnz                     = matrix.n_nonzero_elements();
      nnz                        = id.nnz;
      size_type n_non_zero_local = 0;

      // Get the range of rows owned by this process
      locally_owned_rows        = matrix.locally_owned_range_indices();
      size_type local_non_zeros = 0;

      if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix>)
        {
          const auto &trilinos_matrix = matrix.trilinos_matrix();
          local_non_zeros             = trilinos_matrix.NumMyNonzeros();
        }
      else if constexpr (std::is_same_v<Matrix,
                                        PETScWrappers::MPI::SparseMatrix>)
        {
#  ifdef DEAL_II_WITH_PETSC
          Mat &petsc_matrix =
            const_cast<PETScWrappers::MPI::SparseMatrix &>(matrix)
              .petsc_matrix();
          MatInfo info;
          MatGetInfo(petsc_matrix, MAT_LOCAL, &info);
          local_non_zeros = (size_type)info.nz_used;
#  endif
        }


      // We allocate enough entries for the general, nonsymmetric, case. If case
      // of a symmetric matrix, we will end up with fewer entries.
      irn = std::make_unique<MUMPS_INT[]>(local_non_zeros);
      jcn = std::make_unique<MUMPS_INT[]>(local_non_zeros);
      a   = std::make_unique<double[]>(local_non_zeros);
      irhs_loc.resize(locally_owned_rows.n_elements());

      if (additional_data.symmetric == true)
        {
          if constexpr (std::is_same_v<Matrix,
                                       PETScWrappers::MPI::SparseMatrix>)
            {
#  ifdef DEAL_II_WITH_PETSC
              Mat &petsc_matrix =
                const_cast<PETScWrappers::MPI::SparseMatrix &>(matrix)
                  .petsc_matrix();

              PetscInt rstart, rend;
              MatGetOwnershipRange(petsc_matrix, &rstart, &rend);
              for (PetscInt i = rstart; i < rend; i++)
                {
                  PetscInt           n_cols;
                  const PetscInt    *cols;
                  const PetscScalar *values;
                  MatGetRow(petsc_matrix, i, &n_cols, &cols, &values);

                  for (PetscInt j = 0; j < n_cols; j++)
                    {
                      if (cols[j] >= i)
                        {
                          irn[n_non_zero_local] = i + 1;
                          jcn[n_non_zero_local] = cols[j] + 1;
                          a[n_non_zero_local]   = values[j];

                          // Count local non-zeros
                          n_non_zero_local++;
                        }
                    }
                  MatRestoreRow(petsc_matrix, i, &n_cols, &cols, &values);

                  // Store the row index for the rhs vector
                  const types::global_cell_index local_index =
                    locally_owned_rows.index_within_set(i);
                  irhs_loc[local_index] = i + 1;
                }

              id.a_loc = a.get();
#  endif
            }
          else if constexpr (std::is_same_v<Matrix,
                                            TrilinosWrappers::SparseMatrix>)
            {
              const auto &trilinos_matrix = matrix.trilinos_matrix();
              for (int local_row = 0; local_row < trilinos_matrix.NumMyRows();
                   ++local_row)
                {
                  int     num_entries;
                  double *values;
                  int    *local_cols;
                  int     ierr = trilinos_matrix.ExtractMyRowView(local_row,
                                                              num_entries,
                                                              values,
                                                              local_cols);
                  (void)ierr;
                  Assert(
                    ierr == 0,
                    ExcMessage(
                      "Error extracting global row view from Trilinos matrix. Error code " +
                      std::to_string(ierr) + "."));

                  int global_row = trilinos_matrix.GRID(local_row);

                  for (int j = 0; j < num_entries; ++j)
                    {
                      if (trilinos_matrix.GCID(local_cols[j]) >= global_row)
                        {
                          irn[n_non_zero_local] = global_row + 1;
                          jcn[n_non_zero_local] =
                            trilinos_matrix.GCID(local_cols[j]) + 1;
                          a[n_non_zero_local] = values[j];

                          // Count local non-zeros
                          n_non_zero_local++;
                        }
                    }

                  // Store the row index for the rhs vector
                  irhs_loc[local_row] = global_row + 1;
                }
              id.a_loc = a.get();
            }
          else
            {
              DEAL_II_NOT_IMPLEMENTED();
            }
        }
      else
        {
          // Unsymmetric case
          if constexpr (std::is_same_v<Matrix,
                                       PETScWrappers::MPI::SparseMatrix>)
            {
#  ifdef DEAL_II_WITH_PETSC
              Mat &petsc_matrix =
                const_cast<PETScWrappers::MPI::SparseMatrix &>(matrix)
                  .petsc_matrix();

              PetscInt rstart, rend;
              MatGetOwnershipRange(petsc_matrix, &rstart, &rend);
              for (PetscInt i = rstart; i < rend; i++)
                {
                  PetscInt           n_cols;
                  const PetscInt    *cols;
                  const PetscScalar *values;
                  MatGetRow(petsc_matrix, i, &n_cols, &cols, &values);

                  for (PetscInt j = 0; j < n_cols; j++)
                    {
                      irn[n_non_zero_local] = i + 1;
                      jcn[n_non_zero_local] = cols[j] + 1;
                      a[n_non_zero_local]   = values[j];

                      // Count local non-zeros
                      n_non_zero_local++;
                    }
                  MatRestoreRow(petsc_matrix, i, &n_cols, &cols, &values);

                  // Store the row index for the rhs vector
                  const types::global_cell_index local_index =
                    locally_owned_rows.index_within_set(i);
                  irhs_loc[local_index] = i + 1;
                }

              id.a_loc = a.get();
#  endif
            }
          else if constexpr (std::is_same_v<Matrix,
                                            TrilinosWrappers::SparseMatrix>)
            {
              const auto &trilinos_matrix = matrix.trilinos_matrix();
              for (int local_row = 0; local_row < trilinos_matrix.NumMyRows();
                   ++local_row)
                {
                  int     num_entries;
                  double *values;
                  int    *local_cols;
                  int     ierr = trilinos_matrix.ExtractMyRowView(local_row,
                                                              num_entries,
                                                              values,
                                                              local_cols);
                  (void)ierr;
                  Assert(
                    ierr == 0,
                    ExcMessage(
                      "Error extracting global row view from Trilinos matrix. Error code " +
                      std::to_string(ierr) + "."));

                  int global_row = trilinos_matrix.GRID(local_row);

                  for (int j = 0; j < num_entries; ++j)
                    {
                      irn[n_non_zero_local] = global_row + 1;
                      jcn[n_non_zero_local] =
                        trilinos_matrix.GCID(local_cols[j]) + 1;
                      a[n_non_zero_local] = values[j];

                      // Count local non-zeros
                      n_non_zero_local++;
                    }

                  // Store the row index for the rhs vector
                  irhs_loc[local_row] = global_row + 1;
                }
              id.a_loc = a.get();
            }
          else
            {
              DEAL_II_NOT_IMPLEMENTED();
            }
        }

      // Hand over local arrays to MUMPS
      id.nnz_loc  = n_non_zero_local;
      id.irn_loc  = irn.get();
      id.jcn_loc  = jcn.get();
      id.a_loc    = a.get();
      id.irhs_loc = irhs_loc.data();

      // rhs parameters
      id.icntl[19] = 10; // distributed rhs
      id.icntl[20] = 0;  // centralized solution, stored on rank 0 by MUMPS
      id.nrhs      = 1;
      id.lrhs_loc  = n;
      id.nloc_rhs  = locally_owned_rows.n_elements();
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
    }
}



void
SparseDirectMUMPS::copy_rhs_to_mumps(const Vector<double> &new_rhs) const
{
  Assert(n == new_rhs.size(),
         ExcMessage("Matrix size and rhs length must be equal."));

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      rhs.resize(n);
      for (size_type i = 0; i < n; ++i)
        rhs[i] = new_rhs(i);

      id.rhs = &rhs[0];
    }
}



void
SparseDirectMUMPS::copy_solution(Vector<double> &vector) const
{
  Assert(n == vector.size(),
         ExcMessage("Matrix size and solution vector length must be equal."));
  Assert(n == rhs.size(),
         ExcMessage("Class not initialized with a rhs vector."));

  // Copy solution into the given vector
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      for (size_type i = 0; i < n; ++i)
        vector(i) = rhs[i];

      rhs.resize(0); // remove rhs again
    }
}



template <class Matrix>
void
SparseDirectMUMPS::initialize(const Matrix &matrix)
{
  // Initialize MUMPS instance:
  initialize_matrix(matrix);

  // Start analysis + factorization
  id.job = 4;
  dmumps_c(&id);
}



template <typename VectorType>
void
SparseDirectMUMPS::vmult(VectorType &dst, const VectorType &src) const
{
  // Check that the matrix has at least one nonzero element:
  Assert(nnz != 0, ExcNotInitialized());
  Assert(n == dst.size(), ExcMessage("Destination vector has the wrong size."));
  Assert(n == src.size(), ExcMessage("Source vector has the wrong size."));


  if constexpr (std::is_same_v<VectorType, Vector<double>>)
    {
      // Centralized assembly for serial vectors.

      // Hand over right-hand side
      copy_rhs_to_mumps(src);

      // Start solver
      id.job = 3;
      dmumps_c(&id);
      copy_solution(dst);
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::Vector> ||
                     std::is_same_v<VectorType, PETScWrappers::MPI::Vector>)
    {
      if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
        id.rhs_loc = const_cast<double *>(src.begin());
      else if constexpr (std::is_same_v<VectorType, PETScWrappers::MPI::Vector>)
        {
#  ifdef DEAL_II_WITH_PETSC
          PetscScalar *local_array;
          VecGetArray(
            const_cast<PETScWrappers::MPI::Vector &>(src).petsc_vector(),
            &local_array);
          id.rhs_loc = local_array;
          VecRestoreArray(
            const_cast<PETScWrappers::MPI::Vector &>(src).petsc_vector(),
            &local_array);
#  endif
        }


      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          rhs.resize(id.lrhs_loc);
          id.rhs = rhs.data();
        }

      // Start solver
      id.job = 3;
      dmumps_c(&id);

      // Copy solution into the given vector
      // For MUMPS with centralized solution (icntl[20]=0), the solution is only
      // on the root process (0) and needs to be distributed to all processes
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          // Set all the values in the dst vector
          for (size_type i = 0; i < n; ++i)
            dst[i] = rhs[i];
        }

      // Ensure the distributed vector has the correct values everywhere
      dst.compress(VectorOperation::insert);

      rhs.resize(0); // remove rhs again
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
    }
}


template <typename VectorType>
void
SparseDirectMUMPS::Tvmult(VectorType &dst, const VectorType &src) const
{
  // The matrix has at least one nonzero element:
  Assert(nnz != 0, ExcNotInitialized());
  Assert(n == dst.size(), ExcMessage("Destination vector has the wrong size."));
  Assert(n == src.size(), ExcMessage("Source vector has the wrong size."));

  id.icntl[8] = 2; // transpose
  vmult(dst, src);
  id.icntl[8] = 1; // reset to default
}



int *
SparseDirectMUMPS::get_icntl()
{
  return id.icntl;
}



#else


SparseDirectMUMPS::SparseDirectMUMPS(const AdditionalData &, const MPI_Comm &)
  : mpi_communicator(MPI_COMM_SELF)
{
  AssertThrow(
    false,
    ExcMessage(
      "To call this function you need MUMPS, but you configured deal.II "
      "without passing the necessary switch to 'cmake'. Please consult the "
      "installation instructions at https://dealii.org/current/readme.html"));
}



SparseDirectMUMPS::~SparseDirectMUMPS()
{}


#endif // DEAL_II_WITH_MUMPS



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

// explicit instantiations for SparseDirectMUMPS
#ifdef DEAL_II_WITH_MUMPS

#  define InstantiateMUMPSMatVec(VECTOR)                                    \
    template void SparseDirectMUMPS::vmult(VECTOR &, const VECTOR &) const; \
    template void SparseDirectMUMPS::Tvmult(VECTOR &, const VECTOR &) const;
#  ifdef DEAL_II_WITH_TRILINOS
InstantiateMUMPSMatVec(TrilinosWrappers::MPI::Vector)
#  endif
#  ifdef DEAL_II_WITH_PETSC
  InstantiateMUMPSMatVec(PETScWrappers::MPI::Vector)
#  endif
    InstantiateMUMPSMatVec(Vector<double>)

#  define InstantiateMUMPS(MATRIX) \
    template void SparseDirectMUMPS::initialize(const MATRIX &);

      InstantiateMUMPS(SparseMatrix<double>)
        InstantiateMUMPS(SparseMatrix<float>)
#  ifdef DEAL_II_WITH_TRILINOS
          InstantiateMUMPS(TrilinosWrappers::SparseMatrix)
#  endif
#  ifdef DEAL_II_WITH_PETSC
            InstantiateMUMPS(PETScWrappers::SparseMatrix)
              InstantiateMUMPS(PETScWrappers::MPI::SparseMatrix)
#  endif
  // InstantiateMUMPS(SparseMatrixEZ<double>)
  // InstantiateMUMPS(SparseMatrixEZ<float>)
  InstantiateMUMPS(BlockSparseMatrix<double>)
    InstantiateMUMPS(BlockSparseMatrix<float>)
#endif

      DEAL_II_NAMESPACE_CLOSE
