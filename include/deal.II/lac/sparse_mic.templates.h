// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_sparse_mic_templates_h
#define dealii_sparse_mic_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/sparse_mic.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
SparseMIC<number>::~SparseMIC()
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  SparseMIC<number>::clear();
}


template <typename number>
void
SparseMIC<number>::clear()
{
  {
    std::vector<number> tmp;
    tmp.swap(diag);
  }
  {
    std::vector<number> tmp;
    tmp.swap(inv_diag);
  }
  {
    std::vector<number> tmp;
    tmp.swap(inner_sums);
  }

  SparseLUDecomposition<number>::clear();
}


template <typename number>
template <typename somenumber>
inline void
SparseMIC<number>::initialize(const SparseMatrix<somenumber> &matrix,
                              const AdditionalData           &data)
{
  Assert(matrix.m() == matrix.n(), ExcNotQuadratic());
  Assert(data.strengthen_diagonal >= 0,
         ExcInvalidStrengthening(data.strengthen_diagonal));

  SparseLUDecomposition<number>::initialize(matrix, data);
  this->strengthen_diagonal = data.strengthen_diagonal;
  this->prebuild_lower_bound();
  this->copy_from(matrix);

  Assert(this->m() == this->n(), ExcNotQuadratic());
  Assert(matrix.m() == this->m(), ExcDimensionMismatch(matrix.m(), this->m()));

  if (data.strengthen_diagonal > 0)
    this->strengthen_diagonal_impl();

  // MIC implementation: (S. Margenov lectures)
  // x[i] = a[i][i] - sum(k=1, i-1,
  //              a[i][k]/x[k]*sum(j=k+1, N, a[k][j]))

  // TODO: for sake of simplicity,
  // those are placed here. A better
  // implementation would store this
  // values in the underlying sparse
  // matrix itself.
  diag.resize(this->m());
  inv_diag.resize(this->m());
  inner_sums.resize(this->m());

  // precalc sum(j=k+1, N, a[k][j]))
  for (size_type row = 0; row < this->m(); ++row)
    inner_sums[row] = get_rowsum(row);

  for (size_type row = 0; row < this->m(); ++row)
    {
      const number temp  = this->begin(row)->value();
      number       temp1 = 0;

      // work on the lower left part of the matrix. we know
      // it's symmetric, so we can work with this alone
      for (typename SparseMatrix<somenumber>::const_iterator p =
             matrix.begin(row) + 1;
           (p != matrix.end(row)) && (p->column() < row);
           ++p)
        temp1 += p->value() / diag[p->column()] * inner_sums[p->column()];

      Assert(temp - temp1 > 0, ExcStrengthenDiagonalTooSmall());
      diag[row] = temp - temp1;

      inv_diag[row] = 1.0 / diag[row];
    }
}



template <typename number>
inline number
SparseMIC<number>::get_rowsum(const size_type row) const
{
  Assert(this->m() == this->n(), ExcNotQuadratic());

  number rowsum = 0;
  for (typename SparseMatrix<number>::const_iterator p = this->begin(row) + 1;
       p != this->end(row);
       ++p)
    if (p->column() > row)
      rowsum += p->value();

  return rowsum;
}



template <typename number>
template <typename somenumber>
void
SparseMIC<number>::vmult(Vector<somenumber>       &dst,
                         const Vector<somenumber> &src) const
{
  Assert(dst.size() == src.size(),
         ExcDimensionMismatch(dst.size(), src.size()));
  Assert(dst.size() == this->m(), ExcDimensionMismatch(dst.size(), this->m()));

  const size_type N = dst.size();
  // We assume the underlying matrix A is: A = X - L - U, where -L and -U are
  // strictly lower- and upper- diagonal parts of the system.
  //
  // Solve (X-L)X{-1}(X-U) x = b in 3 steps:
  dst = src;
  for (size_type row = 0; row < N; ++row)
    {
      // Now: (X-L)u = b

      // get start of this row. skip
      // the diagonal element
      for (typename SparseMatrix<number>::const_iterator p =
             this->begin(row) + 1;
           (p != this->end(row)) && (p->column() < row);
           ++p)
        dst(row) -= p->value() * dst(p->column());

      dst(row) *= inv_diag[row];
    }

  // Now: v = Xu
  for (size_type row = 0; row < N; ++row)
    dst(row) *= diag[row];

  // x = (X-U)v
  for (int row = N - 1; row >= 0; --row)
    {
      // get end of this row
      for (typename SparseMatrix<number>::const_iterator p =
             this->begin(row) + 1;
           p != this->end(row);
           ++p)
        if (p->column() > static_cast<size_type>(row))
          dst(row) -= p->value() * dst(p->column());

      dst(row) *= inv_diag[row];
    }
}


// Exists for full compatibility with the LinearOperator class
template <typename number>
template <typename somenumber>
void
SparseMIC<number>::Tvmult(Vector<somenumber> & /*dst*/,
                          const Vector<somenumber> & /*src*/) const
{
  AssertThrow(false, ExcNotImplemented());
}



template <typename number>
std::size_t
SparseMIC<number>::memory_consumption() const
{
  return (SparseLUDecomposition<number>::memory_consumption() +
          MemoryConsumption::memory_consumption(diag) +
          MemoryConsumption::memory_consumption(inv_diag) +
          MemoryConsumption::memory_consumption(inner_sums));
}



DEAL_II_NAMESPACE_CLOSE

#endif // dealii_sparse_mic_templates_h
