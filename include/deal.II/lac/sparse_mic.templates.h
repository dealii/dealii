// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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

#ifndef __deal2__sparse_mic_templates_h
#define __deal2__sparse_mic_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/sparse_mic.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
SparseMIC<number>::SparseMIC ()
  :
  diag(0),
  inv_diag(0),
  inner_sums(0)
{}



template <typename number>
SparseMIC<number>::SparseMIC (const SparsityPattern &sparsity)
  :
  diag(0),
  inv_diag(0),
  inner_sums(0)
{
  SparseMatrix<number>::reinit (sparsity);
}


template <typename number>
SparseMIC<number>::~SparseMIC()
{
  clear();
}


template <typename number>
void SparseMIC<number>::clear()
{
  {
    std::vector<number> tmp;
    tmp.swap (diag);
  }
  {
    std::vector<number> tmp;
    tmp.swap (inv_diag);
  }
  {
    std::vector<number> tmp;
    tmp.swap (inner_sums);
  }

  SparseLUDecomposition<number>::clear();
}


template <typename number>
template <typename somenumber>
inline
void SparseMIC<number>::initialize (const SparseMatrix<somenumber> &matrix,
                                    const AdditionalData &data)
{
  SparseLUDecomposition<number>::initialize(matrix, data);

  decompose(matrix, data.strengthen_diagonal);
}



template <typename number>
void SparseMIC<number>::reinit (const SparsityPattern &sparsity)
{
  {
    std::vector<number> tmp;
    tmp.swap (diag);
  }
  {
    std::vector<number> tmp;
    tmp.swap (inv_diag);
  }
  {
    std::vector<number> tmp;
    tmp.swap (inner_sums);
  }

  SparseMatrix<number>::reinit(sparsity);
  this->decomposed = false;
}



template <typename number>
template <typename somenumber>
void SparseMIC<number>::decompose (const SparseMatrix<somenumber> &matrix,
                                   const double                    strengthen_diagonal)
{
  SparseLUDecomposition<number>::decompose(matrix, strengthen_diagonal);

  Assert (matrix.m()==matrix.n(), ExcNotQuadratic ());
  Assert (this->m()==this->n(),   ExcNotQuadratic ());
  Assert (matrix.m()==this->m(),  ExcDimensionMismatch(matrix.m(), this->m()));

  Assert (strengthen_diagonal>=0, ExcInvalidStrengthening (strengthen_diagonal));

  if (strengthen_diagonal > 0)
    this->strengthen_diagonal_impl ();

  // MIC implementation: (S. Margenov lectures)
  // x[i] = a[i][i] - sum(k=1, i-1,
  //              a[i][k]/x[k]*sum(j=k+1, N, a[k][j]))

  // TODO: for sake of simplicity,
  // those are placed here. A better
  // implementation would store this
  // values in the underlying sparse
  // matrix itself.
  diag.resize (this->m());
  inv_diag.resize (this->m());
  inner_sums.resize (this->m());

  // precalc sum(j=k+1, N, a[k][j]))
  for (size_type row=0; row<this->m(); row++)
    inner_sums[row] = get_rowsum(row);

  for (size_type row=0; row<this->m(); row++)
    {
      const number temp = this->begin(row)->value();
      number temp1 = 0;

      // work on the lower left part of the matrix. we know
      // it's symmetric, so we can work with this alone
      for (typename SparseMatrix<somenumber>::const_iterator
           p = matrix.begin(row)+1;
           (p != matrix.end(row)) && (p->column() < row);
           ++p)
        temp1 += p->value() / diag[p->column()] * inner_sums[p->column()];

      Assert(temp-temp1 > 0, ExcStrengthenDiagonalTooSmall());
      diag[row] = temp - temp1;

      inv_diag[row] = 1.0/diag[row];
    }
}



template <typename number>
inline number
SparseMIC<number>::get_rowsum (const size_type row) const
{
  Assert(this->m()==this->n(), ExcNotQuadratic());

  number rowsum = 0;
  for (typename SparseMatrix<number>::const_iterator
       p = this->begin(row)+1;
       p != this->end(row); ++p)
    if (p->column() > row)
      rowsum += p->value();

  return rowsum;
}



template <typename number>
template <typename somenumber>
void
SparseMIC<number>::vmult (Vector<somenumber>       &dst,
                          const Vector<somenumber> &src) const
{
  SparseLUDecomposition<number>::vmult (dst, src);
  Assert (dst.size() == src.size(), ExcDimensionMismatch(dst.size(), src.size()));
  Assert (dst.size() == this->m(), ExcDimensionMismatch(dst.size(), this->m()));

  const size_type N=dst.size();
  // We assume the underlying matrix A is: A = X - L - U, where -L and -U are
  // strictly lower- and upper- diagonal parts of the system.
  //
  // Solve (X-L)X{-1}(X-U) x = b in 3 steps:
  dst = src;
  for (size_type row=0; row<N; ++row)
    {
      // Now: (X-L)u = b

      // get start of this row. skip
      // the diagonal element
      for (typename SparseMatrix<number>::const_iterator
           p = this->begin(row)+1;
           (p != this->end(row)) && (p->column() < row);
           ++p)
        dst(row) -= p->value() * dst(p->column());

      dst(row) *= inv_diag[row];
    }

  // Now: v = Xu
  for (size_type row=0; row<N; row++)
    dst(row) *= diag[row];

  // x = (X-U)v
  for (int row=N-1; row>=0; --row)
    {
      // get end of this row
      for (typename SparseMatrix<number>::const_iterator
           p = this->begin(row)+1;
           p != this->end(row);
           ++p)
        if (p->column() > static_cast<size_type>(row))
          dst(row) -= p->value() * dst(p->column());

      dst(row) *= inv_diag[row];
    }
}



template <typename number>
std::size_t
SparseMIC<number>::memory_consumption () const
{
  return (SparseLUDecomposition<number>::memory_consumption () +
          MemoryConsumption::memory_consumption(diag) +
          MemoryConsumption::memory_consumption(inv_diag) +
          MemoryConsumption::memory_consumption(inner_sums));
}



DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__sparse_mic_templates_h
