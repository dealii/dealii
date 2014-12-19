// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__sparse_ilu_templates_h
#define __deal2__sparse_ilu_templates_h



#include <deal.II/base/config.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_ilu.h>

#include <algorithm>
#include <cmath>


DEAL_II_NAMESPACE_OPEN

template <typename number>
SparseILU<number>::SparseILU ()
{}



template <typename number>
SparseILU<number>::SparseILU (const SparsityPattern &sparsity)
{
  SparseMatrix<number>::reinit(sparsity);
}




template <typename number>
template <typename somenumber>
void SparseILU<number>::initialize (const SparseMatrix<somenumber> &matrix,
                                    const AdditionalData &data)
{
  SparseLUDecomposition<number>::initialize(matrix, data);

  decompose(matrix, data.strengthen_diagonal);
}



template <typename number>
template <typename somenumber>
void SparseILU<number>::decompose (const SparseMatrix<somenumber> &matrix,
                                   const double strengthen_diagonal)
{
  Assert (matrix.m()==matrix.n(), ExcNotQuadratic ());
  Assert (this->m()==this->n(),   ExcNotQuadratic ());
  Assert (matrix.m()==this->m(),  ExcDimensionMismatch(matrix.m(), this->m()));

  Assert (strengthen_diagonal>=0,
          ExcInvalidStrengthening (strengthen_diagonal));

  SparseLUDecomposition<number>::decompose (matrix, strengthen_diagonal);

  if (strengthen_diagonal>0)
    this->strengthen_diagonal_impl();

  // in the following, we implement algorithm 10.4 in the book by Saad by
  // translating in essence the algorithm given at the end of section 10.3.2,
  // using the names of variables used there
  const SparsityPattern     &sparsity = this->get_sparsity_pattern();
  const std::size_t *const ia    = sparsity.rowstart;
  const size_type *const ja      = sparsity.colnums;

  number *luval = this->SparseMatrix<number>::val;

  const size_type N = this->m();
  size_type jrow = 0;

  std::vector<size_type> iw (N, numbers::invalid_size_type);

  for (size_type k=0; k<N; ++k)
    {
      const size_type j1 = ia[k],
                      j2 = ia[k+1]-1;

      for (size_type j=j1; j<=j2; ++j)
        iw[ja[j]] = j;

      // the algorithm in the book works on the elements of row k left of the
      // diagonal. however, since we store the diagonal element at the first
      // position, start at the element after the diagonal and run as long as
      // we don't walk into the right half
      size_type j = j1+1;

      // pathological case: the current row of the matrix has only the
      // diagonal entry. then we have nothing to do.
      if (j > j2)
        goto label_200;

label_150:

      jrow = ja[j];
      if (jrow >= k)
        goto label_200;

      // actual computations:
      {
        number t1 = luval[j] * luval[ia[jrow]];
        luval[j] = t1;

        // jj runs from just right of the diagonal to the end of the row
        size_type jj = ia[jrow]+1;
        while (ja[jj] < jrow)
          ++jj;
        for (; jj<ia[jrow+1]; ++jj)
          {
            const size_type jw = iw[ja[jj]];
            if (jw != numbers::invalid_size_type)
              luval[jw] -= t1 * luval[jj];
          }

        ++j;
        if (j<=j2)
          goto label_150;
      }

label_200:

      // in the book there is an assertion that we have hit the diagonal
      // element, i.e. that jrow==k. however, we store the diagonal element at
      // the front, so jrow must actually be larger than k or j is already in
      // the next row
      Assert ((jrow > k) || (j==ia[k+1]), ExcInternalError());

      // now we have to deal with the diagonal element. in the book it is
      // located at position 'j', but here we use the convention of storing
      // the diagonal element first, so instead of j we use uptr[k]=ia[k]
      Assert (luval[ia[k]] != 0, ExcInternalError());

      luval[ia[k]] = 1./luval[ia[k]];

      for (size_type j=j1; j<=j2; ++j)
        iw[ja[j]] = numbers::invalid_size_type;
    }
}




template <typename number>
template <typename somenumber>
void SparseILU<number>::vmult (Vector<somenumber>       &dst,
                               const Vector<somenumber> &src) const
{
  Assert (dst.size() == src.size(), ExcDimensionMismatch(dst.size(), src.size()));
  Assert (dst.size() == this->m(), ExcDimensionMismatch(dst.size(), this->m()));

  const size_type N=dst.size();
  const std::size_t *const rowstart_indices
    = this->get_sparsity_pattern().rowstart;
  const size_type *const column_numbers
    = this->get_sparsity_pattern().colnums;

  // solve LUx=b in two steps:
  // first Ly = b, then
  //       Ux = y
  //
  // first a forward solve. since
  // the diagonal values of L are
  // one, there holds
  // y_i = b_i
  //       - sum_{j=0}^{i-1} L_{ij}y_j
  // we split the y_i = b_i off and
  // perform it at the outset of the
  // loop
  dst = src;
  for (size_type row=0; row<N; ++row)
    {
      // get start of this row. skip the
      // diagonal element
      const size_type *const rowstart = &column_numbers[rowstart_indices[row]+1];
      // find the position where the part
      // right of the diagonal starts
      const size_type *const first_after_diagonal = this->prebuilt_lower_bound[row];

      somenumber dst_row = dst(row);
      const number *luval = this->SparseMatrix<number>::val +
                            (rowstart - column_numbers);
      for (const size_type *col=rowstart; col!=first_after_diagonal; ++col, ++luval)
        dst_row -= *luval * dst(*col);
      dst(row) = dst_row;
    }

  // now the backward solve. same
  // procedure, but we need not set
  // dst before, since this is already
  // done.
  //
  // note that we need to scale now,
  // since the diagonal is not equal to
  // one now
  for (int row=N-1; row>=0; --row)
    {
      // get end of this row
      const size_type *const rowend = &column_numbers[rowstart_indices[row+1]];
      // find the position where the part
      // right of the diagonal starts
      const size_type *const first_after_diagonal = this->prebuilt_lower_bound[row];

      somenumber dst_row = dst(row);
      const number *luval = this->SparseMatrix<number>::val +
                            (first_after_diagonal - column_numbers);
      for (const size_type *col=first_after_diagonal; col!=rowend; ++col, ++luval)
        dst_row -= *luval * dst(*col);

      // scale by the diagonal element.
      // note that the diagonal element
      // was stored inverted
      dst(row) = dst_row * this->diag_element(row);
    }
}


template <typename number>
template <typename somenumber>
void SparseILU<number>::Tvmult (Vector<somenumber>       &dst,
                                const Vector<somenumber> &src) const
{
  Assert (dst.size() == src.size(), ExcDimensionMismatch(dst.size(), src.size()));
  Assert (dst.size() == this->m(), ExcDimensionMismatch(dst.size(), this->m()));

  const size_type N=dst.size();
  const std::size_t *const rowstart_indices
    = this->get_sparsity_pattern().rowstart;
  const size_type *const column_numbers
    = this->get_sparsity_pattern().colnums;

  // solve (LU)'x=b in two steps:
  // first U'y = b, then
  //       L'x = y
  //
  // first a forward solve. Due to the
  // fact that the transpose of U'
  // is not easily accessible, a
  // temporary vector is required.
  Vector<somenumber> tmp (N);

  dst = src;
  for (size_type row=0; row<N; ++row)
    {
      dst(row) -= tmp (row);
      // scale by the diagonal element.
      // note that the diagonal element
      // was stored inverted
      dst(row) *= this->diag_element(row);

      // get end of this row
      const size_type *const rowend = &column_numbers[rowstart_indices[row+1]];
      // find the position where the part
      // right of the diagonal starts
      const size_type *const first_after_diagonal = this->prebuilt_lower_bound[row];

      const somenumber dst_row = dst (row);
      const number *luval = this->SparseMatrix<number>::val +
                            (first_after_diagonal - column_numbers);
      for (const size_type *col=first_after_diagonal; col!=rowend; ++col, ++luval)
        tmp(*col) += *luval * dst_row;
    }

  // now the backward solve. same
  // procedure, but we need not set
  // dst before, since this is already
  // done.
  //
  // note that we no scaling is required
  // now, since the diagonal is one
  // now
  tmp = 0;
  for (int row=N-1; row>=0; --row)
    {
      dst(row) -= tmp (row);

      // get start of this row. skip the
      // diagonal element
      const size_type *const rowstart = &column_numbers[rowstart_indices[row]+1];
      // find the position where the part
      // right of the diagonal starts
      const size_type *const first_after_diagonal = this->prebuilt_lower_bound[row];

      const somenumber dst_row = dst (row);
      const number *luval = this->SparseMatrix<number>::val +
                            (rowstart - column_numbers);
      for (const size_type *col=rowstart; col!=first_after_diagonal; ++col, ++luval)
        tmp(*col) += *luval * dst_row;
    }
}


template <typename number>
std::size_t
SparseILU<number>::memory_consumption () const
{
  return SparseLUDecomposition<number>::memory_consumption ();
}



/*----------------------------   sparse_ilu.templates.h     ---------------------------*/

DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   sparse_ilu.templates.h     ---------------------------*/
