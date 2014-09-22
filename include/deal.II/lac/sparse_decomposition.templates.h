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


#ifndef __deal2__sparse_decomposition_templates_h
#define __deal2__sparse_decomposition_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <algorithm>
#include <cstring>

DEAL_II_NAMESPACE_OPEN

template<typename number>
SparseLUDecomposition<number>::SparseLUDecomposition()
  :
  SparseMatrix<number>(),
  decomposed(false),
  own_sparsity(0)
{}



template<typename number>
SparseLUDecomposition<number>::
SparseLUDecomposition (const SparsityPattern &sparsity) :
  SparseMatrix<number>(sparsity),
  decomposed(false),
  own_sparsity(0)
{}



template<typename number>
SparseLUDecomposition<number>::~SparseLUDecomposition()
{
  clear();
}


template<typename number>
void SparseLUDecomposition<number>::clear()
{
  decomposed = false;

  std::vector<const size_type *> tmp;
  tmp.swap (prebuilt_lower_bound);

  SparseMatrix<number>::clear();

  if (own_sparsity)
    {
      delete own_sparsity;
      own_sparsity=0;
    }
}



template<typename number>
template <typename somenumber>
void SparseLUDecomposition<number>::initialize (
  const SparseMatrix<somenumber> &matrix,
  const AdditionalData data)
{
  const SparsityPattern &matrix_sparsity=matrix.get_sparsity_pattern();

  const SparsityPattern *sparsity_pattern_to_use = 0;

  if (data.use_this_sparsity)
    sparsity_pattern_to_use = data.use_this_sparsity;
  else if (data.use_previous_sparsity &&
           !this->empty() &&
           (this->m()==matrix.m()))
    {
      // Use the sparsity that was
      // previously used. This is
      // particularly useful when
      // matrix entries change but
      // not the sparsity, as for the
      // case of several Newton
      // iteration steps on an
      // unchanged grid.
      sparsity_pattern_to_use = &this->get_sparsity_pattern();
    }
  else if (data.extra_off_diagonals==0)
    {
      // Use same sparsity as matrix
      sparsity_pattern_to_use = &matrix_sparsity;
    }
  else
    {
      // Create new sparsity

      // for the case that
      // own_sparsity wasn't deleted
      // before (e.g. by clear()), do
      // it here
      if (own_sparsity)
        {
          // release the sparsity
          SparseMatrix<number>::clear();
          // delete it
          delete own_sparsity;
        }

      // and recreate
      own_sparsity = new SparsityPattern(matrix_sparsity,
                                         matrix_sparsity.max_entries_per_row()
                                         +2*data.extra_off_diagonals,
                                         data.extra_off_diagonals);
      own_sparsity->compress();
      sparsity_pattern_to_use = own_sparsity;
    }

  // now use this sparsity pattern
  Assert (sparsity_pattern_to_use->n_rows()==sparsity_pattern_to_use->n_cols(),
          typename SparsityPattern::ExcDiagonalNotOptimized());
  decomposed = false;
  {
    std::vector<const size_type *> tmp;
    tmp.swap (prebuilt_lower_bound);
  }
  SparseMatrix<number>::reinit (*sparsity_pattern_to_use);
}


template<typename number>
template<typename somenumber>
void
SparseLUDecomposition<number>::
decompose (const SparseMatrix<somenumber> &matrix,
           const double                    strengthen_diagonal)
{
  decomposed = false;

  this->strengthen_diagonal = strengthen_diagonal;
  prebuild_lower_bound ();
  copy_from (matrix);
  decomposed = true;
}



template <typename number>
void SparseLUDecomposition<number>::reinit (const SparsityPattern &sparsity)
{
  Assert (sparsity.n_rows() == sparsity.n_cols(),
          typename SparsityPattern::ExcDiagonalNotOptimized());
  decomposed = false;
  {
    std::vector<const size_type *> tmp;
    tmp.swap (prebuilt_lower_bound);
  }
  SparseMatrix<number>::reinit (sparsity);
}



template<typename number>
void
SparseLUDecomposition<number>::prebuild_lower_bound()
{
  const size_type *const
  column_numbers = this->get_sparsity_pattern().colnums;
  const std::size_t *const
  rowstart_indices = this->get_sparsity_pattern().rowstart;
  const size_type N = this->m();

  prebuilt_lower_bound.resize (N);

  for (size_type row=0; row<N; row++)
    {
      prebuilt_lower_bound[row]
        = Utilities::lower_bound (&column_numbers[rowstart_indices[row]+1],
                                  &column_numbers[rowstart_indices[row+1]],
                                  row);
    }
}

template <typename number>
template <typename somenumber>
void
SparseLUDecomposition<number>::copy_from (const SparseMatrix<somenumber> &matrix)
{
  // check whether we use the same sparsity
  // pattern as the input matrix
  if (&this->get_sparsity_pattern() == &matrix.get_sparsity_pattern())
    {
      const somenumber *input_ptr = matrix.val;
      number *this_ptr = this->val;
      const number *const end_ptr = this_ptr + this->n_nonzero_elements();
      if (types_are_equal<somenumber, number>::value == true)
        std::memcpy (this_ptr, input_ptr, this->n_nonzero_elements()*sizeof(number));
      else
        for ( ; this_ptr != end_ptr; ++input_ptr, ++this_ptr)
          *this_ptr = *input_ptr;
      return;
    }

  // preset the elements by zero. this needs to be written in a slightly
  // awkward way so that we find the corresponding function in the base class.
  SparseMatrix<number>::operator= (number(0));

  // both allow more and less entries in the new matrix
  for (size_type row=0; row<this->m(); ++row)
    {
      typename SparseMatrix<number>::iterator index = this->begin(row);
      typename SparseMatrix<somenumber>::const_iterator
      in_index = matrix.begin(row);
      index->value() = in_index->value();
      ++index, ++in_index;
      while (index < this->end(row) && in_index < matrix.end(row))
        {
          while (index->column() < in_index->column() && index < this->end(row))
            ++index;
          while (in_index->column() < index->column() && in_index < matrix.end(row))
            ++in_index;

          index->value() = in_index->value();
          ++index, ++in_index;
        }
    }
}



template <typename number>
void
SparseLUDecomposition<number>::strengthen_diagonal_impl ()
{
  for (size_type row=0; row<this->m(); ++row)
    {
      // get the global index of the first
      // non-diagonal element in this row
      Assert (this->m() == this->n(),  ExcNotImplemented());
      typename SparseMatrix<number>::iterator
      diagonal_element = this->begin(row);

      number rowsum = 0;
      for (typename SparseMatrix<number>::iterator
           p = diagonal_element + 1;
           p != this->end(row); ++p)
        rowsum += std::fabs(p->value());

      diagonal_element->value() += this->get_strengthen_diagonal (rowsum, row)  *
                                   rowsum;
    }
}



template <typename number>
std::size_t
SparseLUDecomposition<number>::memory_consumption () const
{
  return (SparseMatrix<number>::memory_consumption () +
          MemoryConsumption::memory_consumption(prebuilt_lower_bound));
}


DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__sparse_decomposition_templates_h
