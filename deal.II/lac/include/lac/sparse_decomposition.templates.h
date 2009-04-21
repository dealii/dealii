//---------------------------------------------------------------------------
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2009 by the deal.II authors
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__sparse_decomposition_templates_h
#define __deal2__sparse_decomposition_templates_h

#include <base/memory_consumption.h>
#include <lac/sparse_decomposition.h>
#include <algorithm>

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
SparseLUDecomposition (const SparsityPattern& sparsity) :
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
  
  std::vector<const unsigned int*> tmp;
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
  
  if (data.use_this_sparsity)
    reinit(*data.use_this_sparsity);
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
      reinit(this->get_sparsity_pattern());
    }
  else if (data.extra_off_diagonals==0)
    {
				       // Use same sparsity as matrix
      reinit(matrix_sparsity);
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
      own_sparsity=new SparsityPattern(matrix_sparsity,
				       matrix_sparsity.max_entries_per_row()
				       +2*data.extra_off_diagonals,
				       data.extra_off_diagonals);
      own_sparsity->compress();
      reinit(*own_sparsity);
    }
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
void SparseLUDecomposition<number>::reinit (const SparsityPattern& sparsity)
{
  Assert (sparsity.optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());
  decomposed = false;
  if (true)
    {
      std::vector<const unsigned int*> tmp;
      tmp.swap (prebuilt_lower_bound);
    };
  SparseMatrix<number>::reinit (sparsity);
}



template<typename number>
void
SparseLUDecomposition<number>::prebuild_lower_bound()
{
  const unsigned int * const
    column_numbers = this->get_sparsity_pattern().get_column_numbers();
  const std::size_t * const
    rowstart_indices = this->get_sparsity_pattern().get_rowstart_indices();
  const unsigned int N = this->m();

  prebuilt_lower_bound.resize (N);

  for(unsigned int row=0; row<N; row++) {
    prebuilt_lower_bound[row]
      = std::lower_bound (&column_numbers[rowstart_indices[row]+1],
                          &column_numbers[rowstart_indices[row+1]],
                          row);
  }
}

template <typename number>
template <typename somenumber>
void
SparseLUDecomposition<number>::copy_from (const SparseMatrix<somenumber>& matrix)
{
				   // check whether we use the same sparsity
				   // pattern as the input matrix
  if (&this->get_sparsity_pattern() == &matrix.get_sparsity_pattern())
    {
      const somenumber * input_ptr = matrix.val;
      number * this_ptr = this->val;
      const number * const end_ptr = this_ptr + this->n_nonzero_elements();
      for ( ; this_ptr != end_ptr; ++input_ptr, ++this_ptr)
	*this_ptr = *input_ptr;
      return;
    }

                                   // preset the elements
  std::fill_n (&this->global_entry(0),
               this->n_nonzero_elements(),
               0);

                                   // note: pointers to the sparsity
                                   // pattern of the old matrix!
  const std::size_t * const rowstart_indices
    = matrix.get_sparsity_pattern().get_rowstart_indices();

  const unsigned int * const column_numbers
    = matrix.get_sparsity_pattern().get_column_numbers();

  for (unsigned int row=0; row<this->m(); ++row)
    for (const unsigned int * col = &column_numbers[rowstart_indices[row]];
         col != &column_numbers[rowstart_indices[row+1]]; ++col)
      this->set (row, *col, matrix.global_entry(col-column_numbers));
}



template <typename number>
void
SparseLUDecomposition<number>::strengthen_diagonal_impl ()
{
  for (unsigned int row=0; row<this->m(); ++row)
    {
                                       // get the length of the row
                                       // (without the diagonal element)
      const unsigned int rowlength
        = (this->get_sparsity_pattern().get_rowstart_indices()[row+1]
           -this->get_sparsity_pattern().get_rowstart_indices()[row]
           -1);
	
                                       // get the global index of the first
                                       // non-diagonal element in this row
      const unsigned int rowstart
        = this->get_sparsity_pattern().get_rowstart_indices()[row] + 1;
      number * const diagonal_element = &this->global_entry(rowstart-1);

      number rowsum = 0;
      for (unsigned int global_index=rowstart;
           global_index<rowstart+rowlength; ++global_index)
        rowsum += std::fabs(this->global_entry(global_index));

      *diagonal_element += this->get_strengthen_diagonal (rowsum, row)  *
                           rowsum;
    }
}



template <typename number>
unsigned int
SparseLUDecomposition<number>::memory_consumption () const
{
  unsigned int
    res = (SparseMatrix<number>::memory_consumption () +
           MemoryConsumption::memory_consumption(prebuilt_lower_bound));
  return res;
}


DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__sparse_decomposition_templates_h
