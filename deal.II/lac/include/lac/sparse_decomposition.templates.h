//---------------------  sparse_decomposition.templates.h  ----------------
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------  sparse_decomposition.templates.h  ----------------

#include <base/memory_consumption.h>
#include <lac/sparse_decomposition.h>
#include <algorithm>


template<typename number>
SparseLUDecomposition<number>::SparseLUDecomposition()
                :
                SparseMatrix<number>(),
                decomposed(false)
{}



template<typename number>
SparseLUDecomposition<number>::
SparseLUDecomposition (const SparsityPattern& sparsity) :
                SparseMatrix<number>(sparsity),
                decomposed(false)
{}



template<typename number>
SparseLUDecomposition<number>::~SparseLUDecomposition()
{}



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
void
SparseLUDecomposition<number>::reinit ()
{
  decomposed = false;
  if (true)
    {
      std::vector<const unsigned int*> tmp;
      tmp.swap (prebuilt_lower_bound);
    };
  SparseMatrix<number>::reinit ();
}



template <typename number>
void SparseLUDecomposition<number>::reinit (const SparsityPattern& sparsity)
{
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
  const unsigned int * const
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
                                   // preset the elements
  std::fill_n (&this->global_entry(0),
               this->n_nonzero_elements(),
               0);

                                   // note: pointers to the sparsity
                                   // pattern of the old matrix!
  const unsigned int * const rowstart_indices
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


