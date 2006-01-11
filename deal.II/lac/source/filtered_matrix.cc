//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/filtered_matrix.templates.h>

#define FILT(MM,VV) \
template <> \
void FilteredMatrix<MM , VV >::allocate_tmp_vector() \
{\
    Threads::ThreadMutex::ScopedLock lock (tmp_mutex);  \
    tmp_vector.reinit (matrix->n(), true);  \
}

#define BFILT(MM,VV) \
template <>   \
void   \
FilteredMatrix<MM ,VV >::   \
allocate_tmp_vector ()    \
{   \
  std::vector<unsigned int> block_sizes (matrix->n_block_rows());   \
  for (unsigned int i=0; i<block_sizes.size(); ++i)   \
    block_sizes[i] = matrix->block(i,i).n();   \
     \
  Threads::ThreadMutex::ScopedLock lock (tmp_mutex);   \
  tmp_vector.reinit (block_sizes, true);   \
}

FILT(SparseMatrix<double>, Vector<double>)
BFILT(BlockSparseMatrix<double>, BlockVector<double>)
template class FilteredMatrix<SparseMatrix<double>,Vector<double> >;
template class FilteredMatrix<BlockSparseMatrix<double>,BlockVector<double> >;

FILT(SparseMatrix<float>, Vector<float>)
BFILT(BlockSparseMatrix<float>, BlockVector<float>)
template class FilteredMatrix<SparseMatrix<float>,Vector<float> >;
template class FilteredMatrix<BlockSparseMatrix<float>,BlockVector<float> >;
