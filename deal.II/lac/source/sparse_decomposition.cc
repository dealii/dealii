//---------------------------------------------------------------------------
//    Copyright (C) 2002, 2003, 2005, 2006 by the deal.II authors
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/sparse_decomposition.templates.h>

DEAL_II_NAMESPACE_OPEN


template class SparseLUDecomposition<double>;
template void SparseLUDecomposition<double>::initialize<double> (const SparseMatrix<double> &,
								 const AdditionalData data);
template void SparseLUDecomposition<double>::initialize<float> (const SparseMatrix<float> &,
								 const AdditionalData data);
template void SparseLUDecomposition<double>::decompose<double> (const SparseMatrix<double> &,
                                                                const double);
template void SparseLUDecomposition<double>::decompose<float> (const SparseMatrix<float> &,
                                                               const double);

template void SparseLUDecomposition<double>::copy_from<double> (const SparseMatrix<double> &);
template void SparseLUDecomposition<double>::copy_from<float> (const SparseMatrix<float> &);


template class SparseLUDecomposition<float>;
template void SparseLUDecomposition<float>::initialize<double> (const SparseMatrix<double> &,
								const AdditionalData data);
template void SparseLUDecomposition<float>::initialize<float> (const SparseMatrix<float> &,
							       const AdditionalData data);
template void SparseLUDecomposition<float>::decompose<double> (const SparseMatrix<double> &,
                                                               const double);
template void SparseLUDecomposition<float>::decompose<float> (const SparseMatrix<float> &,
                                                              const double);

template void SparseLUDecomposition<float>::copy_from<double> (const SparseMatrix<double> &);
template void SparseLUDecomposition<float>::copy_from<float> (const SparseMatrix<float> &);

DEAL_II_NAMESPACE_CLOSE
