/* $Id$ */

#include <lac/sparse_ilu.templates.h>


// explicit instantiations
template class SparseILU<double>;
template void SparseILU<double>::decompose (const SparseMatrix<double> &,
					    const double);
template void SparseILU<double>::apply_decomposition (Vector<double> &,
						      const Vector<double> &) const;
template void SparseILU<double>::decompose (const SparseMatrix<float> &,
					    const double);
template void SparseILU<double>::apply_decomposition (Vector<float> &,
						      const Vector<float> &) const;


template class SparseILU<float>;
template void SparseILU<float>::decompose (const SparseMatrix<double> &,
					   const double);
template void SparseILU<float>::apply_decomposition (Vector<double> &,
						     const Vector<double> &) const;
template void SparseILU<float>::decompose (const SparseMatrix<float> &,
					   const double);
template void SparseILU<float>::apply_decomposition (Vector<float> &,
						     const Vector<float> &) const;

