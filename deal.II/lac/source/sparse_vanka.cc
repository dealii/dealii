/* $Id$ */

#include <lac/sparse_vanka.templates.h>


// explicit instantiations
template class SparseVanka<float>;
template class SparseVanka<double>;

template void SparseVanka<double>::forward (Vector<float>       &dst,
					    const Vector<float> &src) const;
template void SparseVanka<double>::forward (Vector<double>       &dst,
					    const Vector<double> &src) const;
template void SparseVanka<float>::forward (Vector<float>       &dst,
					   const Vector<float> &src) const;
template void SparseVanka<float>::forward (Vector<double>       &dst,
					   const Vector<double> &src) const;
