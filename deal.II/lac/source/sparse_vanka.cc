/* $Id$ */

#include <lac/sparse_vanka.templates.h>


// explicit instantiations
template class SparseVanka<float>;
template class SparseVanka<double>;

template void SparseVanka<double>::operator () (Vector<float>       &dst,
						const Vector<float> &src) const;
template void SparseVanka<double>::operator () (Vector<double>       &dst,
						const Vector<double> &src) const;
