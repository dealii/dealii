// $Id$

#include <lac/vector.templates.h>

// explicit instantiations
template class Vector<long double>;
template Vector<long double>& Vector<long double>::operator=(const Vector<float>&);
template Vector<long double>& Vector<long double>::operator=(const Vector<double>&);
template Vector<double>& Vector<double>::operator=(const Vector<long double>&);
template Vector<float>& Vector<float>::operator=(const Vector<long double>&);
// see the .h file for why these functions are disabled.
// template Vector<float>::Vector (const Vector<double>& v);
// template Vector<double>::Vector (const Vector<float>& v);
