// $Id$

#include <lac/vector.templates.h>

// explicit instantiations
template class Vector<double>;
template class Vector<float>;

// see the .h file for why these functions are disabled.
// template Vector<float>::Vector (const Vector<double>& v);
// template Vector<double>::Vector (const Vector<float>& v);
