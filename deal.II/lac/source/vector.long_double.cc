//----------------------------  std::vector.long_double.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  std::vector.long_double.cc  ---------------------------


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
