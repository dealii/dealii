//----------------------------  vector.long_double.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vector.long_double.cc  ---------------------------


#include <lac/vector.templates.h>

// explicit instantiations
template class Vector<long double>;

template Vector<long double>::Vector (const Vector<double> &);
template Vector<long double>::Vector (const Vector<float> &);

template Vector<double>::Vector (const Vector<long double> &);
template Vector<float>::Vector (const Vector<long double> &);

template Vector<long double>& Vector<long double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<double>(const Vector<double>&);
template Vector<long double>& Vector<long double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<float>(const Vector<float>&);
template long double Vector<long double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<long double> (const Vector<long double> &) const;
template long double Vector<long double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<double> (const Vector<double> &) const;
template long double Vector<long double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<float> (const Vector<float> &) const;
template void Vector<long double>::reinit<long double>(const Vector<long double>&, const bool);
template void Vector<long double>::reinit<double>(const Vector<double>&, const bool);
template void Vector<long double>::reinit<float>(const Vector<float>&, const bool);
template void Vector<long double>::equ<long double>(const long double, const Vector<long double>&);
template void Vector<long double>::equ<double>(const long double, const Vector<double>&);
template void Vector<long double>::equ<float>(const long double, const Vector<float>&);

template Vector<double>& Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<long double>(const Vector<long double>&);
template double Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<long double> (const Vector<long double> &) const;
template void Vector<double>::reinit<long double>(const Vector<long double>&, const bool);
template void Vector<double>::equ<long double>(const double, const Vector<long double>&);

template Vector<float>& Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<long double>(const Vector<long double>&);
template float Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<long double> (const Vector<long double> &) const;
template void Vector<float>::reinit<long double>(const Vector<long double>&, const bool);
template void Vector<float>::equ<long double>(const float, const Vector<long double>&);

// see the .h file for why these functions are disabled.
// template Vector<float>::Vector (const Vector<double>& v);
// template Vector<double>::Vector (const Vector<float>& v);
