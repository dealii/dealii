//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector.templates.h>

DEAL_II_NAMESPACE_OPEN

// explicit instantiations for long double
template class Vector<long double>;

#if !defined(DEAL_II_EXPLICIT_CONSTRUCTOR_BUG) && !defined(DEAL_II_LONG_DOUBLE_LOOP_BUG)
template Vector<long double>::Vector (const Vector<double> &);
template Vector<long double>::Vector (const Vector<float> &);

template Vector<double>::Vector (const Vector<long double> &);
template Vector<float>::Vector (const Vector<long double> &);
#endif

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


// explicit instantiations for std::complex<long double>
template class Vector<std::complex<long double> >;

#if !defined(DEAL_II_EXPLICIT_CONSTRUCTOR_BUG) && !defined(DEAL_II_LONG_DOUBLE_LOOP_BUG)
template Vector<std::complex<long double> >::Vector (const Vector<std::complex<double> > &);
template Vector<std::complex<long double> >::Vector (const Vector<std::complex<float> > &);

template Vector<std::complex<double> >::Vector (const Vector<std::complex<long double> > &);
template Vector<std::complex<float> >::Vector (const Vector<std::complex<long double> > &);
#endif

template Vector<std::complex<long double> >& Vector<std::complex<long double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<std::complex<double> >(const Vector<std::complex<double> >&);
template Vector<std::complex<long double> >& Vector<std::complex<long double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<std::complex<float> >(const Vector<std::complex<float> >&);
template std::complex<long double> Vector<std::complex<long double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<std::complex<long double> > (const Vector<std::complex<long double> > &) const;
template std::complex<long double> Vector<std::complex<long double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<std::complex<double> > (const Vector<std::complex<double> > &) const;
template std::complex<long double> Vector<std::complex<long double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<std::complex<float> > (const Vector<std::complex<float> > &) const;
template void Vector<std::complex<long double> >::reinit<std::complex<long double> >(const Vector<std::complex<long double> >&, const bool);
template void Vector<std::complex<long double> >::reinit<std::complex<double> >(const Vector<std::complex<double> >&, const bool);
template void Vector<std::complex<long double> >::reinit<std::complex<float> >(const Vector<std::complex<float> >&, const bool);
template void Vector<std::complex<long double> >::equ<std::complex<long double> >(const std::complex<long double>, const Vector<std::complex<long double> >&);
template void Vector<std::complex<long double> >::equ<std::complex<double> >(const std::complex<long double>, const Vector<std::complex<double> >&);
template void Vector<std::complex<long double> >::equ<std::complex<float> >(const std::complex<long double>, const Vector<std::complex<float> >&);

template Vector<std::complex<double> >& Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<std::complex<long double> >(const Vector<std::complex<long double> >&);
template std::complex<double> Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<std::complex<long double> > (const Vector<std::complex<long double> > &) const;
template void Vector<std::complex<double> >::reinit<std::complex<long double> >(const Vector<std::complex<long double> >&, const bool);
template void Vector<std::complex<double> >::equ<std::complex<long double> >(const std::complex<double>, const Vector<std::complex<long double> >&);

template Vector<std::complex<float> >& Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<std::complex<long double> >(const Vector<std::complex<long double> >&);
template std::complex<float> Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator *<std::complex<long double> > (const Vector<std::complex<long double> > &) const;
template void Vector<std::complex<float> >::reinit<std::complex<long double> >(const Vector<std::complex<long double> >&, const bool);
template void Vector<std::complex<float> >::equ<std::complex<long double> >(const std::complex<float>, const Vector<std::complex<long double> >&);

DEAL_II_NAMESPACE_CLOSE
