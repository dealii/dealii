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

// explicit instantiations
template class Vector<double>;

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
template Vector<double>::Vector (const Vector<float> &);
#endif

template Vector<double>& Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<float>(const Vector<float>&);
template bool Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<double>(const Vector<double>&) const;
template bool Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<float>(const Vector<float>&) const;
template double Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<float>(const Vector<float>&) const;
template double Vector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<double>(const Vector<double>&) const;
template void Vector<double>::reinit<double>(const Vector<double>&, const bool);
template void Vector<double>::reinit<float>(const Vector<float>&, const bool);
template void Vector<double>::equ<double>(const double, const Vector<double>&);
template void Vector<double>::equ<float>(const double, const Vector<float>&);
template void Vector<double>::scale<double>(const Vector<double>&);
template void Vector<double>::scale<float>(const Vector<float>&);

template class Vector<float>;

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
template Vector<float>::Vector (const Vector<double> &);
#endif

template Vector<float>& Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<double>(const Vector<double>&);
template bool Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<double>(const Vector<double>&) const;
template bool Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<float>(const Vector<float>&) const;
template float Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<float>(const Vector<float>&) const;
template float Vector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<double>(const Vector<double>&) const;
template void Vector<float>::reinit<double>(const Vector<double>&, const bool);
template void Vector<float>::reinit<float>(const Vector<float>&, const bool);
template void Vector<float>::equ<double>(const float, const Vector<double>&);
template void Vector<float>::equ<float>(const float, const Vector<float>&);
template void Vector<float>::scale<double>(const Vector<double>&);
template void Vector<float>::scale<float>(const Vector<float>&);



// complex data types
template class Vector<std::complex<double> >;

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
template Vector<std::complex<double> >::Vector (const Vector<std::complex<float> > &);
#endif

template Vector<std::complex<double> >& Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<std::complex<float> >(const Vector<std::complex<float> >&);
template bool Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<std::complex<double> >(const Vector<std::complex<double> >&) const;
template bool Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<std::complex<float> >(const Vector<std::complex<float> >&) const;
template std::complex<double> Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<std::complex<float> >(const Vector<std::complex<float> >&) const;
template std::complex<double> Vector<std::complex<double> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<std::complex<double> >(const Vector<std::complex<double> >&) const;
template void Vector<std::complex<double> >::reinit<std::complex<double> >(const Vector<std::complex<double> >&, const bool);
template void Vector<std::complex<double> >::reinit<std::complex<float> >(const Vector<std::complex<float> >&, const bool);
template void Vector<std::complex<double> >::equ<std::complex<double> >(const std::complex<double>, const Vector<std::complex<double> >&);
template void Vector<std::complex<double> >::equ<std::complex<float> >(const std::complex<double>, const Vector<std::complex<float> >&);
template void Vector<std::complex<double> >::scale<std::complex<double> >(const Vector<std::complex<double> >&);
template void Vector<std::complex<double> >::scale<std::complex<float> >(const Vector<std::complex<float> >&);

template class Vector<std::complex<float> >;

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
template Vector<std::complex<float> >::Vector (const Vector<std::complex<double> > &);
#endif

template Vector<std::complex<float> >& Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<std::complex<double> >(const Vector<std::complex<double> >&);
template bool Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<std::complex<double> >(const Vector<std::complex<double> >&) const;
template bool Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator==<std::complex<float> >(const Vector<std::complex<float> >&) const;
template std::complex<float> Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<std::complex<float> >(const Vector<std::complex<float> >&) const;
template std::complex<float> Vector<std::complex<float> >::DEAL_II_MEMBER_OP_TEMPLATE_INST operator*<std::complex<double> >(const Vector<std::complex<double> >&) const;
template void Vector<std::complex<float> >::reinit<std::complex<double> >(const Vector<std::complex<double> >&, const bool);
template void Vector<std::complex<float> >::reinit<std::complex<float> >(const Vector<std::complex<float> >&, const bool);
template void Vector<std::complex<float> >::equ<std::complex<double> >(const std::complex<float>, const Vector<std::complex<double> >&);
template void Vector<std::complex<float> >::equ<std::complex<float> >(const std::complex<float>, const Vector<std::complex<float> >&);
template void Vector<std::complex<float> >::scale<std::complex<double> >(const Vector<std::complex<double> >&);
template void Vector<std::complex<float> >::scale<std::complex<float> >(const Vector<std::complex<float> >&);

DEAL_II_NAMESPACE_CLOSE
