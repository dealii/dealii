// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
#  include "vector.inst"

#  ifndef DEAL_II_WITH_COMPLEX_VALUES
// instantiate for std::complex<double> since we are using it internally in
// FESeries.
template class Vector<std::complex<double>>;
#  endif

// instantiate for integers:
template class Vector<int>;
template Vector<double> &Vector<double>::
                         operator=<int>(const dealii::Vector<int> &);
template bool Vector<int>::operator==<int>(dealii::Vector<int> const &) const;

template void
Vector<int>::reinit<double>(const Vector<double> &, const bool);

// instantiate for long double manually because we use it in a few places:
template class Vector<long double>;
template long double Vector<long double>::
                     operator*<long double>(const Vector<long double> &) const;

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

#  define TEMPL_COPY_CONSTRUCTOR(S1, S2)             \
    template Vector<S1>::Vector(const Vector<S2> &); \
    template Vector<S1> &Vector<S1>::operator=<S2>(const Vector<S2> &)

TEMPL_COPY_CONSTRUCTOR(double, float);
TEMPL_COPY_CONSTRUCTOR(float, double);

#  ifdef DEAL_II_WITH_COMPLEX_VALUES
TEMPL_COPY_CONSTRUCTOR(std::complex<double>, std::complex<float>);
TEMPL_COPY_CONSTRUCTOR(std::complex<float>, std::complex<double>);
#  endif

#  undef TEMPL_COPY_CONSTRUCTOR


#  define TEMPL_OP_EQ(S1, S2)                            \
    template void Vector<S1>::scale(const Vector<S2> &); \
    template void Vector<S1>::equ(const S1, const Vector<S2> &)

TEMPL_OP_EQ(double, float);
TEMPL_OP_EQ(float, double);


#  ifdef DEAL_II_WITH_COMPLEX_VALUES
TEMPL_OP_EQ(std::complex<double>, std::complex<float>);
TEMPL_OP_EQ(std::complex<float>, std::complex<double>);
#  endif

#  undef TEMPL_OP_EQ
#endif


template <>
Vector<int>::real_type
Vector<int>::lp_norm(const real_type) const
{
  Assert(false, ExcMessage("No lp norm for integer vectors"));
  return -1;
}

DEAL_II_NAMESPACE_CLOSE
