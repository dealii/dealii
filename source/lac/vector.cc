// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
#  include "lac/vector.inst"

#  ifndef DEAL_II_WITH_COMPLEX_VALUES
// instantiate for std::complex<double> since we are using it internally in
// FESeries.
template class Vector<std::complex<double>>;
#  endif

// instantiate for integers:
template class Vector<int>;
template Vector<double> &
Vector<double>::operator=<int>(const dealii::Vector<int> &);
template bool
Vector<int>::operator==<int>(const dealii::Vector<int> &) const;

// instantiate for long double manually because we use it in a few places:
template class Vector<long double>;
template long double
Vector<long double>::operator*<long double>(const Vector<long double> &) const;

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

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef HAVE_TPETRA_INST_FLOAT
template Vector<float>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Host> &);
template Vector<float>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default> &);
template Vector<float> &
Vector<float>::operator=<float>(
  const LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Host> &);
template Vector<float> &
Vector<float>::operator=<float>(
  const LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default> &);
#  endif

#  ifdef HAVE_TPETRA_INST_DOUBLE
template Vector<double>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host> &);
template Vector<double>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &);
template Vector<double> &
Vector<double>::operator=<double>(
  const LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host> &);
template Vector<double> &
Vector<double>::operator=<double>(
  const LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &);
#  endif

#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
template Vector<std::complex<float>>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<float>,
                                              MemorySpace::Host> &);
template Vector<std::complex<float>>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<float>,
                                              MemorySpace::Default> &);
template Vector<std::complex<float>> &
Vector<std::complex<float>>::operator=<std::complex<float>>(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<float>,
                                              MemorySpace::Host> &);
template Vector<std::complex<float>> &
Vector<std::complex<float>>::operator=<std::complex<float>>(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<float>,
                                              MemorySpace::Default> &);
#    endif

#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
template Vector<std::complex<double>>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<double>,
                                              MemorySpace::Host> &);
template Vector<std::complex<double>>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<double>,
                                              MemorySpace::Default> &);
template Vector<std::complex<double>> &
Vector<std::complex<double>>::operator=<std::complex<double>>(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<double>,
                                              MemorySpace::Host> &);
template Vector<std::complex<double>> &
Vector<std::complex<double>>::operator=<std::complex<double>>(
  const LinearAlgebra::TpetraWrappers::Vector<std::complex<double>,
                                              MemorySpace::Default> &);
#    endif
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE
