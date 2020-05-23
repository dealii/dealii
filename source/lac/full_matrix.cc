// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
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

#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "full_matrix.inst"

#ifndef DOXYGEN

#  ifndef DEAL_II_WITH_COMPLEX_VALUES
// instantiate for std::complex<double> because we use it internally in
// FESeries.
template class FullMatrix<std::complex<double>>;
#  endif

// instantiate for long double manually because we use it in a few places
// inside the library
template class FullMatrix<long double>;
template void
FullMatrix<long double>::invert<long double>(const FullMatrix<long double> &);
template void
FullMatrix<long double>::mmult<long double>(FullMatrix<long double> &,
                                            const FullMatrix<long double> &,
                                            const bool) const;
template void
FullMatrix<long double>::Tmmult<long double>(FullMatrix<long double> &,
                                             const FullMatrix<long double> &,
                                             const bool) const;
template void
FullMatrix<long double>::mTmult<long double>(FullMatrix<long double> &,
                                             const FullMatrix<long double> &,
                                             const bool) const;
template void
FullMatrix<long double>::TmTmult<long double>(FullMatrix<long double> &,
                                              const FullMatrix<long double> &,
                                              const bool) const;
template void
FullMatrix<long double>::vmult<long double>(Vector<long double> &,
                                            const Vector<long double> &,
                                            bool) const;
template void
FullMatrix<long double>::Tvmult<long double>(Vector<long double> &,
                                             const Vector<long double> &,
                                             bool) const;
template void
FullMatrix<long double>::add<long double>(const long double,
                                          const FullMatrix<long double> &);


// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

#  define TEMPL_OP_EQ(S1, S2) \
    template FullMatrix<S1> &FullMatrix<S1>::operator=(const FullMatrix<S2> &)

TEMPL_OP_EQ(double, float);
TEMPL_OP_EQ(float, double);

#  ifdef DEAL_II_WITH_COMPLEX_VALUES
TEMPL_OP_EQ(std::complex<double>, std::complex<float>);
TEMPL_OP_EQ(std::complex<float>, std::complex<double>);
TEMPL_OP_EQ(std::complex<double>, double);
TEMPL_OP_EQ(std::complex<float>, float);
#  endif

#  undef TEMPL_OP_EQ

#endif

DEAL_II_NAMESPACE_CLOSE
