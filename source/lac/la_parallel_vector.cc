// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "la_parallel_vector.inst"

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

namespace LinearAlgebra
{
  namespace distributed
  {
#define TEMPL_COPY_CONSTRUCTOR(S1, S2) \
  template Vector<S1>& Vector<S1>::operator=<S2>(const Vector<S2>&)

    TEMPL_COPY_CONSTRUCTOR(double, float);
    TEMPL_COPY_CONSTRUCTOR(float, double);

    TEMPL_COPY_CONSTRUCTOR(std::complex<double>, std::complex<float>);
    TEMPL_COPY_CONSTRUCTOR(std::complex<float>, std::complex<double>);

#undef TEMPL_COPY_CONSTRUCTOR
  } // namespace distributed
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE
