// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/read_write_vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "read_write_vector.inst"

namespace LinearAlgebra
{
// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

#define TEMPL_COPY_CONSTRUCTOR(S1,S2)                           \
  template ReadWriteVector<S1>& ReadWriteVector<S1>::operator=<S2> (const ReadWriteVector<S2> &)

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
  TEMPL_COPY_CONSTRUCTOR(double,float);
  TEMPL_COPY_CONSTRUCTOR(float,double);


  TEMPL_COPY_CONSTRUCTOR(std::complex<double>,std::complex<float>);
  TEMPL_COPY_CONSTRUCTOR(std::complex<float>,std::complex<double>);

#endif

#undef TEMPL_COPY_CONSTRUCTOR
}

DEAL_II_NAMESPACE_CLOSE
