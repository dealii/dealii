// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#include <deal.II/lac/la_sm_vector.h>
#include <deal.II/lac/la_sm_vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "la_sm_vector.inst"

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

namespace LinearAlgebra
{
  namespace SharedMPI
  {
#define TEMPL_COPY_CONSTRUCTOR(S1, S2)                  \
  template Vector<S1, ::dealii::MemorySpace::Host>      \
    &Vector<S1, ::dealii::MemorySpace::Host>::operator= \
      <S2>(const Vector<S2, ::dealii::MemorySpace::Host> &)

    TEMPL_COPY_CONSTRUCTOR(double, float);
    TEMPL_COPY_CONSTRUCTOR(float, double);

#undef TEMPL_COPY_CONSTRUCTOR
  } // namespace SharedMPI
} // namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE
