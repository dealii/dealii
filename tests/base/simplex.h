// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// construct a simplex quadrature, and check that we can get an affine
// transformation out of it.
#ifndef tests_base_simplex_h
#define tests_base_simplex_h

#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>

#include "simplex.h"

#include <numeric>


// Helper functions
template<int dim>
std::array<Point<dim>, dim+1> get_simplex();

template<>
std::array<Point<1>, 2> get_simplex()
{
  return {{Point<1>(3), Point<1>(5)}};
}


template<>
std::array<Point<2>, 3> get_simplex()
{
  return {{Point<2>(4,2), Point<2>(3,3), Point<2>(2,2.5)}};
}


template<>
std::array<Point<3>, 4> get_simplex()
{
  return {{Point<3>(4,2,0), Point<3>(3,3,0), Point<3>(2,2.5,0), Point<3>(4.5, 3, 2)}};
}


#endif
