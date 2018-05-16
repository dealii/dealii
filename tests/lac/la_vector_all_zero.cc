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

// Check the function all_zero

#include "../tests.h"
#include <deal.II/lac/la_vector.h>


template <typename Number>
void
check_all_zero()
{
  LinearAlgebra::Vector<Number> v(10);

  AssertThrow(v.all_zero() == true, ExcInternalError());

  v[0] = 1.;
  AssertThrow(v.all_zero() == false, ExcInternalError());
}

int
main()
{
  initlog();

  check_all_zero<float>();
  check_all_zero<double>();

  deallog << "OK" << std::endl;
}
