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

// Test that mixed floating point type multiplications with std::complex
// are possible.

#include "../tests.h"

#include <deal.II/base/template_constraints.h>

int
main()
{
  initlog();

  double() * std::complex<float>();
  std::complex<float>() * double();

  float() * std::complex<double>();
  std::complex<double>() * float();

  std::complex<double>() * std::complex<float>();
  std::complex<float>() * std::complex<double>();

  deallog << "OK" << std::endl;
}
