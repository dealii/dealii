// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that mixed floating point type multiplications with std::complex
// are possible.

#include <deal.II/base/template_constraints.h>

#include "../tests.h"

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
