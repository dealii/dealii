// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for Point::operator= with different number types

#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  {
    Point<dim, float>                p_float;
    Point<dim, double>               p_double;
    Point<dim, std::complex<float>>  p_complex_float;
    Point<dim, std::complex<double>> p_complex_double;

    for (unsigned int i = 0; i < dim; ++i)
      p_float[i] = i;

    p_double         = p_float;
    p_complex_float  = p_float;
    p_complex_double = p_float;


    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_float[i] << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_double[i] << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_float[i] << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_double[i] << ' ';
    deallog << std::endl;
  }

  {
    Point<dim, double>               p_double;
    Point<dim, std::complex<float>>  p_complex_float;
    Point<dim, std::complex<double>> p_complex_double;

    for (unsigned int i = 0; i < dim; ++i)
      p_double[i] = i + 0.1;

    p_complex_float  = p_double;
    p_complex_double = p_double;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_double[i] << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_float[i] << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_double[i] << ' ';
    deallog << std::endl;
  }
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<1>();
  check<2>();
  check<3>();
}
