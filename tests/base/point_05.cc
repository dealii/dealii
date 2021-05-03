// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
      deallog << p_float(i) << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_double(i) << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_float(i) << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_double(i) << ' ';
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
      deallog << p_double(i) << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_float(i) << ' ';
    deallog << std::endl;

    for (unsigned int i = 0; i < dim; ++i)
      deallog << p_complex_double(i) << ' ';
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
