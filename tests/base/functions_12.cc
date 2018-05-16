// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
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


// Test Polynomial

#include "../tests.h"
#include <deal.II/base/function_lib.h>

template <int dim>
void
check()
{
  unsigned int n_mon = 3;

  Table<2,double> exponents(n_mon,dim);

  for (unsigned int i = 0; i < n_mon; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      exponents[i][d] = i + d;

  std::vector<double> coeffs(n_mon);
  for (unsigned int i = 0; i < n_mon; ++i)
    coeffs[i] = std::pow(-1.0,static_cast<double>(i))*(i+1);

  Functions::Polynomial<dim> poly(exponents, coeffs);

  Point<dim> p;
  for (unsigned int d=0; d<dim; ++d)
    p[d] = d;

  deallog << dim << "-D check" << std::endl;
  deallog << "Polynomial: ";
  for (unsigned int i = 0; i < n_mon; ++i)
    {
      deallog << coeffs[i];
      for (unsigned int d = 0; d < dim; ++d)
        deallog << " x" << d << "^" << exponents[i][d];
      if (i < n_mon-1)
        deallog << " + ";
    }
  deallog << std::endl;
  deallog << "Point: ";
  for (unsigned int d = 0; d < dim; ++d)
    deallog << p[d] << " ";
  deallog << std::endl;

  deallog << "Value: " << poly.value(p) << std::endl;
  deallog << " values checked" << std::endl;

  deallog << "Gradient: " << poly.gradient(p) << std::endl;
  deallog << " gradients checked" << std::endl;
  deallog << std::endl;

}

int
main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);

  check<1>();
  check<2>();
  check<3>();

}



