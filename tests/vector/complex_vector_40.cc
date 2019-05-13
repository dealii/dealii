// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// check Vector<std::complex<double> >::add(s,V,s,V)

#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"


void
test(Vector<std::complex<double>> &v,
     Vector<std::complex<double>> &w,
     Vector<std::complex<double>> &x)
{
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      v(i) = i;
      w(i) = std::complex<double>(i + 1., i + 2.);
      x(i) = i + 2.;
    }

  v.compress();
  w.compress();
  x.compress();

  v.add(2, w, 3, x);

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      AssertThrow(w(i) == std::complex<double>(i + 1., i + 2.),
                  ExcInternalError());
      AssertThrow(x(i) == i + 2., ExcInternalError());
      AssertThrow(v(i) == 1. * i + 2. * std::complex<double>(i + 1., i + 2.) +
                            3 * (i + 2.),
                  ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      Vector<std::complex<double>> v(100);
      Vector<std::complex<double>> w(100);
      Vector<std::complex<double>> x(100);
      test(v, w, x);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
