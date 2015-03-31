// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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



// check Vector<std::complex<double> >::sadd(s,s,V,s,V,s,V)

#include "../tests.h"
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>
#include <vector>


void test (Vector<std::complex<double> > &v,
           Vector<std::complex<double> > &w,
           Vector<std::complex<double> > &x,
           Vector<std::complex<double> > &y)
{
  for (unsigned int i=0; i<v.size(); ++i)
    {
      v(i) = i;
      w(i) = std::complex<double> (i+1., i+2.);
      x(i) = i+2.;
      y(i) = std::complex<double> (i+3., i+4.);
    }

  v.compress ();
  w.compress ();
  x.compress ();
  y.compress ();

  v.sadd (1.5, 2, w, 3, x, 4, y);

  // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      AssertThrow (w(i) == std::complex<double> (i+1., i+2.),
                   ExcInternalError());
      AssertThrow (x(i) == i+2., ExcInternalError());
      AssertThrow (y(i) == std::complex<double> (i+3., i+4.),
                   ExcInternalError());
      AssertThrow (v(i) ==
                   1.5*i+2.*std::complex<double> (i+1., i+2.)+
                   3.*(i+2.)+4.*std::complex<double> (i+3., i+4.),
                   ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Vector<std::complex<double> > v (100);
      Vector<std::complex<double> > w (100);
      Vector<std::complex<double> > x (100);
      Vector<std::complex<double> > y (100);
      test (v,w,x,y);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
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
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
