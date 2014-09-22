// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// a stress test using repeated multiplication of tensors

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>



template <int dim>
void test ()
{
  const double lambda = 1.5,
               mu     = 1.7;
  Tensor<2,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      t[i][j] = (1. + (i+lambda)*(mu+13));

  Tensor<2,dim> a;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      a[i][j] = (1. + (i+1)*(j+1));

  // stress test the whole thing many
  // times. normalize in each step to
  // make sure the result remains
  // representable in floating point
  // arithmetic. essentially, this
  // invokes the power method to
  // compute the largest eigenvector
  // (eigentensor in this case)
  for (unsigned int i=0; i<1000000; ++i)
    {
      a = t*a;
      a /= a.norm();
    }

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      deallog << i << ' ' << j << ' ' << a[i][j] << std::endl;
}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
