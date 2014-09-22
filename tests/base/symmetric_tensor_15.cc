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


// test outer_product

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<2,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      t[i][j] = (1.+(i+1)*(j*2));

  // test 1: check trace-like operator
  {
    const SymmetricTensor<4,dim> T
      = outer_product<dim> (unit_symmetric_tensor<dim>(),
                            unit_symmetric_tensor<dim>());

    // T*t should yield a diagonal tensor
    // where the diagonal elements are the
    // traces of t
    SymmetricTensor<2,dim> x = T * t;
    Assert ((x-trace(t)*unit_symmetric_tensor<dim>()).norm()
            < 1e-15*t.norm(), ExcInternalError());

    deallog << "x=" << std::endl;
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        deallog << i << ' ' << j << ' ' << x[i][j] << std::endl;
  }

  // test 2: check outer product of t with
  // itself
  {
    const SymmetricTensor<4,dim> T = outer_product<dim> (t,t);

    // T*t should yield norm(t)^2*t
    SymmetricTensor<2,dim> x = T * t;
    Assert ((x-(t*t)*t).norm()
            < 1e-15*t.norm(), ExcInternalError());

    deallog << "x=" << std::endl;
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        deallog << i << ' ' << j << ' ' << x[i][j] << std::endl;
  }

}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
