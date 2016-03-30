// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// test Tensor<1,dim> * SymmetricTensor<2,dim>

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<2,dim> s;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      s[i][j] = (i+1) * (j+1);

  Tensor<1,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    t[i] = (i==0 ? 1 : 2);

  deallog << s *t << std::endl
          << t *s << std::endl;

  Assert (s*t == t*s, ExcInternalError());
}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<3> ();
  test<5> ();

  deallog << "OK" << std::endl;
}
