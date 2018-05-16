// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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


// test construction, indexing, and conversion of SymmetricTensor<2,dim>
// from/to Tensor<2,dim> with dim > 3

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>


template <int dim>
void
test ()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<2,dim> s;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      s[i][j] = (i+1) * (j+1);

  Tensor<2,dim> t = s;
  SymmetricTensor<2,dim> u (t);

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      deallog << (i+1)*(j+1) << " " << (int)s[i][j] << " "
              << (int)t[i][j] << " " << (int)u[i][j] << std::endl;
}




int
main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  test<3> ();
  test<5> ();

  deallog << "OK" << std::endl;
}
