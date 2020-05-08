// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// check   Tensor<rank,dim,std::complex<double> > operator * (const
// Tensor<rank,dim>     &t,
//                                                            const
//                                                            std::complex<double>
//                                                            factor)
// and
// check   Tensor<rank,dim,std::complex<double> > operator * (const
// std::complex<double>  factor,
//                                                            const
//                                                            Tensor<rank,dim>
//                                                            &t)
// by multiplying on the left and on the right.

#include <deal.II/base/tensor.h>

#include "../tests.h"


template <int dim>
void
test_tensor()
{
  // a real tensor
  Tensor<1, dim, double> t;

  for (unsigned int i = 0; i < dim; ++i)
    {
      t[i] = 2 * i + dim + 1;
    }

  // multiply on the right by a complex<double>
  const Tensor<1, dim, std::complex<double>> right =
    t * std::complex<double>(1, 2);

  // multiply on the left by a complex<double>
  const Tensor<1, dim, std::complex<double>> left =
    std::complex<double>(1, 2) * t;

  // they should yield the same result
  Assert(left == right, ExcInternalError());

  deallog << "dim = " << dim << std::endl
          << left << " : " << right << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_tensor<1>();
  test_tensor<2>();
  test_tensor<3>();
  test_tensor<4>();

  deallog << "OK" << std::endl;
}
