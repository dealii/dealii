// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


// check serialization for Tensor<1,dim>

#include <deal.II/base/tensor.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>

#include "../tests.h"

#include "serialization.h"



void
test()
{
  const unsigned int dim  = 3;
  const unsigned int rank = 2;

  double            a1[3][3] = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
  Tensor<rank, dim> t1(a1);

  double a2[3][3] = {{10., 11., 12.}, {13., 14., 15.}, {16., 17., 18.}};
  Tensor<rank, dim> t2(a2);

  verify(t1, t2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
