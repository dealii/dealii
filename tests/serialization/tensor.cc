// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
