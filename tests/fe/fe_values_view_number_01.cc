// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we get the right type of CurlType with Custom number type.

#include <deal.II/fe/fe_values.h>

#include <fstream>

#include "../tests.h"



template <typename Number>
void
test()
{
  if (typeid(typename internal::CurlType<1, Number>::type) !=
      typeid(Tensor<1, 1, Number>))
    deallog << "NOT OK!" << std::endl;
  if (typeid(typename internal::CurlType<2, Number>::type) !=
      typeid(Tensor<1, 1, Number>))
    deallog << "NOT OK!" << std::endl;
  if (typeid(typename internal::CurlType<3, Number>::type) !=
      typeid(Tensor<1, 3, Number>))
    deallog << "NOT OK!" << std::endl;
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  test<double>();
  test<float>();
}
