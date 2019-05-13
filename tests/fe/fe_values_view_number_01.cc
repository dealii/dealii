// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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
