// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2023 by the deal.II authors
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


// Verify that DerivativeForm can be initialized from a
// Tensor<1,dim,Tensor<order,spacedim>. We used to have a bug where the
// copy constructor offered Tensor<order,dim,Tensor<1,spacedim> instead,
// but that was just wrong.

#include <deal.II/base/derivative_form.h>

#include <boost/type_traits.hpp>

#include <complex>

#include "../tests.h"

template <typename Number>
void
test()
{
  // Choose orders and dimensions so that they are all not equal:
  constexpr int order    = 2;
  constexpr int dim      = 1;
  constexpr int spacedim = 3;


  Tensor<order, dim, Number> inner;
  inner[0][0] = 1;
  Tensor<1, spacedim, Tensor<order, dim, Number>> t;
  t[0] = t[1] = t[2] = inner;

  const DerivativeForm<order, dim, spacedim, Number> d(t);

  deallog << "norm: " << d.norm() << std::endl;
}

int
main()
{
  initlog();

  deallog << "testing float" << std::endl;
  test<float>();

  deallog << "testing double" << std::endl;
  test<double>();
}
