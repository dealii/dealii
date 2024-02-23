// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
