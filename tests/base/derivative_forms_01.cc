// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/derivative_form.h>

#include "../tests.h"

#include "Sacado.hpp"

// Compute the derivative of F: R^dim -> R^spacedim
// We construct the function F_i, i=[0,spacedim)
// and we compute the first and second derivatives using
// Sacado

using Sdouble  = typename Sacado::Fad::DFad<double>;
using SSdouble = typename Sacado::Fad::DFad<Sdouble>;

template <int dim, int spacedim>
void
test()
{
  Tensor<1, dim, SSdouble> x;
  for (unsigned int j = 0; j < dim; ++j)
    {
      // This is for the first derivative
      Sdouble jv = numbers::PI * (j + 1) / 2.0;
      jv.diff(j, dim);
      // Now the second
      x[j] = jv;
      x[j].diff(j, dim);
    }

  Tensor<1, spacedim, SSdouble> F;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      F[i] += std::sin(x[j] + numbers::PI * i / 2);

  DerivativeForm<1, dim, spacedim, Sdouble> dF;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      dF[i][j] = F[i].dx(j);

  DerivativeForm<2, dim, spacedim, double> ddF;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        ddF[i][j][k] = dF[i][j].dx(k);

  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;
  deallog << "x  : " << x << std::endl;
  for (unsigned int i = 0; i < spacedim; ++i)
    deallog << "F[" << i << "] : " << F[i] << std::endl;
  for (unsigned int i = 0; i < spacedim; ++i)
    deallog << "dF[" << i << "] : " << dF[i] << std::endl;
  for (unsigned int i = 0; i < spacedim; ++i)
    deallog << "ddF[" << i << "] : " << ddF[i] << std::endl;
}

int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
