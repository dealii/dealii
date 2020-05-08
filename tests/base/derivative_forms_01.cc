// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>

#include "../tests.h"

#include "Sacado.hpp"

// Compute the derivative of F: R^dim -> R^spacedim
// We construct the function F_i, i=[0,spacedim)
// and we compute the first and second derivatives using
// Sacado

typedef typename Sacado::Fad::DFad<double>  Sdouble;
typedef typename Sacado::Fad::DFad<Sdouble> SSdouble;

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
