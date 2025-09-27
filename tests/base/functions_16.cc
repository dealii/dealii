// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test VectorFunctionFromTensorFunctionObject

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
check1()
{
  const auto x = [](const Point<dim> &p) -> Tensor<1, dim> { return p; };

  VectorFunctionFromTensorFunctionObject<dim> object(x, 1, dim + 2);

  AssertThrow(object.n_components == dim + 2, ExcInternalError());

  for (unsigned int i = 0; i < 10; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = i + d;

      for (unsigned int c = 0; c < dim + 2; ++c)
        if (c == 0 || c == dim + 1)
          {
            AssertThrow(object.value(p, c) == 0, ExcInternalError());
          }
        else
          {
            AssertThrow(object.value(p, c) == p[c - 1], ExcInternalError());
          }

      Vector<double> v(dim + 2);
      object.vector_value(p, v);
      for (unsigned int c = 0; c < dim + 2; ++c)
        if (c == 0 || c == dim + 1)
          {
            AssertThrow(v(c) == 0, ExcInternalError());
          }
        else
          {
            AssertThrow(v(c) == p[c - 1], ExcInternalError());
          }
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  check1<1>();
  check1<2>();
  check1<3>();
}
