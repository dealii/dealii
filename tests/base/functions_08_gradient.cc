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


// Test VectorFunctionFromTensorFunction gradient

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
class X : public TensorFunction<1, dim>
{
public:
  virtual Tensor<2, dim>
  gradient(const Point<dim> &p) const
  {
    Tensor<2, dim> T;
    for (unsigned int d = 0; d < dim; ++d)
      T[d] = d * p;
    return T;
  }
};


template <int dim>
void
check1()
{
  X<dim>                                x;
  VectorFunctionFromTensorFunction<dim> object(x, 1, dim + 2);

  AssertThrow(object.n_components == dim + 2, ExcInternalError());

  Tensor<1, dim> T;

  for (unsigned int i = 0; i < 10; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = i + d;

      for (unsigned int c = 0; c < dim + 2; ++c)
        if (c == 0 || c == dim + 1)
          {
            AssertThrow(object.gradient(p, c) == T, ExcInternalError());
          }
        else
          {
            AssertThrow(object.gradient(p, c) == (c - 1) * p,
                        ExcInternalError());
          }

      std::vector<Tensor<1, dim>> v(dim + 2);
      object.vector_gradient(p, v);
      for (unsigned int c = 0; c < dim + 2; ++c)
        if (c == 0 || c == dim + 1)
          {
            AssertThrow(v[c] == T, ExcInternalError());
          }
        else
          {
            AssertThrow(v[c] == (c - 1) * p, ExcInternalError());
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
