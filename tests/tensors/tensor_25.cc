// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check Tensor<2,dim>::operator[] (TableIndices)

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  using T = Tensor<3, dim>;
  T t;
  for (unsigned int i = 0; i < T::n_independent_components; ++i)
    t[T::unrolled_to_component_indices(i)] = (i + 1) * (i + 2);

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        {
          AssertThrow(TableIndices<3>(i, j, k) ==
                        T::unrolled_to_component_indices(i * dim * dim +
                                                         j * dim + k),
                      ExcInternalError());
          AssertThrow(T::component_to_unrolled_index(
                        TableIndices<3>(i, j, k)) ==
                        i * dim * dim + j * dim + k,
                      ExcInternalError());
          AssertThrow(t[TableIndices<3>(i, j, k)] == t[i][j][k],
                      ExcInternalError());
        }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<1>();
  check<2>();
  check<3>();

  deallog << "OK" << std::endl;
}
