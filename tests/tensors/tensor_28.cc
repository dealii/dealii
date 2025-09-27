// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Additional tests for tensor contraction
//
// Specifically, its to clarify the order of arguments and tensor
// indices when performing contractions on the outermost indices
// of tensors. The test here implements the contravarient push-forward
// operation that is common in nonlinear continuum mechanics.
//
// The example case leading to the test involved contraction with the
// identity tensor: the result of having the arguments transposed
// was that some of the tensor elements became transposed.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"


template <int dim>
void
test_tensor(const Tensor<2, dim> &F)
{
  // Contraction with the a tensor on the outermost
  // indices of standard tensors:

  // Rank-2 Tensors
  {
    Tensor<2, dim> tmp1;
    unsigned int   c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        {
          tmp1[i][j] = c++;
        }

    const Tensor<2, dim> tmp2 = contract<1, 0>(
      F, contract<1, 1>(tmp1, F)); // Note: Order of arguments is important
    const Tensor<2, dim> tmp3 = F * tmp1 * transpose(F);

    Assert((tmp2 - tmp3).norm() < 1e-9,
           ExcMessage("Contraction using contract() function is incorrect."));
  }

  // Rank-3 Tensors
  {
    Tensor<3, dim> tmp1;
    unsigned int   c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          {
            tmp1[i][j][k] = c++;
          }

    const Tensor<3, dim> tmp2 = contract<1, 0>(
      F, contract<2, 1>(tmp1, F)); // Note: Order of arguments is important
    const Tensor<3, dim> tmp3 = F * tmp1 * transpose(F);

    Assert((tmp2 - tmp3).norm() < 1e-9,
           ExcMessage("Contraction using contract() function is incorrect."));
  }

  // Rank-4 Tensors
  {
    Tensor<4, dim> tmp1;
    unsigned int   c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            {
              tmp1[i][j][k][l] = c++;
            }

    const Tensor<4, dim> tmp2 = contract<1, 0>(
      F, contract<3, 1>(tmp1, F)); // Note: Order of arguments is important
    const Tensor<4, dim> tmp3 = F * tmp1 * transpose(F);

    Assert((tmp2 - tmp3).norm() < 1e-9,
           ExcMessage("Contraction using contract() function is incorrect."));
  }
}

template <int dim>
void
test()
{
  // Test with unit tensor
  test_tensor<dim>(unit_symmetric_tensor<dim>());

  // Test with non-trivial tensor
  Tensor<2, dim> F = unit_symmetric_tensor<dim>();
  double         c = 0.1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        F[i][j] += c;
        c += 0.05;
      }
  test_tensor<dim>(F);
}

int
main(int argc, char *argv[])
{
  initlog();

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;

  return 0;
}
