// ------------------------------------------------------------------------
//
// Copyright (C) 2023 - 2024 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Author: Tao Jin
 *         University of Ottawa, Ottawa, Ontario, Canada
 *         August 2023
 *
 * Test the positive-negative split of a 2nd-order symmetric tensor
 * and the positive and negative 4th-order projectors
 */


#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

#include "../tests.h"

template <int dim>
void
positive_negative_split_test()
{
  SymmetricTensor<2, dim> random_tensor;
  srand(time(0));

  for (unsigned int i = 0; i < dim; i++)
    for (unsigned int j = 0; j <= i; j++)
      {
        random_tensor[i][j] = ((double)rand() / (RAND_MAX));
        if (j != i)
          random_tensor[j][i] = random_tensor[i][j];
      }

  const auto &[positive_part_tensor, negative_part_tensor] =
    positive_negative_split(random_tensor);

  bool positive_negative_split_success = true;

  // test: (A^+) + (A^-) = A
  if ((positive_part_tensor + negative_part_tensor - random_tensor).norm() >
      1.0e-12 * random_tensor.norm())
    positive_negative_split_success = false;

  if (!positive_negative_split_success)
    Assert(false, ExcMessage("Positive-negative split failed!"));

  const auto &[positive_part_tensor_1,
               negative_part_tensor_1,
               positive_projector,
               negative_projector] =
    positive_negative_projectors(random_tensor);

  bool                    positive_projector_success = true;
  SymmetricTensor<2, dim> projected_positive_tensor;
  projected_positive_tensor = positive_projector * random_tensor;

  // test: (P^+) : A = (A^+)
  if ((projected_positive_tensor - positive_part_tensor_1).norm() >
      1.0e-12 * random_tensor.norm())
    positive_projector_success = false;

  // test: (P^+) : (A^+) = (A^+)
  if ((positive_projector * projected_positive_tensor - positive_part_tensor_1)
        .norm() > 1.0e-12 * random_tensor.norm())
    positive_projector_success = false;

  // test: (P^+) : (A^-) = 0
  if ((positive_projector * negative_part_tensor_1).norm() >
      1.0e-12 * random_tensor.norm())
    positive_projector_success = false;

  bool                    negative_projector_success = true;
  SymmetricTensor<2, dim> projected_negative_tensor;
  projected_negative_tensor = negative_projector * random_tensor;

  // test: (P^-) : A = (A^-)
  if ((projected_negative_tensor - negative_part_tensor_1).norm() >
      1.0e-12 * random_tensor.norm())
    negative_projector_success = false;

  // test: (P^-) : (A^-) = (A^-)
  if ((negative_projector * projected_negative_tensor - negative_part_tensor_1)
        .norm() > 1.0e-12 * random_tensor.norm())
    negative_projector_success = false;

  // test: (P^-) : (A^+) = 0
  if ((negative_projector * positive_part_tensor_1).norm() >
      1.0e-12 * random_tensor.norm())
    negative_projector_success = false;

  // test: (P^+) + (P^-) = S (S is 4th-order symmetric identity tensor)
  if ((positive_projector + negative_projector -
       Physics::Elasticity::StandardTensors<dim>::S)
        .norm() > 1.0e-12)
    {
      positive_projector_success = false;
      negative_projector_success = false;
    }

  if (!positive_projector_success)
    Assert(false, ExcMessage("Positive projector failed!"));

  if (!negative_projector_success)
    Assert(false, ExcMessage("Negative projector failed!"));
}


int
main()
{
  initlog();

  positive_negative_split_test<1>();
  positive_negative_split_test<2>();
  positive_negative_split_test<3>();

  deallog << "OK" << std::endl;

  return 0;
}
