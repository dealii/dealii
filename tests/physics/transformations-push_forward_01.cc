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


// Check that the push-forward operation in the Physics::Transformations
// namespace works as expected for all tensor types.
// This gives particular scrutiny to the use of the contract() function
// of the tensor class.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/transformations.h>

#include "../tests.h"


template <int dim>
void
test_tensor(const Tensor<2, dim> &F)
{
  // Rank-1 Tensors
  {
    Tensor<1, dim> T;
    unsigned int   c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      T[i] = c++;

    // Hand calculation: Push forward index 0
    Tensor<1, dim> T_calc;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int I = 0; I < dim; ++I)
        T_calc[i] += F[i][I] * T[I];

    const Tensor<1, dim> T_trans =
      Physics::Transformations::Contravariant::push_forward(T, F);
    Assert(
      (T_calc - T_trans).norm() < 1e-9,
      ExcMessage(
        "Rank 1 tensor: Contraction using push_forward() function is incorrect."));
  }

  // Rank-2 Tensors
  {
    Tensor<2, dim> T;
    unsigned int   c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        T[i][j] = c++;

    // Hand calculation: Push forward index 0
    Tensor<2, dim> tmp0;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int J = 0; J < dim; ++J)
        for (unsigned int I = 0; I < dim; ++I)
          tmp0[i][J] += F[i][I] * T[I][J];

    // Hand calculation: Push forward index 1
    Tensor<2, dim> T_calc;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int J = 0; J < dim; ++J)
          T_calc[i][j] += F[j][J] * tmp0[i][J];

    const Tensor<2, dim> T_trans =
      Physics::Transformations::Contravariant::push_forward(T, F);
    Assert(
      (T_calc - T_trans).norm() < 1e-9,
      ExcMessage(
        "Rank 2 tensor: Contraction using push_forward() function is incorrect."));
  }

  // Rank-4 Tensors
  {
    Tensor<4, dim> T;
    unsigned int   c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            T[i][j][k][l] = c++;

    // Hand calculation: Push forward index 0
    Tensor<4, dim> tmp0;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int J = 0; J < dim; ++J)
        for (unsigned int K = 0; K < dim; ++K)
          for (unsigned int L = 0; L < dim; ++L)
            for (unsigned int I = 0; I < dim; ++I)
              tmp0[i][J][K][L] += F[i][I] * T[I][J][K][L];

    // Hand calculation: Push forward index 1
    Tensor<4, dim> tmp1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int K = 0; K < dim; ++K)
          for (unsigned int L = 0; L < dim; ++L)
            for (unsigned int J = 0; J < dim; ++J)
              tmp1[i][j][K][L] += F[j][J] * tmp0[i][J][K][L];

    // Hand calculation: Push forward index 2
    Tensor<4, dim> tmp2;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int L = 0; L < dim; ++L)
            for (unsigned int K = 0; K < dim; ++K)
              tmp2[i][j][k][L] += F[k][K] * tmp1[i][j][K][L];

    // Hand calculation: Push forward index 3
    Tensor<4, dim> T_calc;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            for (unsigned int L = 0; L < dim; ++L)
              T_calc[i][j][k][l] += F[l][L] * tmp2[i][j][k][L];

    const Tensor<4, dim> T_trans =
      Physics::Transformations::Contravariant::push_forward(T, F);
    Assert(
      (T_calc - T_trans).norm() < 1e-9,
      ExcMessage(
        "Rank 4 tensor: Contraction using push_forward() function is incorrect."));
  }
}

template <int dim>
void
test_symmetric_tensor(const Tensor<2, dim> &F)
{
  // Rank-2 Symmetric tensors
  {
    SymmetricTensor<2, dim> T;
    unsigned int            c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j) // symmetry
        T[i][j] = c++;

    // Hand calculation: Push forward index 0
    // Note: The symmetry of the initial tensor is lost after the
    //       transformation of the first index
    Tensor<2, dim> tmp0;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int J = 0; J < dim; ++J)
        for (unsigned int I = 0; I < dim; ++I)
          tmp0[i][J] += F[i][I] * T[I][J];

    // Hand calculation: Push forward index 1
    // Note: At this point we recover the symmetry of the
    //       transformed tensor.
    SymmetricTensor<2, dim> T_calc;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j) // symmetry
        for (unsigned int J = 0; J < dim; ++J)
          T_calc[i][j] += F[j][J] * tmp0[i][J];

    const Tensor<2, dim> T_trans =
      Physics::Transformations::Contravariant::push_forward(T, F);
    Assert(
      (T_calc - T_trans).norm() < 1e-9,
      ExcMessage(
        "Rank 2 symmetric tensor: Contraction using push_forward() function is incorrect."));
  }

  // Rank-4 Symmetric tensors
  {
    SymmetricTensor<4, dim> T;
    unsigned int            c = 1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j) // symmetry
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = k; l < dim; ++l)
            T[i][j][k][l] = c++;

    // Hand calculation: Push forward index 0
    // Note: The symmetry of the initial tensor is lost after the
    //       transformation of the first index
    Tensor<4, dim> tmp0;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int J = 0; J < dim; ++J)
        for (unsigned int K = 0; K < dim; ++K)
          for (unsigned int L = 0; L < dim; ++L)
            for (unsigned int I = 0; I < dim; ++I)
              tmp0[i][J][K][L] += F[i][I] * T[I][J][K][L];

    // Hand calculation: Push forward index 1
    Tensor<4, dim> tmp1;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int K = 0; K < dim; ++K)
          for (unsigned int L = 0; L < dim; ++L)
            for (unsigned int J = 0; J < dim; ++J)
              tmp1[i][j][K][L] += F[j][J] * tmp0[i][J][K][L];

    // Hand calculation: Push forward index 2
    Tensor<4, dim> tmp2;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int L = 0; L < dim; ++L)
            for (unsigned int K = 0; K < dim; ++K)
              tmp2[i][j][k][L] += F[k][K] * tmp1[i][j][K][L];

    // Hand calculation: Push forward index 3
    // Note: At this point we recover the symmetry of the
    //       transformed tensor.
    SymmetricTensor<4, dim> T_calc;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j) // symmetry
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = k; l < dim; ++l) // symmetry
            for (unsigned int L = 0; L < dim; ++L)
              T_calc[i][j][k][l] += F[l][L] * tmp2[i][j][k][L];

    const SymmetricTensor<4, dim> T_trans =
      Physics::Transformations::Contravariant::push_forward(T, F);

    std::cout << "T_calc: " << T_calc << std::endl;
    std::cout << "T_trans: " << T_trans << std::endl;

    Assert(
      (T_calc - T_trans).norm() < 1e-9,
      ExcMessage(
        "Rank 4 symmetric tensor: Contraction using push_forward() function is incorrect."));
  }
}

template <int dim>
void
test()
{
  // Test with unit tensor
  test_tensor<dim>(unit_symmetric_tensor<dim>());
  test_symmetric_tensor<dim>(unit_symmetric_tensor<dim>());

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
  test_symmetric_tensor<dim>(F);
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
