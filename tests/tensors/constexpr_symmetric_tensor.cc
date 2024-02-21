// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// create and manipulate constexpr SymmetricTensor objects

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int dim, typename Number>
void
test_symmetric_tensor()
{
  deallog << "*** Test constexpr SymmetricTensor functions, "
          << "dim = " << Utilities::to_string(dim) << std::endl;

  constexpr Number                 a = 1.0;
  constexpr Tensor<1, dim, Number> v{};
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> A(
    unit_symmetric_tensor<dim, Number>());
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> B(
    unit_symmetric_tensor<dim, Number>());
  DEAL_II_CONSTEXPR const Tensor<2, dim, Number> A_ns(
    unit_symmetric_tensor<dim, Number>());
  DEAL_II_CONSTEXPR const SymmetricTensor<4, dim, Number> HH(
    identity_tensor<dim, Number>());

  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> C1 = A + B;
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> C2 = A - B;
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> C4 = a * A;
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> C5 = A * a;
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> C6 = A / a;

  DEAL_II_CONSTEXPR const Number det_A = determinant(A);
  DEAL_II_CONSTEXPR const Number tr_A  = trace(A);
  DEAL_II_CONSTEXPR const Number I1_A  = first_invariant(A);
  DEAL_II_CONSTEXPR const Number I2_A  = second_invariant(A);
  DEAL_II_CONSTEXPR const Number I3_A  = third_invariant(A);

  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> A_inv = invert(A);
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> A_T   = transpose(A);
  DEAL_II_CONSTEXPR const SymmetricTensor<2, dim, Number> A_dev = deviator(A);

  DEAL_II_CONSTEXPR const Number A_ddot_B_1             = A * B;
  DEAL_II_CONSTEXPR const Number sp_A_B                 = scalar_product(A, B);
  DEAL_II_CONSTEXPR const Tensor<4, dim, Number> op_A_B = outer_product(A, B);

  DEAL_II_CONSTEXPR const Tensor<1, dim, Number> v3 = A * v;
  DEAL_II_CONSTEXPR const Tensor<1, dim, Number> v4 = v * A;
  DEAL_II_CONSTEXPR const Tensor<2, dim, Number> C7 = A * A_ns;
  DEAL_II_CONSTEXPR const Tensor<2, dim, Number> C8 = A_ns * A;
}

DEAL_II_CONSTEXPR SymmetricTensor<2, 2>
                  get_tensor_2()
{
  SymmetricTensor<2, 2> A;
  A[0][0] = 1.;
  A[1][1] = 3.;
  A[0][1] = -5.;
  return A;
}


DEAL_II_CONSTEXPR SymmetricTensor<4, 2>
                  get_tensor_4()
{
  SymmetricTensor<4, 2> B;
  B[0][0][0][0] = 1.;
  B[1][1][1][1] = 2.5;
  B[0][1][0][1] = 0.2;
  return B;
}


int
main()
{
  initlog();

  {
    LogStream::Prefix p("float");
    test_symmetric_tensor<1, float>();
    test_symmetric_tensor<2, float>();
    test_symmetric_tensor<3, float>();
  }

  {
    LogStream::Prefix p("double");
    test_symmetric_tensor<1, double>();
    test_symmetric_tensor<2, double>();
    test_symmetric_tensor<3, double>();
  }

  DEAL_II_CONSTEXPR const auto A = get_tensor_2();
  deallog << "SymmetricTensor<2,2> = " << A << std::endl;

  DEAL_II_CONSTEXPR const auto B = get_tensor_4();
  deallog << "SymmetricTensor<4,2> = " << B << std::endl;

  {
    constexpr double a_init[3][3] = {{1., 0., 0.}, {2., 1., 0.}, {3., 2., 1.}};
    constexpr Tensor<2, 3>  dummy_a{a_init};
    DEAL_II_CONSTEXPR const SymmetricTensor<2, 3> a = symmetrize(dummy_a);
    DEAL_II_CONSTEXPR const auto                  inverted = invert(a);
    constexpr double        ref_init[3][3]                 = {{0., -2., 2.},
                                                              {-2., 5., -2.},
                                                              {2., -2., 0.}};
    constexpr Tensor<2, 3>  dummy_ref{ref_init};
    DEAL_II_CONSTEXPR const SymmetricTensor<2, 3> ref = symmetrize(dummy_ref);
    Assert(inverted == ref, ExcInternalError());
  }
  {
    constexpr double a_init[3][3] = {{1., 2., 3.}, {2., 1., 2.}, {3., 2., 1.}};
    constexpr Tensor<2, 3>  dummy_a{a_init};
    DEAL_II_CONSTEXPR const SymmetricTensor<2, 3> a = symmetrize(dummy_a);
    DEAL_II_CONSTEXPR const auto                  transposed = transpose(a);
    Assert(transposed == a, ExcInternalError());
    DEAL_II_CONSTEXPR const auto dummy   = scalar_product(a, a);
    DEAL_II_CONSTEXPR const auto dummy_5 = a * a;
    DEAL_II_CONSTEXPR const auto middle  = outer_product(a, a);
    DEAL_II_CONSTEXPR const auto dummy_6 = contract3(a, middle, a);
  }

  deallog << "OK" << std::endl;
  return 0;
}
