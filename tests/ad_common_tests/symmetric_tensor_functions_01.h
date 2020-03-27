// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Header file:
// Test to check that the various SymmetricTensor functions compile
// with the various auto-differentiable number types

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <iostream>

#include "../tests.h"

namespace AD = dealii::Differentiation::AD;

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_symmetric_tensor()
{
  typedef
    typename AD::NumberTraits<number_t, ad_type_code>::ad_type ADNumberType;

  std::cout << "*** Test SymmetricTensor functions, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  const ADNumberType                          a = 1.0;
  const Tensor<1, dim, ADNumberType>          v{};
  const SymmetricTensor<2, dim, ADNumberType> A(
    unit_symmetric_tensor<dim, ADNumberType>());
  const SymmetricTensor<2, dim, ADNumberType> B(
    unit_symmetric_tensor<dim, ADNumberType>());
  const Tensor<2, dim, ADNumberType> A_ns(
    unit_symmetric_tensor<dim, ADNumberType>());
  const SymmetricTensor<4, dim, ADNumberType> HH(
    identity_tensor<dim, ADNumberType>());

  const SymmetricTensor<2, dim, ADNumberType> C1 = A + B;
  const SymmetricTensor<2, dim, ADNumberType> C2 = A - B;
  const SymmetricTensor<2, dim, ADNumberType> C4 = a * A;
  const SymmetricTensor<2, dim, ADNumberType> C5 = A * a;
  const SymmetricTensor<2, dim, ADNumberType> C6 = A / a;

  const ADNumberType det_A = determinant(A);
  const ADNumberType tr_A  = trace(A);
  const ADNumberType I1_A  = first_invariant(A);
  const ADNumberType I2_A  = second_invariant(A);
  const ADNumberType I3_A  = third_invariant(A);

  const SymmetricTensor<2, dim, ADNumberType> A_inv  = invert(A);
  const SymmetricTensor<2, dim, ADNumberType> A_T    = transpose(A);
  const SymmetricTensor<2, dim, ADNumberType> A_dev  = deviator(A);
  const SymmetricTensor<2, dim, ADNumberType> A_symm = symmetrize(A_ns);
  // const ADNumberType A_l1_norm = l1_norm(A);
  // const ADNumberType A_linf_norm = linfty_norm(A);

  const ADNumberType                 A_ddot_B_1 = A * B;
  const ADNumberType                 sp_A_B     = scalar_product(A, B);
  const Tensor<4, dim, ADNumberType> op_A_B     = outer_product(A, B);


  SymmetricTensor<2, dim, ADNumberType> C_dc;
  double_contract(C_dc, A, HH);
  double_contract(C_dc, HH, A);

  const Tensor<1, dim, ADNumberType> v3 = A * v;
  const Tensor<1, dim, ADNumberType> v4 = v * A;
  const Tensor<2, dim, ADNumberType> C7 = A * A_ns;
  const Tensor<2, dim, ADNumberType> C8 = A_ns * A;
}
