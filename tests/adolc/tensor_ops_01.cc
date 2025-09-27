// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// compilation check constructor and operations acting on tensors of ADOL-C
// number types

#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"


template <int dim,
          typename number_t,
          enum Differentiation::AD::NumberTypes ad_number_enum>
void
test_tensor()
{
  using ad_number_t =
    typename Differentiation::AD::NumberTraits<number_t,
                                               ad_number_enum>::ad_type;
  using AD_Tensor    = Tensor<2, dim, ad_number_t>;
  using NonAD_Tensor = Tensor<2, dim, number_t>;

  // Constructors
  NonAD_Tensor t1;
  AD_Tensor    adt1, adt2;
  AD_Tensor    adt3(adt1);

  // Assignment
  for (unsigned int i = 0; i < dim; ++i)
    {
      adt1[i][i] += 1.0;
      for (unsigned int j = 0; j < dim; ++j)
        {
          t1[i][j] = dim * i + j;
          adt1[i][j] += 2 * dim * i + j + 1;
          adt2[i][j] += 3 * dim * i + j + 2;
        }
    }
  AD_Tensor adt4 = adt1;

  // Comparison
  const bool res_1 = (adt1 == adt2);
  const bool res_2 = (adt1 != adt2);

  // Scalar operations
  adt4 *= 2.0;
  adt4 /= 4.0;
  adt4[0][0] += 1.0;
  adt4[0][0] *= 3.0;

  // Vector space operations
  adt4 = 2.0 * adt1;
  adt4 += adt2;
  adt4 -= adt2;
  adt4 += t1;
  adt4 -= t1;
  const AD_Tensor adt5 = adt1 + adt2;
  const AD_Tensor adt6 = adt1 - adt2;

  // Contractions
  const ad_number_t ad_res1 = scalar_product(adt1, adt2);
  const ad_number_t ad_res2 = double_contract<0, 0, 1, 1>(adt1, adt2);
  // TODO: Not defined. Conflicting number types
  // ad_res1 = double_contract(adt1,t1);
  // TODO: Not defined. Conflicting number types
  // ad_res1 = double_contract(t1,adt2);
  AD_Tensor adt7 = adt1 * adt2;
  adt7           = t1 * adt2;
  adt7           = adt1 * t1;

  // Outer product
  const Tensor<1, dim, number_t>    v1{};
  const Tensor<1, dim, ad_number_t> av1{};
  const AD_Tensor                   adt8  = outer_product(v1, av1);
  const AD_Tensor                   adt9  = outer_product(av1, v1);
  const AD_Tensor                   adt10 = outer_product(av1, av1);

  // Special operations
  const ad_number_t det = determinant(adt1);
  const ad_number_t tr  = trace(adt1);
  //  const ad_number_t l1 = l1_norm(adt1);       // TODO: Defined with std::
  //  math operator const ad_number_t linf = linfty_norm(adt1); // TODO: Defined
  //  with std:: math operator
  const AD_Tensor inv   = invert(adt1);
  const AD_Tensor trans = transpose(adt1);
}

template <int dim,
          typename number_t,
          enum Differentiation::AD::NumberTypes ad_number_enum>
void
test_symmetric_tensor()
{
  using ad_number_t =
    typename Differentiation::AD::NumberTraits<number_t,
                                               ad_number_enum>::ad_type;
  using AD_STensor     = SymmetricTensor<2, dim, ad_number_t>;
  using AD_STensor4    = SymmetricTensor<4, dim, ad_number_t>;
  using NonAD_STensor  = SymmetricTensor<2, dim, number_t>;
  using NonAD_STensor4 = SymmetricTensor<4, dim, number_t>;

  // Constructors
  NonAD_STensor  t1;
  NonAD_STensor4 t2;
  AD_STensor     adt1, adt2;
  AD_STensor     adt3(adt1);

  // Assignment
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = i; j < dim; ++j)
        {
          t1[i][j]   = dim * i + j;
          adt1[i][j] = 2 * dim * i + j + 1;
          adt2[i][j] += 3 * dim * i + j + 2;
        }
    }
  AD_STensor adt4 = adt1;

  // Comparison
  bool res_1 = (adt1 == adt2);
  bool res_2 = (adt1 != adt2);

  // Scalar operations
  adt4 *= 2.0;
  adt4 /= 4.0;
  adt4 = 2.0 * adt1;
  adt4[0][0] += 1.0;
  adt4[0][0] -= 0.5;
  adt4[0][0] *= 3.0;
  adt4[0][0] /= 1.5;

  // Vector space operations
  adt4 += adt2;
  adt4 -= adt2;
  adt4 += t1;
  adt4 -= t1;
  const AD_STensor adt5 = adt1 + adt2;
  const AD_STensor adt6 = adt1 - adt2;

  // Contractions
  const ad_number_t ad_res1 = scalar_product(adt1, adt2);
  const ad_number_t ad_res2 = adt1 * adt2;
  //  const ad_number_t ad_res3 = scalar_product(t1,adt2);
  //  const ad_number_t ad_res4 = scalar_product(adt1,t1);
  // const ad_number_t ad_res5 = t1*adt2; // TODO: Overload issue
  //  const ad_number_t ad_res6 = adt1*t1; // TODO: Overload issue

  // Outer product
  const Tensor<1, dim, number_t>    v1{};
  const Tensor<1, dim, ad_number_t> av1{};
  const AD_STensor                  adt8  = symmetrize(outer_product(v1, av1));
  const AD_STensor                  adt9  = symmetrize(outer_product(av1, v1));
  const AD_STensor                  adt10 = symmetrize(outer_product(av1, av1));

  // Tensor outer product
  const AD_STensor4 adt11 = outer_product(adt1, adt2);
  // const AD_STensor4 adt12 = outer_product(t1,adt1); // TODO: Cannot deduce
  // return type

  // Special operations
  const ad_number_t det  = determinant(adt1);
  const ad_number_t tr   = trace(adt1);
  const ad_number_t inv1 = first_invariant(adt1);
  const ad_number_t inv2 = second_invariant(adt1);
  const ad_number_t inv3 = third_invariant(adt1); // TODO: Return type incorrect
  // const ad_number_t l1 = l1_norm(adt1); // TODO: Defined with std:: math
  // operator const ad_number_t linf = linfty_norm(adt1); // TODO: Defined with
  // std:: math operator
  const AD_STensor inv   = invert(adt1);
  const AD_STensor trans = transpose(adt1);
  const AD_STensor dev   = deviator(adt1);

  // Special tensors
  const SymmetricTensor<2, dim, ad_number_t> ust =
    unit_symmetric_tensor<dim, ad_number_t>();
  const SymmetricTensor<4, dim, ad_number_t> dt =
    deviator_tensor<dim, ad_number_t>();
  //  SymmetricTensor<4,dim,ad_number_t> idt =
  //  identity_tensor<4,dim,ad_number_t> (); // TODO: No number type
  const SymmetricTensor<4, dim, ad_number_t> op = outer_product(adt1, adt2);
  //  SymmetricTensor<4,dim,ad_number_t> inv_op = invert(op); // TODO: Not
  //  implemented!
}

int
main()
{
  initlog();

  // It is assumed that whatever works for rank-2
  // tensors will work for the other ranks

  // --- Taped ---
  // double type
  test_tensor<3, double, Differentiation::AD::NumberTypes::adolc_taped>();
  test_symmetric_tensor<3,
                        double,
                        Differentiation::AD::NumberTypes::adolc_taped>();
  // complex double
  //  test_tensor<3,std::complex<double>,Differentiation::AD::NumberTypes::adolc_taped>();
  //  test_symmetric_tensor<3,std::complex<double>,Differentiation::AD::NumberTypes::adolc_taped>();

  // --- Tapeless ---
  // double type
  test_tensor<3, double, Differentiation::AD::NumberTypes::adolc_tapeless>();
  test_symmetric_tensor<3,
                        double,
                        Differentiation::AD::NumberTypes::adolc_tapeless>();
  // complex double
  //  test_tensor<3,std::complex<double>,Differentiation::AD::NumberTypes::adolc_tapeless>();
  //  test_symmetric_tensor<3,std::complex<double>,Differentiation::AD::NumberTypes::adolc_tapeless>();

  deallog << "OK" << std::endl;
}
