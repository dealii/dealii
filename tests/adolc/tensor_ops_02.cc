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


// Check ability to cast from a tensor of ADOL-C numbers to floats
// In particular, this checks that the specializations for
// internal::NumberType<Number>::value() work as expected when
// Number == double and the input to value() is an ADOL-C number.

#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

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

  AD_Tensor adt1;
  for (unsigned int i = 0; i < NonAD_Tensor::n_independent_components; ++i)
    adt1[NonAD_Tensor::unrolled_to_component_indices(i)] = i + 1;

  const NonAD_Tensor t1 = NonAD_Tensor(adt1);

  Assert(t1.norm() > 0.0, ExcMessage("Cast and copy unsuccessful"));
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
  using AD_STensor2    = SymmetricTensor<2, dim, ad_number_t>;
  using AD_STensor4    = SymmetricTensor<4, dim, ad_number_t>;
  using NonAD_STensor2 = SymmetricTensor<2, dim, number_t>;
  using NonAD_STensor4 = SymmetricTensor<4, dim, number_t>;

  // Constructors
  AD_STensor2 adt1;
  AD_STensor4 adt2;

  for (unsigned int i = 0; i < NonAD_STensor2::n_independent_components; ++i)
    adt1[NonAD_STensor2::unrolled_to_component_indices(i)] = i + 1;

  for (unsigned int i = 0; i < NonAD_STensor4::n_independent_components; ++i)
    adt2.access_raw_entry(i) = i + 1;

  const NonAD_STensor2 t1 = NonAD_STensor2(adt1);
  const NonAD_STensor4 t2 = NonAD_STensor4(adt2);

  Assert(t1.norm() > 0.0, ExcMessage("Cast and copy unsuccessful"));
  Assert(t2.norm() > 0.0, ExcMessage("Cast and copy unsuccessful"));
}

int
main()
{
  initlog();

  // It is assumed that whatever works for rank-2
  // tensors will work for the other ranks

  // --- Taped ---
  test_tensor<3, double, Differentiation::AD::NumberTypes::adolc_taped>();
  test_symmetric_tensor<3,
                        double,
                        Differentiation::AD::NumberTypes::adolc_taped>();

  // --- Tapeless ---
  test_tensor<3, double, Differentiation::AD::NumberTypes::adolc_tapeless>();
  test_symmetric_tensor<3,
                        double,
                        Differentiation::AD::NumberTypes::adolc_tapeless>();

  deallog << "OK" << std::endl;
}
