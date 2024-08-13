// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Create and copy tensors in constexpr setting

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int rank, int dim, typename Number>
void
test_constexpr_tensor_constructors()
{
  constexpr dealii::Tensor<rank, dim, Number> a;
  constexpr dealii::Tensor<rank, dim, Number> b(a);
  constexpr dealii::Tensor<rank, dim, Number> c = a;
  deallog << " Tensor<" << rank << ',' << dim << '>' << std::endl;
  deallog << a << std::endl;
  deallog << b << std::endl;
  deallog << c << std::endl;

  DEAL_II_CONSTEXPR const auto d = a * 2.;
  DEAL_II_CONSTEXPR const auto e = 2. * a;
  DEAL_II_CONSTEXPR const auto f = a / 2.;
  DEAL_II_CONSTEXPR const auto g = d + e;
  DEAL_II_CONSTEXPR const auto h = d - e;
  deallog << d << std::endl;
  deallog << e << std::endl;
  deallog << f << std::endl;
  deallog << g << std::endl;
  deallog << h << std::endl;
}

DEAL_II_CONSTEXPR double
tensor_rank0_constexpr(const double value)
{
  const Tensor<0, 2> a{value};
  return 2. * a;
}

DEAL_II_CONSTEXPR double
tensor_rank2_constexpr()
{
  constexpr double a_init[3][3] = {{4., 2., 0.}, {2., 3., 0.}, {0., 0., 5.}};
  constexpr Tensor<2, 3>       a(a_init);
  DEAL_II_CONSTEXPR const auto det   = determinant(a);
  DEAL_II_CONSTEXPR const auto tr    = trace(transpose(a) * a);
  DEAL_II_CONSTEXPR const auto dummy = symmetrize(a);

  auto d = a;
  d /= 2.;
  const auto test_1 = d[0][0];
  d *= 2.;
  const auto test_2 = d[0][0];
  d += a;
  const auto test_3 = d[0][0];
  d -= a;
  const auto test_4 = d[0][0];
  d.clear();
  const auto test_5 = d[0][0];
  d                 = 0.;
  const auto test_6 = d[0][0];
  const auto test_7 = d.norm_square();

  return det + tr + test_1 + test_2 + test_3 + test_4 + test_5 + test_6 +
         test_7;
}

int
main()
{
  initlog();

  deallog << "Checking constexpr default constructor of Tensor<rank,dim,Number>"
          << std::endl;
  {
    dealii::LogStream::Prefix p("float");
    test_constexpr_tensor_constructors<0, 1, float>();
    test_constexpr_tensor_constructors<0, 2, float>();
    test_constexpr_tensor_constructors<0, 3, float>();
    test_constexpr_tensor_constructors<1, 1, float>();
    test_constexpr_tensor_constructors<1, 2, float>();
    test_constexpr_tensor_constructors<1, 3, float>();
    test_constexpr_tensor_constructors<2, 1, float>();
    test_constexpr_tensor_constructors<2, 2, float>();
    test_constexpr_tensor_constructors<2, 3, float>();
  }
  {
    dealii::LogStream::Prefix p("double");
    test_constexpr_tensor_constructors<0, 1, double>();
    test_constexpr_tensor_constructors<0, 2, double>();
    test_constexpr_tensor_constructors<0, 3, double>();
    test_constexpr_tensor_constructors<1, 1, double>();
    test_constexpr_tensor_constructors<1, 2, double>();
    test_constexpr_tensor_constructors<1, 3, double>();
    test_constexpr_tensor_constructors<2, 1, double>();
    test_constexpr_tensor_constructors<2, 2, double>();
    test_constexpr_tensor_constructors<2, 3, double>();
  }

  deallog << "Using Tensor within constexpr functions" << std::endl;
  constexpr double             number = 7.6;
  DEAL_II_CONSTEXPR const auto tensor_rank0_result =
    tensor_rank0_constexpr(number);
  deallog << tensor_rank0_result << std::endl;
  DEAL_II_CONSTEXPR const auto tensor_rank2_result = tensor_rank2_constexpr();
  deallog << tensor_rank2_result << std::endl;

  {
    constexpr double        initializer[2] = {1., -1.};
    constexpr Tensor<1, 2>  c{initializer};
    DEAL_II_CONSTEXPR const Tensor<1, 2> c_cross = cross_product_2d(c);
    constexpr double                     initializer_ref[2] = {-1., -1.};
    constexpr Tensor<1, 2>               ref{initializer_ref};
    DEAL_II_CONSTEXPR const bool         is_same     = (c_cross == ref);
    DEAL_II_CONSTEXPR const bool         is_not_same = (c_cross != ref);
    Assert(is_same && !is_not_same, ExcInternalError());
  }
  {
    constexpr double       initializer_1[3] = {1., 0., 0.};
    constexpr Tensor<1, 3> c_1{initializer_1};
    constexpr double       initializer_2[3] = {0., 1., 0.};
    constexpr Tensor<1, 3> c_2{initializer_2};

    DEAL_II_CONSTEXPR const auto c_cross = cross_product_3d(c_1, c_2);
    constexpr double             initializer_ref[3] = {0., 0., 1.};
    constexpr Tensor<1, 3>       ref{initializer_ref};
    DEAL_II_CONSTEXPR const bool is_same     = (c_cross == ref);
    DEAL_II_CONSTEXPR const bool is_not_same = (c_cross != ref);
    Assert(is_same && !is_not_same, ExcInternalError());
  }

  DEAL_II_CONSTEXPR const auto table_indices =
    Tensor<2, 3>::unrolled_to_component_indices(0);
  DEAL_II_CONSTEXPR const auto index =
    Tensor<2, 3>::component_to_unrolled_index(TableIndices<2>{});
  Assert(index == 0, ExcInternalError());

  DEAL_II_CONSTEXPR const auto used_memory = Tensor<2, 3>::memory_consumption();
  deallog << "Used memory: " << used_memory << std::endl;

  {
    constexpr double a_init[3][3] = {{1., 0., 0.}, {2., 1., 0.}, {3., 2., 1.}};
    constexpr Tensor<2, 3>       a{a_init};
    DEAL_II_CONSTEXPR const auto inverted       = invert(a);
    constexpr double             ref_init[3][3] = {{1., 0., 0.},
                                                   {-2., 1., 0.},
                                                   {1., -2., 1}};
    constexpr Tensor<2, 3>       ref{ref_init};
    Assert(inverted == ref, ExcInternalError());
  }
  {
    constexpr double a_init[3][3] = {{1., 0., 0.}, {2., 1., 0.}, {3., 2., 1.}};
    constexpr Tensor<2, 3>       a{a_init};
    DEAL_II_CONSTEXPR const auto transposed = transpose(a);
    constexpr double ref_init[3][3] = {{1., 2., 3.}, {0., 1., 2.}, {0., 0., 1}};
    constexpr Tensor<2, 3> ref{ref_init};
    Assert(transposed == ref, ExcInternalError());
    DEAL_II_CONSTEXPR const auto dummy   = scalar_product(a, ref);
    DEAL_II_CONSTEXPR const auto dummy_2 = contract<0, 0>(a, ref);
    DEAL_II_CONSTEXPR const auto dummy_3 = double_contract<0, 0, 1, 1>(a, ref);
    DEAL_II_CONSTEXPR const auto dummy_4 = schur_product(a, ref);
    DEAL_II_CONSTEXPR const auto dummy_5 = a * ref;
    DEAL_II_CONSTEXPR const auto middle  = outer_product(a, a);
    DEAL_II_CONSTEXPR const auto dummy_6 = contract3(a, middle, a);
    DEAL_II_CONSTEXPR const auto dummy_7 = adjugate(a);
    DEAL_II_CONSTEXPR const auto dummy_8 = cofactor(a);
  }

  {
    constexpr Tensor<1, 3> dummy_1;
    constexpr Tensor<0, 3> dummy_0;
    DEAL_II_CONSTEXPR auto product_result = dummy_1 * dummy_1;
    DEAL_II_CONSTEXPR auto contraction_result =
      contract<0, 0>(dummy_1, dummy_1);
    DEAL_II_CONSTEXPR auto outer_product_result =
      outer_product(dummy_0, dummy_0);
  }

  deallog << "OK" << std::endl;
  return 0;
}
