// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the general number capability of symmetric and normal tensors

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <type_traits>

#include "../tests.h"

template <int rank, int dim, typename NumberType>
void
fill_tensor(SymmetricTensor<rank, dim, NumberType> &t)
{
  for (unsigned int i = 0; i != t.n_independent_components; ++i)
    {
      t.access_raw_entry(i) = i + 1;
    }
}


template <int rank, int dim, typename NumberType>
void
fill_tensor(Tensor<rank, dim, NumberType> &t)
{
  for (unsigned int i = 0; i != t.n_independent_components; ++i)
    {
      t[t.unrolled_to_component_indices(i)] = i + 1;
    }
}

template <int dim,
          template <int, int, typename>
          class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
void
print(const TensorType1<2, dim, NumberType1> &t2,
      const TensorType2<4, dim, NumberType2> &t4)
{
  deallog << t2 << std::endl;
  deallog << t4 << std::endl;
}

template <int dim,
          template <int, int, typename>
          class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
void
print(const TensorType1<2, dim, NumberType1> &t2_1,
      const TensorType2<2, dim, NumberType2> &t2_2,
      const TensorType1<4, dim, NumberType1> &t4_1,
      const TensorType2<4, dim, NumberType2> &t4_2)
{
  deallog << t2_1 << std::endl;
  deallog << t2_2 << std::endl;
  deallog << t4_1 << std::endl;
  deallog << t4_2 << std::endl;
}


template <template <int, int, typename> class TensorType1,
          template <int, int, typename>
          class TensorType2>
struct AreSame : std::false_type
{};

template <template <int, int, typename> class TensorType1>
struct AreSame<TensorType1, TensorType1> : std::true_type
{};


template <template <int, int, typename> class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
void
test_one()
{
  const unsigned int               dim = 2;
  TensorType1<2, dim, NumberType1> t2_1;
  TensorType2<2, dim, NumberType2> t2_2;
  TensorType1<4, dim, NumberType1> t4_1;
  TensorType2<4, dim, NumberType2> t4_2;

  deallog << "Fill" << std::endl;
  fill_tensor(t2_1);
  fill_tensor(t2_2);
  fill_tensor(t4_1);
  fill_tensor(t4_2);
  print(t2_1, t2_2, t4_1, t4_2);

  deallog << "Product / Division" << std::endl;
  t2_1 *= 2.0;
  t2_2 /= 2.0;
  t4_1 *= 2.0;
  t4_2 /= 2.0;
  print(t2_1, t2_2, t4_1, t4_2);

  deallog << "Operator * (single / double contraction)" << std::endl;
  deallog << (t2_1 * t2_2) << std::endl;
  deallog << (t4_1 * t4_2) << std::endl;

  deallog << "Scalar product" << std::endl;
  deallog << scalar_product(t2_1, t2_2) << std::endl;
}

template <template <int, int, typename> class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
std::enable_if_t<AreSame<TensorType1, TensorType2>::value>
test_two()
{
  const unsigned int               dim = 2;
  TensorType1<2, dim, NumberType1> t2_1;
  TensorType1<2, dim, NumberType2> t2_2;
  TensorType1<4, dim, NumberType1> t4_1;
  TensorType1<4, dim, NumberType2> t4_2;

  deallog << "Fill" << std::endl;
  fill_tensor(t2_1);
  fill_tensor(t4_1);
  print(t2_1, t2_2, t4_1, t4_2);

  deallog << "Equals" << std::endl;
  t2_2 = t2_1;
  t4_2 = t4_1;
  deallog << t2_2 << std::endl;
  deallog << t4_2 << std::endl;

  deallog << "Zero" << std::endl;
  t2_1 = 0;
  t4_1 = 0;
  deallog << t2_1 << std::endl;
  deallog << t4_1 << std::endl;

  deallog << "Cast" << std::endl;
  t2_2 = static_cast<decltype(t2_2)>(t2_1);
  t4_2 = static_cast<decltype(t4_2)>(t4_1);
  deallog << t2_2 << std::endl;
  deallog << t4_2 << std::endl;

  deallog << "In place addition" << std::endl;
  t2_1 += t2_2;
  t4_1 += t4_2;
  deallog << t2_1 << std::endl;
  deallog << t4_1 << std::endl;
  print(t2_1, t4_1);

  deallog << "In place subtraction" << std::endl;
  t2_1 -= t2_2;
  t4_1 -= t4_2;
  print(t2_1, t4_1);

  deallog << "Operator * (double contraction)" << std::endl;
  deallog << (t2_1 * t2_2) << std::endl;
  deallog << (t4_1 * t4_2) << std::endl;
}

template <template <int, int, typename> class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
std::enable_if_t<!AreSame<TensorType1, TensorType2>::value>
test_two()
{}

template <template <int, int, typename> class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
std::enable_if_t<(AreSame<TensorType1, SymmetricTensor>::value &&
                  AreSame<TensorType2, SymmetricTensor>::value)>
test_three()
{
  const unsigned int               dim = 2;
  TensorType1<2, dim, NumberType1> t2_1;
  TensorType1<4, dim, NumberType2> t4_2;

  fill_tensor(t2_1);
  fill_tensor(t4_2);

  deallog << "Double contract" << std::endl;
  TensorType1<2, dim, typename ProductType<NumberType1, NumberType2>::type>
    res_dc;
  double_contract(res_dc, t4_2, t2_1);
  deallog << res_dc << std::endl;
  double_contract(res_dc, t2_1, t4_2);
  deallog << res_dc << std::endl;
}

template <template <int, int, typename> class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
std::enable_if_t<!(AreSame<TensorType1, SymmetricTensor>::value &&
                   AreSame<TensorType2, SymmetricTensor>::value)>
test_three()
{}


template <template <int, int, typename> class TensorType1,
          typename NumberType1,
          template <int, int, typename>
          class TensorType2,
          typename NumberType2>
void
test_all()
{
  test_one<TensorType1, NumberType1, TensorType2, NumberType2>();
  test_two<TensorType1, NumberType1, TensorType2, NumberType2>();
  test_three<TensorType1, NumberType1, TensorType2, NumberType2>();
}

template <typename Number1, typename Number2>
void
test_T()
{
  deallog.push("ST,ST");
  test_all<SymmetricTensor, double, SymmetricTensor, double>();
  deallog.pop();

  deallog.push("T,T");
  test_all<Tensor, double, Tensor, double>();
  deallog.pop();

  deallog.push("ST,T");
  test_all<SymmetricTensor, double, Tensor, double>();
  deallog.pop();

  deallog.push("T,ST");
  test_all<Tensor, double, SymmetricTensor, double>();
  deallog.pop();
}

int
main()
{
  initlog();

  // Benchmark results for same types
  deallog.push("d,d");
  test_T<double, double>();
  deallog.pop();

  // Differing number types:
  deallog.push("d,f");
  test_T<double, float>();
  deallog.pop();
  deallog.push("f,d");
  test_T<float, double>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
