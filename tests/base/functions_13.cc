// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check for support of std::complex<double> as Number template argument

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <complex>

#include "../tests.h"

template <typename Number>
class MyFunctionTime : public dealii::FunctionTime<Number>
{};

template <int dim, typename Number>
class MyFunction : public dealii::Function<dim, Number>
{};

template <int rank, int dim, typename Number>
class MyTensorFunction : public dealii::TensorFunction<rank, dim, Number>
{};

int
main()
{
  initlog();

  MyFunctionTime<std::complex<double>> function_time_1;
  static_assert(std::is_same<decltype(function_time_1.get_time()),
                             std::complex<double>>::value,
                "Wrong return type for FunctionTime<double>::get_time().");
  MyFunctionTime<std::complex<float>> function_time_2;
  static_assert(std::is_same<decltype(function_time_2.get_time()),
                             std::complex<float>>::value,
                "Wrong return type for FunctionTime<float>::get_time().");

  MyFunction<1, std::complex<double>> function_1;
  static_assert(std::is_same_v<decltype(function_1.get_time()), double>,
                "Wrong return type for Function<double>::get_time().");
  MyFunction<2, std::complex<float>> function_2;
  static_assert(std::is_same_v<decltype(function_2.get_time()), float>,
                "Wrong return type for Function<float>::get_time().");

  MyTensorFunction<1, 1, std::complex<double>> tensor_function_1;
  static_assert(std::is_same_v<decltype(tensor_function_1.get_time()), double>,
                "Wrong return type for TensorFunction<double>::get_time().");
  MyTensorFunction<2, 2, std::complex<float>> tensor_function_2;
  static_assert(std::is_same_v<decltype(tensor_function_2.get_time()), float>,
                "Wrong return type for TensorFunction<float>::get_time().");

  deallog << "OK" << std::endl;
}
