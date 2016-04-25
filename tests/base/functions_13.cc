// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
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


// Check for suport of std::complex<double> as Number template argument

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <complex>

using namespace dealii;

template <typename Number>
class MyFunctionTime : dealii::FunctionTime<Number>
{
};

template <int dim, typename Number>
class MyFunction : dealii::Function<dim, Number>
{
};

template <int rank, int dim, typename Number>
class MyTensorFunction : dealii::TensorFunction<rank, dim, Number>
{
};

int main()
{
  initlog();

  MyFunction<1, std::complex<double> > function_1;
  MyFunction<2, std::complex<float> > function_2;

  MyTensorFunction<1, 1, std::complex<double> > tensor_function_1;
  MyTensorFunction<2, 2, std::complex<float> > tensor_function_2;

  deallog << "OK" << std::endl;
}
