// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check methods of the Q1WedgeFunction
// 1-dim does not exist
// 2-dim Q1WedgeFunction = xy
// 3-dim Q1WedgeFunction = xy

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/flow_function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"


template <int dim>
void
check_values(const Function<dim> &f)
{
  Point<dim> p;
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = i + 1;
  deallog << f.value(p) << std::endl;
  deallog << " values checked" << std::endl;
}

template <int dim>
void
check_value_list(const Function<dim> &f)
{
  const unsigned int      max_points = 3;
  std::vector<Point<dim>> point_vector(max_points);
  for (unsigned int j = 0; j < max_points; ++j)
    {
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = j + 1;
      point_vector[j] = p;
    }
  std::vector<double> values(max_points);
  f.value_list(point_vector, values);
  for (unsigned int j = 0; j < max_points; ++j)
    deallog << values[j] << std::endl;
  deallog << " value_list checked" << std::endl;
}

template <int dim>
void
check_vector_value_list(const Function<dim> &f)
{
  const unsigned int      max_points = 3;
  std::vector<Point<dim>> point_vector(max_points);
  for (unsigned int j = 0; j < max_points; ++j)
    {
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = j + 1;
      point_vector[j] = p;
    }
  std::vector<Vector<double>> values(max_points, Vector<double>(1));
  f.vector_value_list(point_vector, values);
  for (unsigned int j = 0; j < max_points; ++j)
    deallog << values[j](0) << std::endl;
  deallog << " vector_value_list checked" << std::endl;
}

template <int dim>
void
check_gradients(const Function<dim> &f)
{
  Point<dim> p;
  for (unsigned int i = 0; i < dim; ++i)
    {
      p[i] = i + 3;
    }
  Tensor<1, dim> grads = f.gradient(p);
  for (unsigned int i = 0; i < dim; ++i)
    {
      deallog << i << "-der: " << grads[i] << std::endl;
    }
  deallog << " gradients checked" << std::endl;
}

template <int dim>
void
check_gradient_list(const Function<dim> &f)
{
  const unsigned int      max_points = 3;
  std::vector<Point<dim>> point_vector(max_points);
  for (unsigned int j = 0; j < max_points; ++j)
    {
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = j + 1;
      point_vector[j] = p;
    }
  std::vector<Tensor<1, dim>> grads(max_points);
  f.gradient_list(point_vector, grads);
  for (unsigned int j = 0; j < max_points; ++j)
    deallog << grads[j] << std::endl;
  deallog << " gradient_list checked" << std::endl;
}

template <int dim>
void
check_vector_gradient_list(const Function<dim> &f)
{
  const unsigned int      max_points = 3;
  std::vector<Point<dim>> point_vector(max_points);
  for (unsigned int j = 0; j < max_points; ++j)
    {
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = j + 1;
      point_vector[j] = p;
    }
  std::vector<std::vector<Tensor<1, dim>>> gradients(
    max_points, std::vector<Tensor<1, dim>>(1));
  f.vector_gradient_list(point_vector, gradients);
  for (unsigned int j = 0; j < max_points; ++j)
    deallog << gradients[j][0] << std::endl;
  deallog << " vector_gradient_list checked" << std::endl;
}

template <int dim>
void
check_laplacian(const Function<dim> &f)
{
  Point<dim> p;
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = i + 1;
  Assert(std::fabs(f.laplacian(p)) < 1e-12, ExcInternalError());
}

template <int dim>
void
check_laplacian_list(const Function<dim> &f)
{
  const unsigned int      max_points = 3;
  std::vector<Point<dim>> point_vector(max_points);
  for (unsigned int j = 0; j < max_points; ++j)
    {
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = j + 1;
      point_vector[j] = p;
    }
  std::vector<double> values(max_points);
  f.laplacian_list(point_vector, values);
  for (unsigned int j = 0; j < max_points; ++j)
    Assert(std::fabs(values[j]) < 1e-12, ExcInternalError());
}

int
main()
{
  initlog();

  if (true)
    {
      deallog << " Functions::Q1WedgeFunction<2>" << std::endl;
      check_values(Functions::Q1WedgeFunction<2>());
      check_values(Functions::Q1WedgeFunction<3>());
      check_value_list(Functions::Q1WedgeFunction<2>());
      check_value_list(Functions::Q1WedgeFunction<3>());
      check_vector_value_list(Functions::Q1WedgeFunction<2>());
      check_vector_value_list(Functions::Q1WedgeFunction<3>());
      check_gradients(Functions::Q1WedgeFunction<2>());
      check_gradients(Functions::Q1WedgeFunction<3>());
      check_gradient_list(Functions::Q1WedgeFunction<2>());
      check_gradient_list(Functions::Q1WedgeFunction<3>());
      check_vector_gradient_list(Functions::Q1WedgeFunction<2>());
      check_vector_gradient_list(Functions::Q1WedgeFunction<3>());
      check_laplacian(Functions::Q1WedgeFunction<2>());
      check_laplacian(Functions::Q1WedgeFunction<3>());
    }
}
