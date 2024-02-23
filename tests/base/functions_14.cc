// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check IdentityFunction

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"



template <int dim>
void
check_value(const Function<dim> &f)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << " d = " << d << std::endl;
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = i;

      deallog << f.value(p, d) << std::endl;
    }
  deallog << " values checked" << std::endl;
}

template <int dim>
void
check_value_list(const Function<dim> &f)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << " d = " << d << std::endl;
      const unsigned int      max_number_of_points = 5;
      std::vector<Point<dim>> points(max_number_of_points);

      for (unsigned int i = 0; i < max_number_of_points; ++i)
        {
          Point<dim> p;
          for (unsigned int j = 0; j < dim; ++j)
            p[j] = i + 1;

          points[i] = p;
        }
      std::vector<double> values(max_number_of_points);
      f.value_list(points, values, d);

      for (unsigned int j = 0; j < max_number_of_points; ++j)
        deallog << values[j] << std::endl;
    }

  deallog << " value_list checked" << std::endl;
}


template <int dim>
void
check_gradient(const Function<dim> &f)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << " d = " << d << std::endl;
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = i;

      deallog << f.gradient(p, d) << std::endl;
    }
  deallog << " gradients checked" << std::endl;
}

template <int dim>
void
check_gradient_list(const Function<dim> &f)
{
  const unsigned int      max_number_of_points = 5;
  std::vector<Point<dim>> points(max_number_of_points);

  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << " d = " << d << std::endl;
      for (unsigned int i = 0; i < max_number_of_points; ++i)
        {
          Point<dim> p;
          for (unsigned int j = 0; j < dim; ++j)
            p[j] = i + 1;

          points[i] = p;
        }
      std::vector<Tensor<1, dim>> tensors(max_number_of_points);
      f.gradient_list(points, tensors, d);

      for (unsigned int j = 0; j < max_number_of_points; ++j)
        deallog << tensors[j] << std::endl;
    }

  deallog << " gradient_list checked" << std::endl;
}


template <int dim>
void
check_laplacian(const Function<dim> &f)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << " d = " << d << std::endl;
      Point<dim> p;
      for (unsigned int i = 0; i < dim; ++i)
        p[i] = i;

      deallog << f.laplacian(p, d) << std::endl;
    }
  deallog << " laplacians checked" << std::endl;
}


template <int dim>
void
check_laplacian_list(const Function<dim> &f)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << " d = " << d << std::endl;
      const unsigned int      max_number_of_points = 5;
      std::vector<Point<dim>> points(max_number_of_points);

      for (unsigned int i = 0; i < max_number_of_points; ++i)
        {
          Point<dim> p;
          for (unsigned int j = 0; j < dim; ++j)
            p[j] = i + 1;

          points[i] = p;
        }
      std::vector<double> values(max_number_of_points);
      f.laplacian_list(points, values, d);

      for (unsigned int j = 0; j < max_number_of_points; ++j)
        deallog << values[j] << std::endl;
    }

  deallog << " laplacian_list checked" << std::endl;
}


int
main()
{
  initlog();

  deallog << "Functions IdentityFunction" << std::endl;
  check_value(Functions::IdentityFunction<1>());
  check_value(Functions::IdentityFunction<2>());
  check_value(Functions::IdentityFunction<3>());

  check_value_list(Functions::IdentityFunction<1>());
  check_value_list(Functions::IdentityFunction<2>());
  check_value_list(Functions::IdentityFunction<3>());

  check_gradient(Functions::IdentityFunction<1>());
  check_gradient(Functions::IdentityFunction<2>());
  check_gradient(Functions::IdentityFunction<3>());

  check_gradient_list(Functions::IdentityFunction<1>());
  check_gradient_list(Functions::IdentityFunction<2>());
  check_gradient_list(Functions::IdentityFunction<3>());

  check_gradient(Functions::IdentityFunction<1>());
  check_gradient(Functions::IdentityFunction<2>());
  check_gradient(Functions::IdentityFunction<3>());

  check_gradient_list(Functions::IdentityFunction<1>());
  check_gradient_list(Functions::IdentityFunction<2>());
  check_gradient_list(Functions::IdentityFunction<3>());

  check_laplacian(Functions::IdentityFunction<1>());
  check_laplacian(Functions::IdentityFunction<2>());
  check_laplacian(Functions::IdentityFunction<3>());

  check_laplacian_list(Functions::IdentityFunction<1>());
  check_laplacian_list(Functions::IdentityFunction<2>());
  check_laplacian_list(Functions::IdentityFunction<3>());
}
