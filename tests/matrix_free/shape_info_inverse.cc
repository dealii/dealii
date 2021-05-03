// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// test the correctness of the inverse_shape_values field of
// internal::MatrixFreeFunctions::ShapeInfo

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/matrix_free/shape_info.h>

#include <iostream>

#include "../tests.h"



template <int dim>
void
test(const FiniteElement<dim> &fe,
     const Quadrature<1> &     quad,
     const std::string &       quadrature_name)
{
  internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;
  shape_info.reinit(quad, fe, 0);
  deallog << "Testing " << fe.get_name() << " with " << quadrature_name << "("
          << quad.size() << ")" << std::endl;
  deallog << "shape values: " << std::endl;
  const auto &univariate_shape_data = shape_info.get_shape_data(0, 0);
  for (unsigned int i = 0; i < univariate_shape_data.fe_degree + 1; ++i)
    {
      for (unsigned int q = 0; q < quad.size(); ++q)
        deallog << std::setw(15)
                << univariate_shape_data.shape_values[i * quad.size() + q]
                << " ";
      deallog << std::endl;
    }
  deallog << "inverse shape values: " << std::endl;
  for (unsigned int i = 0; i < univariate_shape_data.fe_degree + 1; ++i)
    {
      for (unsigned int q = 0; q < quad.size(); ++q)
        deallog
          << std::setw(15)
          << univariate_shape_data.inverse_shape_values[i * quad.size() + q]
          << " ";
      deallog << std::endl;
    }
  deallog << "inverse shapes' * shapes: " << std::endl;
  for (unsigned int i = 0; i < quad.size(); ++i)
    {
      for (unsigned int j = 0; j < quad.size(); ++j)
        {
          double sum = 0;
          for (unsigned int k = 0; k < univariate_shape_data.fe_degree + 1; ++k)
            sum +=
              univariate_shape_data.inverse_shape_values[k * quad.size() + i] *
              univariate_shape_data.shape_values[k * quad.size() + j];
          deallog << std::setw(15) << sum << " ";
        }
      deallog << std::endl;
    }
  deallog << "inverse shapes * shapes': " << std::endl;
  for (unsigned int i = 0; i < univariate_shape_data.fe_degree + 1; ++i)
    {
      for (unsigned int j = 0; j < univariate_shape_data.fe_degree + 1; ++j)
        {
          double sum = 0;
          for (unsigned int k = 0; k < quad.size(); ++k)
            sum +=
              univariate_shape_data.inverse_shape_values[i * quad.size() + k] *
              univariate_shape_data.shape_values[j * quad.size() + k];
          deallog << std::setw(15) << sum << " ";
        }
      deallog << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(9);
  test(FE_DGQ<1>(1), QGauss<1>(2), "QGauss");
  test(FE_DGQ<1>(1), QGauss<1>(3), "QGauss");
  test(FE_DGQ<1>(3), QGauss<1>(2), "QGauss");
  test(FE_DGQ<1>(3), QGauss<1>(4), "QGauss");
  test(FE_DGQ<1>(3), QGauss<1>(5), "QGauss");
  test(FE_DGQ<1>(8), QGauss<1>(7), "QGauss");
  test(FE_DGQ<1>(8), QGauss<1>(9), "QGauss");
  test(FE_DGQ<1>(8), QGauss<1>(12), "QGauss");
  test(FE_DGQHermite<1>(1), QGauss<1>(2), "QGauss");
  test(FE_DGQHermite<1>(1), QGauss<1>(3), "QGauss");
  test(FE_DGQHermite<1>(3), QGauss<1>(2), "QGauss");
  test(FE_DGQHermite<1>(3), QGauss<1>(4), "QGauss");
  test(FE_DGQHermite<1>(3), QGauss<1>(5), "QGauss");
  test(FE_DGQHermite<1>(8), QGauss<1>(7), "QGauss");
  test(FE_DGQHermite<1>(8), QGauss<1>(9), "QGauss");
  test(FE_DGQHermite<1>(8), QGauss<1>(12), "QGauss");
  test(FE_DGQArbitraryNodes<1>(QGauss<1>(9)), QGauss<1>(9), "QGauss");
  test(FE_DGQArbitraryNodes<1>(QGauss<1>(9)), QGauss<1>(12), "QGauss");
  test(FE_DGQ<1>(3), QGaussLobatto<1>(4), "QGaussLobatto");
  test(FE_DGQ<1>(3), QGaussLobatto<1>(5), "QGaussLobatto");
  test(FE_Q<1>(3), QGauss<1>(4), "QGauss");
  test(FE_Q<1>(3), QGauss<1>(5), "QGauss");

  test(FE_DGQ<2>(3), QGauss<1>(4), "QGauss");
  test(FE_DGQ<2>(3), QGauss<1>(5), "QGauss");
  test(FE_Q<2>(3), QGauss<1>(4), "QGauss");
  test(FE_Q<2>(3), QGauss<1>(5), "QGauss");
}
