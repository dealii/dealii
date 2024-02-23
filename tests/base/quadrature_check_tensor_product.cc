// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the boolean is_tensor_product for all the quadrature classes


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

template <int dim>
void
check_tensor_product(const std::vector<Quadrature<dim>> &quadratures,
                     const std::vector<std::string>     &quadrature_names)
{
  DEAL_II_NOT_IMPLEMENTED();
}

template <>
void
check_tensor_product(const std::vector<Quadrature<1>> &quadratures,
                     const std::vector<std::string>   &quadrature_names)
{
  for (unsigned int i = 0; i < quadratures.size(); ++i)
    {
      const Quadrature<1> &quadrature = quadratures[i];
      if (quadrature.is_tensor_product())
        {
          deallog << "1D " << quadrature_names[i];
          const auto &q_basis = quadrature.get_tensor_basis();
          AssertThrow(q_basis.size() == 1, ExcInternalError());
          const auto &q_points  = quadrature.get_points();
          const auto &q_weights = quadrature.get_weights();
          AssertThrow(q_basis[0].size() == q_points.size(), ExcInternalError());
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              std::cout << q_points[q] << ' ' << q_basis[0].get_points()[q]
                        << std::endl;
              AssertThrow(std::abs(
                            (q_points[q] - q_basis[0].get_points()[q]).norm()) <
                            1.e-10,
                          ExcInternalError());
              AssertThrow(std::abs(q_weights[q] - q_basis[0].get_weights()[q]) <
                            1.e-10,
                          ExcInternalError());
            }
          deallog << " OK" << std::endl;
        }
    }
}

template <>
void
check_tensor_product(const std::vector<Quadrature<2>> &quadratures,
                     const std::vector<std::string>   &quadrature_names)
{
  for (unsigned int i = 0; i < quadratures.size(); ++i)
    {
      const Quadrature<2> &quadrature = quadratures[i];
      if (quadrature.is_tensor_product())
        {
          deallog << "2D " << quadrature_names[i];
          const auto &q_basis = quadrature.get_tensor_basis();
          AssertThrow(q_basis.size() == 2, ExcInternalError());
          const auto &q_points  = quadrature.get_points();
          const auto &q_weights = quadrature.get_weights();
          AssertThrow(q_basis[0].size() * q_basis[1].size() == q_points.size(),
                      ExcInternalError());
          unsigned int q = 0;
          for (unsigned int q2 = 0; q2 < q_basis[0].size(); ++q2)
            for (unsigned int q1 = 0; q1 < q_basis[1].size(); ++q1)
              {
                AssertThrow(std::abs(q_points[q][0] -
                                     q_basis[0].get_points()[q1][0]) < 1.e-10,
                            ExcInternalError());
                AssertThrow(std::abs(q_points[q][1] -
                                     q_basis[1].get_points()[q2][0]) < 1.e-10,
                            ExcInternalError());
                AssertThrow(std::abs((q_weights[q] -
                                      q_basis[0].get_weights()[q1] *
                                        q_basis[1].get_weights()[q2])) < 1.e-10,
                            ExcInternalError());
                ++q;
              }
          deallog << " OK" << std::endl;
        }
    }
}

template <>
void
check_tensor_product(const std::vector<Quadrature<3>> &quadratures,
                     const std::vector<std::string>   &quadrature_names)
{
  for (unsigned int i = 0; i < quadratures.size(); ++i)
    {
      const Quadrature<3> &quadrature = quadratures[i];
      if (quadrature.is_tensor_product())
        {
          deallog << "3D " << quadrature_names[i];
          const auto &q_basis = quadrature.get_tensor_basis();
          AssertThrow(q_basis.size() == 3, ExcInternalError());
          const auto &q_points  = quadrature.get_points();
          const auto &q_weights = quadrature.get_weights();
          AssertThrow(q_basis[0].size() * q_basis[1].size() *
                          q_basis[2].size() ==
                        q_points.size(),
                      ExcInternalError());
          unsigned int q = 0;
          for (unsigned int q3 = 0; q3 < q_basis[2].size(); ++q3)
            for (unsigned int q2 = 0; q2 < q_basis[1].size(); ++q2)
              for (unsigned int q1 = 0; q1 < q_basis[0].size(); ++q1)
                {
                  AssertThrow(std::abs(q_points[q][0] -
                                       q_basis[0].get_points()[q1][0]) < 1.e-10,
                              ExcInternalError());
                  AssertThrow(std::abs(q_points[q][1] -
                                       q_basis[1].get_points()[q2][0]) < 1.e-10,
                              ExcInternalError());
                  AssertThrow(std::abs(q_points[q][2] -
                                       q_basis[2].get_points()[q3][0]) < 1.e-10,
                              ExcInternalError());
                  AssertThrow(std::abs(q_weights[q] -
                                       q_basis[0].get_weights()[q1] *
                                         q_basis[1].get_weights()[q2] *
                                         q_basis[2].get_weights()[q3]) < 1.e-10,
                              ExcInternalError());
                  ++q;
                }
          deallog << " OK" << std::endl;
        }
    }
}

template <int dim>
void
fill_quadrature_vector(std::vector<Quadrature<dim>> &quadratures,
                       std::vector<std::string>     &quadrature_names)
{
  quadratures.push_back(Quadrature<dim>());
  quadrature_names.push_back("Quadrature");

  quadratures.push_back(QIterated<dim>(QGauss<1>(2), 2));
  quadrature_names.push_back("QIterated");

  quadratures.push_back(QGauss<dim>(2));
  quadrature_names.push_back("QGauss");

  quadratures.push_back(QGaussLobatto<dim>(2));
  quadrature_names.push_back("QGaussLobatto");

  quadratures.push_back(QMidpoint<dim>());
  quadrature_names.push_back("QMidPoint");

  quadratures.push_back(QSimpson<dim>());
  quadrature_names.push_back("QSimpson");

  quadratures.push_back(QTrapezoid<dim>());
  quadrature_names.push_back("QTrapezoid");

  quadratures.push_back(QMilne<dim>());
  quadrature_names.push_back("QMilne");

  quadratures.push_back(QWeddle<dim>());
  quadrature_names.push_back("QWeddle");

  quadratures.push_back(QGaussChebyshev<dim>(3));
  quadrature_names.push_back("QGaussChebyshev");

  quadratures.push_back(QGaussRadauChebyshev<dim>(2));
  quadrature_names.push_back("QGaussRadauChebyshev");

  quadratures.push_back(QGaussLobattoChebyshev<dim>(2));
  quadrature_names.push_back("QGaussLobattoChebyshev");

  quadratures.push_back(QSorted<dim>(Quadrature<dim>()));
  quadrature_names.push_back("QSorted");

  quadratures.push_back(QTelles<dim>(1, Point<dim>()));
  quadrature_names.push_back("QTelles");
}

int
main()
{
  initlog();
  deallog << std::boolalpha;

  Quadrature<1> q;

  std::vector<Quadrature<1>> quadratures_1d;
  std::vector<std::string>   quadrature_names_1d;
  fill_quadrature_vector(quadratures_1d, quadrature_names_1d);
  quadratures_1d.push_back(QAnisotropic<1>(q));
  quadrature_names_1d.push_back("QAnisotropic");
  quadratures_1d.push_back(QGaussLog<1>(1));
  quadrature_names_1d.push_back("QGaussLog");
  quadratures_1d.push_back(QGaussLogR<1>(1));
  quadrature_names_1d.push_back("QGaussLogR");
  check_tensor_product(quadratures_1d, quadrature_names_1d);

  std::vector<Quadrature<2>> quadratures_2d;
  std::vector<std::string>   quadrature_names_2d;
  fill_quadrature_vector(quadratures_2d, quadrature_names_2d);
  quadratures_2d.push_back(QAnisotropic<2>(q, q));
  quadrature_names_2d.push_back("QAnisotropic");
  quadratures_2d.push_back(QGaussOneOverR<2>(1, Point<2>()));
  quadrature_names_2d.push_back("QGaussOneOverR");
  check_tensor_product(quadratures_2d, quadrature_names_2d);

  std::vector<Quadrature<3>> quadratures_3d;
  std::vector<std::string>   quadrature_names_3d;
  fill_quadrature_vector(quadratures_3d, quadrature_names_3d);
  quadratures_3d.push_back(QAnisotropic<3>(q, q, q));
  quadrature_names_3d.push_back("QAnisotropic");
  check_tensor_product(quadratures_3d, quadrature_names_3d);
}
